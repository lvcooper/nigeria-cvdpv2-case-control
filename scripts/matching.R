
matching<-function(define.distance,define.age,define.onset,match.on.RI=F,match.limit=NA,under=NA){
  d<-read.csv("data/eligible.csv")
  
  d$date<-as.Date(d$date)
  if(match.on.RI) d<-d[!is.na(d$doses_opv_ri)&!is.na(d$doses_ipv012),]
  if(!is.na(under)) d<-d[d$age<under,]
  cases<-d[d$case,]
  controls<-d[!d$case,]
  
  distFlat<-function(mat,vec){
    mat[,1]<-(mat[,1]-vec[1])^2
    mat[,2]<-(mat[,2]-vec[2])^2
    sqrt(rowSums(mat))
  }
  
  d$date<-as.Date(d$date)
  controls$keep<-F
  if(define.distance>1){
  for(i in 1:nrow(cases)){
    controls$keep[abs(controls$age-cases$age[i])<=define.age&
                    abs(as.numeric(controls$date-cases$date[i]))<=define.onset&
                    controls$adm1==cases$adm1[i]&
                    distHaversine(controls[,c("CENTER_LON","CENTER_LAT")],c(cases$CENTER_LON[i],cases$CENTER_LAT[i]))<=define.distance]<-T
  }
  controls<-controls[controls$keep,]
  }
  if(define.distance>=1){
  for(i in 1:nrow(cases)){
    controls$keep[abs(controls$age-cases$age[i])<=define.age&
                    abs(as.numeric(controls$date-cases$date[i]))<=define.onset&
                    controls$adm1==cases$adm1[i]&
                    controls$adm2==cases$adm2[i]]<-T
  }
  controls<-controls[controls$keep,]
  }
  # Closest method
  x<-cases[,c("date","age","doses_opv_ri","doses_ipv","CENTER_LON","CENTER_LAT","GUID","EPID","adm1","adm2")]
  y<-controls[,c("date","age","doses_opv_ri","doses_ipv","CENTER_LON","CENTER_LAT","GUID","EPID","adm1","adm2")]
  names(x)<-paste0("case_",c("date","age","ri","ipv","lon","lat","GUID","EPID","adm1","adm2"))
  names(y)<-paste0("control_",c("date","age","ri","ipv","lon","lat","GUID","EPID","adm1","adm2"))
  z<-merge(x,y)
  z$date_dist<-as.numeric(abs(z$case_date-z$control_date))/define.onset
  z$age_dist<-abs(z$case_age-z$control_age)/define.age
  # All possible matches
  if(match.on.RI) {
    # If matching on RI, case and control must have same IPV history (0 versus 1 or more doses) and same OPV RI history (0-2 doses versus 3+ doses)
    z<-z[as.numeric(z$case_ri>=3)==as.numeric(z$control_ri>=3)&as.numeric(z$case_ipv>0)==as.numeric(z$control_ipv>0),]}
  if(define.distance>1){
  k<-unique(z[,c("case_lon","case_lat","control_lon","control_lat")])
  k$space_dist<-ifelse(k$case_lon==k$control_lon&k$case_lat==k$control_lat,0,NA)
  k$space_dist[which(is.na(k$space_dist))]<-sapply(which(is.na(k$space_dist)),function(i){distHaversine(c(k$case_lon[i],k$case_lat[i]),c(k$control_lon[i],k$control_lat[i]))})
  z<-merge(z,k)
  z$space_dist<-z$space_dist/define.distance
  z$space_dist[which(z$case_adm1!=z$control_adm1)]<-2 #adding to constrain within same adm1 if relax spatial distance
  }
  if(define.distance<=1){z$space_dist<-ifelse(z$case_adm1==z$control_adm1&z$case_adm2==z$control_adm2,0,2)}
  z<-z[z$date_dist<=1&z$age_dist<=1&z$space_dist<=1,]
  # Aggregate distance in age and onset
  z$dist<-z$age_dist+z$date_dist+z$space_dist
  round<-0
  r<-NULL

  while(nrow(z)>0){
    round<-round+1
    # print(round)
    # For each case, select the closest control
    a<-mutate(group_by(as_tibble(z),case_EPID),min_dist=min(dist))
    a<-a[a$dist==a$min_dist,]
    # If there are two controls that are the same distance, choose one at random
    if(!all(table(a$case_EPID)==1)){
      a<-mutate(group_by(as_tibble(a),case_EPID),N=n(),n=sample(1:n()))
      a<-a[a$n==1,]}
    # If two cases have been assigned the same control, choose the closer case
    if(!all(table(a$control_EPID)==1)){
      a<-mutate(group_by(as_tibble(a),control_EPID),min_dist=min(dist))
      a<-a[a$dist==a$min_dist,]
      # If there are two cases that are the same distance, choose one at random
      if(!all(table(a$control_EPID)==1)){
        a<-mutate(group_by(as_tibble(a),control_EPID),N=n(),n=sample(1:n()))
        a<-a[a$n==1,]
      }
    }
    # Store result
    r[[round]]<-a
    # Repeat with remaining controls
    z<-z[!z$control_EPID%in%a$control_EPID,]
  }

  res1<-do.call("rbind",r)
  if(!is.na(match.limit)){
    res1<-res1[order(res1$dist,decreasing=F),]
    res1<-mutate(group_by(as_tibble(res1),case_EPID),rank=1:n())
    res1<-res1[res1$rank<=match.limit,]
  }
  # Closest match
  cases1<-cases[cases$EPID%in%res1$case_EPID,]
  controls1<-controls[controls$EPID%in%res1$control_EPID,names(cases)]
  key<-res1$case_EPID
  names(key)<-res1$control_EPID
  cases1$stratum<-cases1$EPID
  controls1$stratum<-key[controls1$EPID]
  m1<-rbind(cases1,controls1)
  m1$case_epid<-m1$stratum
  m1$stratum<-as.numeric(as.factor(m1$stratum))
  
  
  write.csv(m1,paste0("data/matched_distance_",define.distance,
                      "_age_",define.age,
                      "_onset_",define.onset,
                      ifelse(match.on.RI,"_RI",""),
                      ifelse(!is.na(match.limit),paste0("_limit_",match.limit),""),
                      ifelse(!is.na(under),paste0("_under_",under),""),".csv"),row.names=F)
  
  flow<-read.csv("data/flow.csv")
  flow<-flow[flow$X!="doses_ipv",]
  flow<-rbind(flow,cbind.data.frame(X="matched",Control=sum(m1$case==F),Case=sum(m1$case)))
  m1<-m1[!is.na(m1$doses_ipv),]
  m1<-m1[m1$stratum%in%m1$stratum[m1$case==1],]
  m1<-m1[m1$stratum%in%m1$stratum[m1$case==0],]
  flow<-rbind(flow,cbind.data.frame(X="matched_ipv",Control=sum(m1$case==F),Case=sum(m1$case)))
  write.csv(flow,file=paste0("data/flow_match_distance_",define.distance,
                             "_age_",define.age,
                             "_onset_",define.onset,
                             ifelse(match.on.RI,"_RI",""),
                             ifelse(!is.na(match.limit),paste0("_limit_",match.limit),""),
                             ifelse(!is.na(under),paste0("_under_",under),""),".csv"),row.names=F)
  
  
  
  
}
