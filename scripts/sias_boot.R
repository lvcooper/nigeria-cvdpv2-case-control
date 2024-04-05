sias_boot<-function(define.distance,define.age,define.onset,match.on.RI=F,define.offset,use.date.of.last.opv=T,boot){
  m1<-read.csv(paste0("data/matched_distance_",define.distance,
                      "_age_",define.age,
                      "_onset_",define.onset,
                      ifelse(match.on.RI,"_RI",""),".csv"))
  # Add SIA data
  sia<-read.csv("data/SIAs_Nigeria_case_control.csv")
  sia$date<-as.Date(sia$date)
  
  if(use.date.of.last.opv==F) m1$last_opv_date<-NA
  for(i in c("dob","date","investigation_date","last_opv_date")){m1[,i]<-as.Date(m1[,i])}
  
  # New sampler
  library(arrangements)
  
  possible_histories<-function(v,s){
    v_received=combinations(1:length(s),v,freq=s)
    x=apply(v_received,1,function(x){tabulate(x,nbins=length(s))})
    q=sapply(1:ncol(x),function(k){
      d=x[,k]
      (factorial(s[1])/factorial(s[1]-d[1]))*
        (factorial(s[2])/factorial(s[2]-d[2]))*
        (factorial(s[3])/factorial(s[3]-d[3]))*
        (factorial(s[4])/factorial(s[4]-d[4]))*
        (factorial(sum(s)-v)/factorial(sum(s)))*
        (factorial(v)/(factorial(d[1])*factorial(d[2])*factorial(d[3])*factorial(d[4])))
    })
    r<-cbind.data.frame(t(x),q)
    names(r)<-c(names(s),"prob")
    idx<-sample(nrow(r),boot,replace=T,prob=r$prob)
    return(cbind(r[idx,"nOPV2"],r[idx,"mOPV2"]))
  }
  
  sampler<-function(doses,sias,before){
    # doses=3
    # sias=c("mOPV2","bOPV","bOPV","nOPV2")
    # before=c(T,T,T,F)
    sias_before=sias
    sias_before[!before]=""
    m<-matrix(data=0,nrow=boot,ncol=4)
    m[,1:2]<-possible_histories(v=doses,s=table(factor(sias,levels=c("","bOPV","mOPV2","nOPV2"))))
    m[,3:4]<-possible_histories(v=doses,s=table(factor(sias_before,levels=c("","bOPV","mOPV2","nOPV2"))))
    return(m)
  }
  
  
  sample_doses<-function(EPID,guid,dob,onset_date,investigation_date,date_last_opv=NA,doses,offset){
    # EPID=m1$EPID[i]; dob=m1$dob[i];onset_date=m1$date[i];investigation_date=m1$investigation_date[i];guid=m1$GUID[i];date_last_opv=m1$last_opv_date[i];doses=m1$doses_opv_sia[i];offset=28
    guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
    sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
    # if give date of last dose but report 0 doses, give 1 SIA dose
    if(!is.na(date_last_opv)&doses==0) doses<-1 # this represents 1 case/control
    # All eligible SIAs between birth and investigation
    sub<-sia[which(sia$date>dob&
                     sia$date<investigation_date&
                     sia$admin2guid%in%guid& 
                     sia$age<=sia$age_max&
                     sia$age>=sia$age_min),]
    sub<-unique(sub[,c("date","vaccine_type")])
    OPV_SIAs<-sum(grepl("OPV",sub$vaccine_type))
    IPV_SIAs<-sum(grepl("IPV",sub$vaccine_type))
    sub$vaccine_type[sub$vaccine_type%in%c("IPV + bOPV")]<-"bOPV"
    sub<-sub[!sub$vaccine_type%in%c("IPV","f-IPV"),]
    exp<-sub
    exp_ons<-exp[exp$date<=(onset_date-offset),]
    # If zero SIAs or all SIAs are bOPV, return
    if(nrow(sub)==0|all(sub$vaccine_type=="bOPV")){
      return(cbind.data.frame(nOPV2=0,
                              mOPV2=0,
                              nOPV2_ons=0,
                              mOPV2_ons=0,
                              mOPV2_SIAs=0,
                              nOPV2_SIAs=0,
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=0,
                              nOPV2_SIAs_ons=0,
                              OPV_SIAs=OPV_SIAs,
                              IPV_SIAs=IPV_SIAs,
                              date_used=F,type="",calendar_date=NA,boot=1:boot,EPID=EPID))
    }
    # If date of last OPV consistent with any SIA, assign as "known"
    if(any(abs(sub$date-date_last_opv)<offset)&!is.na(date_last_opv)){
      known<-sub[which.min(abs(sub$date-date_last_opv)),]
      doses_remaining<-doses-1
      date_used<-T
      type<-known$vaccine_type
      calendar_date<-known$date} else{
        known<-NULL
        doses_remaining<-doses
        date_used<-F
        type<-""
        calendar_date<-NA
      }
    # If zero doses remaining, return
    if(doses_remaining==0){
      return(cbind.data.frame(nOPV2=sum(known$vaccine_type=="nOPV2"),
                              mOPV2=sum(known$vaccine_type=="mOPV2"),
                              nOPV2_ons=sum(known$vaccine_type=="nOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_ons=sum(known$vaccine_type=="mOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_SIAs=sum(exp$vaccine_type=="mOPV2"),
                              nOPV2_SIAs=sum(exp$vaccine_type=="nOPV2"),
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="mOPV2"),
                              nOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="nOPV2"),
                              OPV_SIAs=OPV_SIAs,
                              IPV_SIAs=IPV_SIAs,
                              date_used=date_used,type=type,calendar_date=calendar_date,boot=1:boot,EPID=EPID))
    }
    # Remove any SIAs falling after the date of last reported SIA
    # If date of last reported SIA not close to any known SIA, ignore
    sias_remaining<-sub
    if(!is.null(known)){sub<-sub[sub$date<=known$date,]
    sias_remaining<-sub[sub$date<known$date,]}
    # If exactly same number of doses as SIAs
    if(doses_remaining==nrow(sias_remaining)){
      return(cbind.data.frame(nOPV2=sum(sub$vaccine_type=="nOPV2"),
                              mOPV2=sum(sub$vaccine_type=="mOPV2"),
                              nOPV2_ons=sum(sub$vaccine_type=="nOPV2"&sub$date<=(onset_date-offset)),
                              mOPV2_ons=sum(sub$vaccine_type=="mOPV2"&sub$date<=(onset_date-offset)),
                              mOPV2_SIAs=sum(exp$vaccine_type=="mOPV2"),
                              nOPV2_SIAs=sum(exp$vaccine_type=="nOPV2"),
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="mOPV2"),
                              nOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="nOPV2"),
                              OPV_SIAs=OPV_SIAs,
                              IPV_SIAs=IPV_SIAs,
                              date_used=date_used,type=type,calendar_date=calendar_date,boot=1:boot,EPID=EPID))
    }
    # Otherwise sample remaining doses
    if(doses_remaining<nrow(sias_remaining)){
      m<-sampler(doses=doses_remaining,sias=sias_remaining$vaccine_type,before=sias_remaining$date<=(onset_date-offset))
      return(cbind.data.frame(nOPV2=m[,1]+sum(known$vaccine_type=="nOPV2"),
                              mOPV2=m[,2]+sum(known$vaccine_type=="mOPV2"),
                              nOPV2_ons=m[,3]+sum(known$vaccine_type=="nOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_ons=m[,4]+sum(known$vaccine_type=="mOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_SIAs=sum(exp$vaccine_type=="mOPV2"),
                              nOPV2_SIAs=sum(exp$vaccine_type=="nOPV2"),
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="mOPV2"),
                              nOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="nOPV2"),
                              OPV_SIAs=OPV_SIAs,
                              IPV_SIAs=IPV_SIAs,
                              date_used=date_used,type=type,calendar_date=calendar_date,boot=1:boot,EPID=EPID))
    }
    # If over-reporting
    if(doses_remaining>nrow(sias_remaining)){
      doses_remaining<-doses_remaining-nrow(sias_remaining)
      known<-rbind(known,sias_remaining)
      while(doses_remaining>nrow(sub)){
        doses_remaining=doses_remaining-nrow(sub)
        known<-rbind(known,sub)
      }
      m<-sampler(doses=doses_remaining,sias=sub$vaccine_type,before=sub$date<=(onset_date-offset))
      return(cbind.data.frame(nOPV2=m[,1]+sum(known$vaccine_type=="nOPV2"),
                              mOPV2=m[,2]+sum(known$vaccine_type=="mOPV2"),
                              nOPV2_ons=m[,3]+sum(known$vaccine_type=="nOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_ons=m[,4]+sum(known$vaccine_type=="mOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_SIAs=sum(exp$vaccine_type=="mOPV2"),
                              nOPV2_SIAs=sum(exp$vaccine_type=="nOPV2"),
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="mOPV2"),
                              nOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="nOPV2"),
                              OPV_SIAs=OPV_SIAs,
                              IPV_SIAs=IPV_SIAs,
                              date_used=date_used,type=type,calendar_date=calendar_date,boot=1:boot,EPID=EPID))
    }
  }
  
  r<-NULL
  for(i in 1:nrow(m1)){
    r[[i]]<-sample_doses(EPID=m1$EPID[i],dob=m1$dob[i],onset_date=m1$date[i],investigation_date=m1$investigation_date[i],
                         guid=m1$GUID[i],date_last_opv=m1$last_opv_date[i],doses=m1$doses_opv_sia[i],offset=define.offset)
  }
  x<-do.call("rbind.data.frame",r)
  x<-as.data.frame(x)
  for(i in 1:ncol(x)){x[,i]<-unlist(x[,i])}
  m1_boot<-merge(m1,x,by="EPID")
  
  write.csv(m1_boot,paste0("data/matched_distance_",define.distance,
                      "_age_",define.age,
                      "_onset_",define.onset,
                      ifelse(match.on.RI,"_RI_",""),
                      "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),"_boot_",boot,".csv"),row.names=F)}