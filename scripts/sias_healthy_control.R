define.offset=28
use.date.of.last.opv=T

# Add SIA data
sia<-read.csv("data/SIAs_Nigeria_case_control.csv")
sia$date<-as.Date(sia$date)


# Read in matched combined case and healthy control data
 
 m1<-read.csv("data/matched_HC.csv")
 m1$date<-as.Date(m1$date)
 m1$investigation_date<-as.Date(m1$investigation_date)
 m1$dob<-as.Date(m1$dob)
 m1$last_opv_date<-as.Date(m1$last_opv_date)
 
 m1$OPV_SIAs<-sapply(1:nrow(m1),function(i){
   dob=m1$dob[i];investigation_date=m1$investigation_date[i];guid=m1$GUID[i]
   guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
   sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
   # All eligible SIAs between birth and investigation
   sub<-sia[which(sia$date>dob&
                    sia$date<investigation_date&
                    sia$admin2guid%in%guid& 
                    sia$age<=sia$age_max&
                    sia$age>=sia$age_min),]
   sub<-unique(sub[!sub$vaccine_type%in%c("IPV","f-IPV"),c("date","vaccine_type")])
   return(nrow(sub))})
 
 m1$fIPV_SIAs<-sapply(1:nrow(m1),function(i){
   dob=m1$dob[i];investigation_date=m1$investigation_date[i];guid=m1$GUID[i]
   guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
   sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
   # All eligible SIAs between birth and investigation
   sub<-sia[which(sia$date>dob&
                    sia$date<investigation_date&
                    sia$admin2guid%in%guid& 
                    sia$age<=sia$age_max&
                    sia$age>=sia$age_min),]
   sub<-unique(sub[sub$vaccine_type%in%c("f-IPV"),c("date","vaccine_type")])
   return(nrow(sub))})
 
 m1$IPV_SIAs<-sapply(1:nrow(m1),function(i){
   dob=m1$dob[i];investigation_date=m1$investigation_date[i];guid=m1$GUID[i]
   guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
   sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
   # All eligible SIAs between birth and investigation
   sub<-sia[which(sia$date>dob&
                    sia$date<investigation_date&
                    sia$admin2guid%in%guid& 
                    sia$age<=sia$age_max&
                    sia$age>=sia$age_min),]
   sub<-unique(sub[sub$vaccine_type%in%c("IPV"),c("date","vaccine_type")])
   return(nrow(sub))})
 
 m1$bOPV_IPV_SIAs<-sapply(1:nrow(m1),function(i){
   dob=m1$dob[i];investigation_date=m1$investigation_date[i];guid=m1$GUID[i]
   guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
   sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
   # All eligible SIAs between birth and investigation
   sub<-sia[which(sia$date>dob&
                    sia$date<investigation_date&
                    sia$admin2guid%in%guid& 
                    sia$age<=sia$age_max&
                    sia$age>=sia$age_min),]
   sub<-unique(sub[sub$vaccine_type%in%c("IPV + bOPV"),c("date","vaccine_type")])
   return(nrow(sub))})
 
 x<-lapply(1:nrow(m1),function(i){
   dob=m1$dob[i];investigation_date=m1$investigation_date[i];guid=m1$GUID[i];date_last_opv=m1$last_opv_date[i]
   guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
   sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
   # All eligible SIAs between birth and investigation
   sub<-sia[which(sia$date>dob&
                    sia$date<investigation_date&
                    sia$admin2guid%in%guid& 
                    sia$age<=sia$age_max&
                    sia$age>=sia$age_min),]
   sub$vaccine_type[sub$vaccine_type%in%c("IPV + bOPV")]<-"bOPV"
   sub<-unique(sub[!sub$vaccine_type%in%c("IPV","f-IPV"),c("date","vaccine_type","activity_parent_code")])
   # Most recent SIA 
   a=ifelse(nrow(sub)==0,"No SIAs",sub$activity_parent_code[which.max(sub$date)])
   b=ifelse(nrow(sub)==0,"No SIAs",sub$vaccine_type[which.max(sub$date)])
   # SIA consistent with reported date
   c=ifelse(nrow(sub)==0,"No SIAs",ifelse(as.numeric(min(abs(sub$date-date_last_opv)))<=30,sub$activity_parent_code[which.min(abs(sub$date-date_last_opv))],"Inaccurate date"))
   d=ifelse(nrow(sub)==0,"No SIAs",ifelse(as.numeric(min(abs(sub$date-date_last_opv)))<=30,sub$vaccine_type[which.min(abs(sub$date-date_last_opv))],"Inaccurate date"))
   e=ifelse(nrow(sub)==0,"No SIAs",ifelse(as.numeric(min(abs(sub$date-date_last_opv)))<=30,sub$date[which.min(abs(sub$date-date_last_opv))],"Inaccurate date"))
   return(cbind.data.frame(a,b,c,d,e))})
 x<-do.call("rbind",x)
 names(x)<-c("most_recent_SIA_code","most_recent_SIA_type","reported_recent_SIA_code","reported_recent_SIA_type","reported_recent_SIA_date")
 m1<-cbind(m1,x)
 rm(x)
 
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
    return(c(sum(r[,"nOPV2"]*r$prob),sum(r[,"mOPV2"]*r$prob)))
  }
  
  sampler<-function(doses,sias,before){
    # doses=3
    # sias=c("mOPV2","bOPV","bOPV","nOPV2")
    # before=c(T,T,T,F)
    sias_before=sias
    sias_before[!before]=""
    m<-c(0,0,0,0)
    m[1:2]<-possible_histories(v=doses,s=table(factor(sias,levels=c("","bOPV","mOPV2","nOPV2"))))
    m[3:4]<-possible_histories(v=doses,s=table(factor(sias_before,levels=c("","bOPV","mOPV2","nOPV2"))))
    return(m)
  }
  
  
  sample_doses<-function(guid,dob,onset_date,investigation_date,date_last_opv=NA,doses,offset){
    # dob=m1$dob[i];onset_date=m1$date[i];investigation_date=m1$survey_date[i];guid=m1$GUID[i];date_last_opv=m1$last_opv_date[i];doses=m1$doses_opv_sia[i];offset=28
    guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
    sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
    # if give date of last dose but report 0 doses, give 1 SIA dose
    if(!is.na(date_last_opv)&doses==0) doses<-1 # this represents 1 case/control
    # All eligible SIAs between birth and onset of linked case
    sub<-sia[which(sia$date>dob&
                     sia$date<investigation_date&
                     sia$admin2guid%in%guid& 
                     sia$age<=sia$age_max&
                     sia$age>=sia$age_min),]
    sub$vaccine_type[sub$vaccine_type%in%c("IPV + bOPV")]<-"bOPV"
    sub<-unique(sub[!sub$vaccine_type%in%c("IPV","f-IPV"),c("date","vaccine_type")])
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
                              OPV_SIAs_ons=0,
                              mOPV2_SIAs_ons=0,
                              nOPV2_SIAs_ons=0,
                              date_used=F,type="",calendar_date=NA))
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
                              date_used=date_used,type=type,calendar_date=calendar_date))
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
                              date_used=date_used,type=type,calendar_date=calendar_date))
    }
    # Otherwise sample remaining doses
    if(doses_remaining<nrow(sias_remaining)){
      m<-sampler(doses=doses_remaining,sias=sias_remaining$vaccine_type,before=sias_remaining$date<=(onset_date-offset))
      return(cbind.data.frame(nOPV2=m[1]+sum(known$vaccine_type=="nOPV2"),
                              mOPV2=m[2]+sum(known$vaccine_type=="mOPV2"),
                              nOPV2_ons=m[3]+sum(known$vaccine_type=="nOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_ons=m[4]+sum(known$vaccine_type=="mOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_SIAs=sum(exp$vaccine_type=="mOPV2"),
                              nOPV2_SIAs=sum(exp$vaccine_type=="nOPV2"),
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="mOPV2"),
                              nOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="nOPV2"),
                              date_used=date_used,type=type,calendar_date=calendar_date))
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
      return(cbind.data.frame(nOPV2=m[1]+sum(known$vaccine_type=="nOPV2"),
                              mOPV2=m[2]+sum(known$vaccine_type=="mOPV2"),
                              nOPV2_ons=m[3]+sum(known$vaccine_type=="nOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_ons=m[4]+sum(known$vaccine_type=="mOPV2"&known$date<=(onset_date-offset)),
                              mOPV2_SIAs=sum(exp$vaccine_type=="mOPV2"),
                              nOPV2_SIAs=sum(exp$vaccine_type=="nOPV2"),
                              OPV_SIAs_ons=nrow(exp_ons),
                              mOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="mOPV2"),
                              nOPV2_SIAs_ons=sum(exp_ons$vaccine_type=="nOPV2"),
                              date_used=date_used,type=type,calendar_date=calendar_date))
    }
  }
  
  r<-NULL
  for(i in 1:nrow(m1)){
    r[[i]]<-sample_doses(dob=m1$dob[i],onset_date=m1$date[i],investigation_date=m1$investigation_date[i],guid=m1$GUID[i],date_last_opv=m1$last_opv_date[i],doses=m1$doses_opv_sia[i],offset=define.offset)
  }
  x<-do.call("rbind.data.frame",r)
  x<-as.data.frame(x)
  for(i in 1:ncol(x)){x[,i]<-unlist(x[,i])}
  m1<-cbind(m1,x)
  write.csv(m1,paste0("data/matched_HC_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv"),row.names=F)



