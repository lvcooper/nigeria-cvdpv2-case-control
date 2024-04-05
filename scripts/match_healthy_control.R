# Match healthy controls

# Exclusion criteria
# q<-read.csv("data/2017_2022_AFP information_230223.csv")
# q<-q[which(q$Finalclassification=="cVDPV2"),]
# q$date<-as.Date(q$Onset.Date,"%d-%b-%y")
# q<-q[which(!is.na(q$date)),]
# q<-q[which(q$date>=as.Date("2021-03-01")),]
# 
# e<-q$EPID[as.Date(q$Onset.Date,"%d-%b-%y")>=as.Date("2021-03-01")&q$Finalclassification=="cVDPV2"]
# f<-q$EPID[q$Finalclassification=="cVDPV2"]

# AFP database
a<-read.csv("data/all.csv")
a<-a[which(a$class=="Case"),]
a<-a[which(!is.na(a$date)),]
a<-a[which(a$date>=as.Date("2021-03-01")),]

# Healthy controls
h<-read.csv("data/healthy_control_data.csv")
h$class<-"Community control"

# Calculate DOB
h$dob<-as.Date(h$dob)
h$dob[is.na(h$dob)]<-(as.Date(h$survey_date)-h$age_month*365/12)[is.na(h$dob)]

# Calculate IPV
# IPV SIAs between birth and survey
# Add SIA data
sia<-read.csv("data/SIAs_Nigeria_case_control.csv")
sia$date<-as.Date(sia$date)
sia<-sia[sia$vaccine_type%in%c("IPV","f-IPV",'IPV + bOPV'),]

h$IPV_SIAs<-sapply(1:nrow(h),function(i){
  dob=h$dob[i];investigation_date=h$survey_date[i];guid=h$GUID[i]
  if(is.na(dob)|is.na(investigation_date)|is.na(guid)) return(NA) else{
    guid<-gsub("{","",gsub("}","",tolower(guid),fixed=T),fixed=T)
    sia$age<-floor(as.numeric(sia$date-dob)/(365/12)) # Age at time of SIA
    # All eligible SIAs between birth and investigation
    sub<-sia[which(sia$date>dob&
                     sia$date<investigation_date&
                     sia$admin2guid%in%guid& 
                     sia$age<=sia$age_max&
                     sia$age>=sia$age_min),]
    sub<-unique(sub[,c("date","vaccine_type")])
    return(nrow(sub))}})

h$doses_ipv<-(h$doses_ipv_ri+ifelse(is.na(h$doses_ipv_sia)&h$IPV_SIAs==0,0,h$doses_ipv_sia))
table(is.na(h$doses_ipv_sia),h$IPV_SIAs==0)

# Calculate OPV from birth to case onset
# We haven't asked to distinguish between NA doses OPV SIA since because no SIAs or because can't remember
h$doses_opv_sia<-(h$doses_opv_sia_overall-ifelse(is.na(h$doses_opv_sia_since),0,h$doses_opv_sia_since))
h$doses_opv_sia[h$doses_opv_sia<0]<-NA

a$dob<-as.Date(a$dob)

# Match variable names between case and HC
names(h)[names(h)=="doses_opv_sia"]<-"doses_opv_sia_before"
names(h)[names(h)=="doses_opv_sia_overall"]<-"doses_opv_sia"
names(h)[names(h)=="case.epid"]<-"epid"
names(h)[names(h)=="survey_date"]<-"investigation_date"
names(h)[names(h)=="date_last_opv_sia"]<-"last_opv_date"
names(h)[names(h)=="age_month"]<-"age"

# Get case onset date for HC
# a$epid<-gsub("-","",a$EPID,fixed=T)
# h<-merge(h,a[,c("epid","date","dob")],by="epid",all.x=T,suffixes=c("","_case"))
get_epid_root<-function(epid){
  first_half<-substring(epid,first=1,last=nchar(epid)-5)
  second_half<-substring(epid,first=nchar(epid)-4,last=nchar(epid))
  second_half<-strsplit(second_half,"C")[[1]][1]
  return(paste0(first_half,second_half))
}
a$epid<-gsub("-","",a$EPID,fixed=T)
a$epid_root<-sapply(a$epid,get_epid_root)
h$epid_root<-sapply(h$epid,get_epid_root)
h<-merge(h,a[,c("epid_root","date","dob")],by="epid_root",all.x=T,suffixes=c("","_case"))

# Check date quality
# Onset date should come before investigation date
h$diff1<-as.numeric(as.Date(h$investigation_date)-as.Date(h$date))
# Date last OPV should come before investigation date
h$diff2<-as.numeric(as.Date(h$investigation_date)-as.Date(h$last_opv_date))
# Date of birth should come before onset date
h$diff3<-as.numeric(as.Date(h$date)-as.Date(h$dob))
# Date of birth should come before investigation date
h$diff4<-as.numeric(as.Date(h$investigation_date)-as.Date(h$dob))

h$dates_correct<-h$diff1>0&h$diff3>0&h$diff4>0
h$dates_correct[which(is.na(h$date))]<-T

# Exclusion flow chart
library(tidyverse)
d<-bind_rows(a,h)
d$case<-d$class=="Case"

x<-unique(d$epid[d$case])
y<-unique(d$epid[!d$case])
z<-unique(d$File.Name[!d$case])
# How many cases since March 2021?
length(x)
# How many healthy control surveys?
length(y)
length(z)
# How many controls represented in these surveys?
table(d$case)
# How many cases with no survey?
length(x[!x%in%y])
# How many surveys with no case?
length(y[!y%in%x])
# How many cases with a survey?
length(intersect(x,y))
# How many of the surveys with no case are around contacts?
table(nchar(y[!y%in%x]))
length(unique(substring(y,1,14)))


afp<-read.csv("data/all.csv")
afp$epid<-gsub("-","",afp$EPID,fixed=T)
(y[!y%in%x])[nchar(y[!y%in%x])<15]
afp[afp$epid%in%(y[!y%in%x])[nchar(y[!y%in%x])<15],]

total<-table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""])))))

d<-d[!is.na(d$investigation_date),]
flow<-rbind(total,inv_date=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))

d<-d[!is.na(d$GUID),]
flow<-rbind(flow,guid=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))


d<-d[!is.na(d$date)|d$case==0,]
flow<-rbind(flow,ons_date=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))

d<-d[!is.na(d$dob),]
flow<-rbind(flow,known_age=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))



head(d[!d$dates_correct,c("case","date","investigation_date","last_opv_date","dob")])

d<-d[d$dates_correct,]
flow<-rbind(flow,dates_correct=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))

d<-d[!is.na(d$doses_opv_sia),]
flow<-rbind(flow,doses_opv_sia=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))

d$doses_ipv012<-d$doses_ipv
d$doses_ipv012[d$doses_ipv>2]<-2

d$age<-floor(as.numeric(as.Date(d$investigation_date)-as.Date(d$dob))/365*12)

d<-d[!is.na(d$doses_ipv),]
flow<-rbind(flow,doses_ipv=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))

d<-d[d$dob>as.Date("2016-05-01"),]
flow<-rbind(flow,dob_switch=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))
d<-d[d$age<60,]
flow<-rbind(flow,under_5=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))


write.csv(d,"data/eligible_HC.csv",row.names=F)

x<-unique(d$epid[d$case])
y<-unique(d$epid[!d$case])
z<-unique(d$File.Name[!d$case])
# How many eligible cases since March 2021?
length(x)
# How many eligible healthy control surveys?
length(y)
length(z)
# How many controls represented in these surveys?
table(d$case)
# How many eligible cases with no survey?
length(x[!x%in%y])
# How many surveys with no eligible case?
length(y[!y%in%x])
# How many eligible cases with a survey?
length(intersect(x,y))
# How many of the surveys with no case are around controls?
table(nchar(y[!y%in%x]))
length(unique(substring(y,1,14)))

afp<-read.csv("data/all.csv")
afp$epid<-gsub("-","",afp$EPID,fixed=T)
afp[afp$epid%in%(y[!y%in%x])[nchar(y[!y%in%x])<15],]

# Subset to cases and HC with matches
# Remove contact suffix if required
get_epid_root<-function(epid){
  first_half<-substring(epid,first=1,last=nchar(epid)-5)
  second_half<-substring(epid,first=nchar(epid)-4,last=nchar(epid))
  second_half<-strsplit(second_half,"C")[[1]][1]
  return(paste0(first_half,second_half))
}
d$epid_root<-sapply(d$epid,get_epid_root)

x<-unique(d$epid_root[d$case])
y<-unique(d$epid_root[!d$case])
length(x[!x%in%y])
length(y[!y%in%x])
length(intersect(x,y))

d<-d[d$epid_root%in%intersect(x,y),]
flow<-rbind(flow,match=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))

# Subset to HC within 12 months of age of cases
d<-d[-which(abs(as.numeric(d$dob-d$dob_case))>(365/12*13)),]
x<-unique(d$epid_root[d$case])
y<-unique(d$epid_root[!d$case])
d<-d[d$epid_root%in%intersect(x,y),]
flow<-rbind(flow,match_age=table(c(d$class,rep("Surveys",length(unique(d$File.Name[d$File.Name!=""]))))))
flow

write.csv(d,"data/matched_HC.csv")

write.csv(flow,"data/flow_HC.csv")

x<-read.csv("data/flow_HC.csv")
x
rev(diff(rev(x$Case)))
round(100*rev(diff(rev(x$Case)))/x$Case[1:(nrow(x)-1)])

rev(diff(rev(x$Community.control)))
round(100*rev(diff(rev(x$Community.control)))/x$Community.control[1:(nrow(x)-1)])

rev(diff(rev(x$Surveys)))
round(100*rev(diff(rev(x$Surveys)))/x$Surveys[1:(nrow(x)-1)])
