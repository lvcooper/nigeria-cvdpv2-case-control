define.offset=28
use.date.of.last.opv=T
IPV=T

filename<-paste0("data/matched_HC_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv")
m1<-read.csv(filename)
m1$stratum<-as.numeric(as.factor(m1$epid_root))

if(IPV==T) {
  m1<-m1[!is.na(m1$doses_ipv012),]
  m1<-m1[m1$stratum%in%m1$stratum[m1$case==1],]
  m1<-m1[m1$stratum%in%m1$stratum[m1$case==0],]
}
library(dplyr)
library(survival)

if(IPV==T) mod<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
if(IPV==F) mod<-clogit(case~mOPV2_ons+nOPV2_ons+strata(stratum),data=m1)

m1$doses_ipv01<-as.numeric(m1$doses_ipv012>0)
m1$bOPV_SIAs_ons<-m1$OPV_SIAs_ons-m1$mOPV2_SIAs_ons-m1$nOPV2_SIAs_ons
m1$age<-floor(as.numeric(as.Date(m1$investigation_date)-as.Date(m1$dob))/365*12)
m1$age[!m1$case]<-floor(as.numeric(as.Date(m1$date)-as.Date(m1$dob))/365*12)[!m1$case]
m1$age_cat<-as.character(cut(m1$age,c(0,12,24,36,60),right=F))
tab1<-pivot_longer(m1,cols=c("doses_opv_ri","doses_ipv012","doses_opv_sia","mOPV2_ons","nOPV2_ons","bOPV_SIAs_ons","mOPV2_SIAs_ons","nOPV2_SIAs_ons","OPV_SIAs","age"))%>%group_by(class,name)%>%summarise(mean=round(mean(value,na.rm=T),2),sd=round(sd(value,na.rm=T),2),median=median(value,na.rm=T))%>%pivot_wider(names_from="class",values_from=c("mean","sd","median"))

m1$doses_opv_ri<-as.character(ifelse(m1$doses_opv_ri>3,3,m1$doses_opv_ri))
m1$doses_ipv012<-as.character(m1$doses_ipv012)
tab2<-pivot_longer(m1,cols=c("doses_opv_ri","doses_ipv012","sex","age_cat"))%>%group_by(class,name,value)%>%summarise(n=n())%>%group_by(class,name)%>%
  mutate(total=sum(n),p=paste0("(",round(100*n/total),"%)"))%>%pivot_wider(names_from="class",values_from=c("p","n","total"))


results<-list(model=mod,table1=tab1,table2=tab2)
View(cbind(results$table1[c(5,4,3,8,6,9,7,2,1),1],
           round(results$table1[c(5,4,3,8,6,9,7,2,1),c(6,2,4,7,3,5)],2)))
View(results$table2[c(13,14,1:4,8:12,5:7),c(1,2,5,3,6,4)])

summary(mod)

m1$mOPV2_round<-ifelse(floor(m1$mOPV2_ons)>4,4,floor(m1$mOPV2_ons))
m1$nOPV2_round<-ifelse(floor(m1$nOPV2_ons)>4,4,floor(m1$nOPV2_ons))
mod<-clogit(case~as.factor(mOPV2_round)+as.factor(nOPV2_round)+as.factor(doses_ipv012)+strata(stratum),data=m1)
summary(mod)

# Analysis of reported date of last OPV SIA
temp<-m1

x<-table(ifelse(is.na(temp$last_opv_date),"Date of last OPV dose not reported","Date of last OPV dose reported"),ifelse(temp$case,"Case","Control"))
x
paste0(x, " (",round(100*x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x)))),"%)")

x<-table(ifelse(!temp$date_used[!is.na(temp$last_opv_date)],"Date of last OPV reported, but > 30 days from SIA calendar",
                "Date of last OPV reported within 30 days from SIA calendar"),ifelse(temp$case[!is.na(temp$last_opv_date)],"Case","Control"))
x
paste0(x, " (",round(100*x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x)))),"%)")

temp$after<-(as.Date(temp$calendar_date))>(as.Date(temp$date)-28)

x<-table(paste0("Date of last OPV reported, consistent with ",
                temp$type[!is.na(temp$last_opv_date)&temp$date_used]),ifelse(temp$case[!is.na(temp$last_opv_date)&temp$date_used],"Case","Control"))
x
paste0(x, " (",round(100*x/sum(x)),"%)")

y<-table(paste0("Date of last OPV reported before onset, consistent with ",temp$type[!is.na(temp$last_opv_date)&temp$date_used&!temp$after]),
         ifelse(temp$case[!is.na(temp$last_opv_date)&temp$date_used&!temp$after],"Case","Control"))
y
y/x

