regression<-function(define.distance,define.age,define.onset,define.offset,match.on.RI=F,match.limit=NA,IPV,use.date.of.last.opv=T,define.ceiling_OPV=F){
  filename<-paste0("data/matched_distance_",define.distance,
                   "_age_",define.age,
                   "_onset_",define.onset,
                   ifelse(match.on.RI,"_RI",""),
                   ifelse(!is.na(match.limit),paste0("_limit_",match.limit),""),
                   "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),ifelse(define.ceiling_OPV,"_ceiling",""),".csv")
  m1<-read.csv(filename)
  
  if(IPV==T) {
  m1<-m1[!is.na(m1$doses_ipv),]
  m1<-m1[m1$stratum%in%m1$stratum[m1$case==1],]
  m1<-m1[m1$stratum%in%m1$stratum[m1$case==0],]
  }

if(IPV==T) mod<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
if(IPV==F) mod<-clogit(case~mOPV2_ons+nOPV2_ons+strata(stratum),data=m1)

m1$doses_ipv01<-as.numeric(m1$doses_ipv012>0)
m1$bOPV_SIAs_ons<-m1$OPV_SIAs_ons-m1$mOPV2_SIAs_ons-m1$nOPV2_SIAs_ons
m1$age_cat<-as.character(cut(m1$age,c(0,12,24,36,60),right=F))
tab1<-pivot_longer(m1,cols=c("doses_opv_ri","doses_ipv012","doses_opv_sia","mOPV2_ons","nOPV2_ons","bOPV_SIAs_ons","mOPV2_SIAs_ons","nOPV2_SIAs_ons","IPV_SIAs","age"))%>%group_by(class,name)%>%summarise(mean=round(mean(value,na.rm=T),2),sd=round(sd(value,na.rm=T),2),median=median(value,na.rm=T))%>%pivot_wider(names_from="class",values_from=c("mean","sd","median"))

m1$doses_opv_ri<-as.character(m1$doses_opv_ri)
m1$doses_ipv012<-as.character(m1$doses_ipv012)
tab2<-pivot_longer(m1,cols=c("doses_opv_ri","doses_ipv012","sex","age_cat"))%>%group_by(class,name,value)%>%summarise(n=n())%>%group_by(class,name)%>%
  mutate(total=sum(n),p=paste0("(",round(100*n/total),"%)"))%>%pivot_wider(names_from="class",values_from=c("p","n","total"))


return(list(model=mod,table1=tab1,table2=tab2))}