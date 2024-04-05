
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")

# Anyone reported zero OPV SIA doses but a date of last OPV from SIA?
sum(m1$doses_opv_sia==0&!is.na(m1$last_opv_date))

# How many date of last OPV from SIA might be consistent with RI? <2%
table(cut(as.numeric(as.Date(m1$last_opv_date)-as.Date(m1$dob)),c(-Inf,0,14*7,Inf),right=F),m1$age<=floor(14*7/30),useNA="ifany")
table(cut(as.numeric(as.Date(m1$last_opv_date)-as.Date(m1$dob)),c(-Inf,0,18*7,Inf),right=F),m1$age<=floor(18*7/30),useNA="ifany")

x<-table(m1$source[m1$doses_ipv==0],m1$class[m1$doses_ipv==0],useNA="ifany")
y<-table(m1$source[m1$doses_ipv>0],m1$class[m1$doses_ipv>0],useNA="ifany")
x
round(100*x/rbind(colSums(x),colSums(x),colSums(x)))
y
round(100*y/rbind(colSums(y),colSums(y),colSums(y)))

# Histograms
ymax<-0.8

m1$mOPV2_cat<-cut(m1$mOPV2_ons,c(seq(-0.5,5.5,1),100),right=F,labels=c(0:5,">=6"))
p8<-mutate(group_by(ungroup(summarise(group_by(m1,mOPV2_cat,case),n=n())),case),total=sum(n))%>%ggplot()+
  geom_col(aes(x=mOPV2_cat,y=n/total,fill=ifelse(case,"cVDPV2 case","Matched test-negative control")),position=position_dodge(width=0.75))+labs(x="Modelled mOPV2 doses",fill="",y="Proportion")+scale_y_continuous(limits=c(0,ymax))

m1$nOPV2_cat<-cut(m1$nOPV2_ons,c(seq(-0.5,5.5,1),100),right=F,labels=c(0:5,">=6"))
p9<-mutate(group_by(ungroup(summarise(group_by(m1,nOPV2_cat,case),n=n())),case),total=sum(n))%>%ggplot()+
  geom_col(aes(x=nOPV2_cat,y=n/total,fill=ifelse(case,"cVDPV2 case","Matched test-negative control")),position=position_dodge(width=0.75))+labs(x="Modelled nOPV2 doses",fill="",y="Proportion")+scale_y_continuous(limits=c(0,ymax))



m1$doses_opv_ri_cat<-ifelse(is.na(m1$doses_opv_ri),"Unknown",m1$doses_opv_ri)
m1$mOPV2_SIAs_ons_cat<-cut(m1$mOPV2_SIAs_ons,c(0:6,100),right=F,labels=c(0:5,">=6"))
p4<-mutate(group_by(ungroup(summarise(group_by(m1,mOPV2_SIAs_ons_cat,case),n=n())),case),total=sum(n))%>%ggplot()+
  geom_col(aes(x=mOPV2_SIAs_ons_cat,y=n/total,fill=ifelse(case,"cVDPV2 case","Matched test-negative control")),position=position_dodge(width=0.75))+labs(x="mOPV2 SIAs",fill="",y="Proportion")+scale_y_continuous(limits=c(0,ymax))
m1$nOPV2_SIAs_ons_cat<-cut(m1$nOPV2_SIAs_ons,c(0:6,100),right=F,labels=c(0:5,">=6"))
p5<-mutate(group_by(ungroup(summarise(group_by(m1,nOPV2_SIAs_ons_cat,case),n=n())),case),total=sum(n))%>%ggplot()+
  geom_col(aes(x=nOPV2_SIAs_ons_cat,y=n/total,fill=ifelse(case,"cVDPV2 case","Matched test-negative control")),position=position_dodge(width=0.75))+labs(x="nOPV2 SIAs",fill="",y="Proportion")+scale_y_continuous(limits=c(0,ymax))
m1$OPV<-cut(m1$OPV_SIAs_ons,c(seq(0,16,by=2),100),right=F,labels=c(paste0(seq(0,14,by=2),"-",seq(1,15,by=2)),">=16"))

m1$doses_opv_sia_cat<-cut(m1$doses_opv_sia,c(0:15,100),right=F,labels=c(0:14,">=15"))
p1_alt<-mutate(group_by(ungroup(summarise(group_by(m1,doses_opv_sia_cat,case),n=n())),case),total=sum(n))%>%ggplot()+
  geom_col(aes(x=doses_opv_sia_cat,y=n/total,fill=ifelse(case,"cVDPV2 case","Matched test-negative control")),position=position_dodge(width=0.75))+labs(x="OPV doses from SIAs",fill="",y="Proportion")+scale_y_continuous(limits=c(0,0.2))
m1$OPV<-cut(m1$OPV_SIAs_ons,c(0:15,100),right=F,labels=c(0:14,">=15"))
p6_alt<-mutate(group_by(ungroup(summarise(group_by(m1,OPV,case),n=n())),case),total=sum(n))%>%ggplot()+
  geom_col(aes(x=OPV,y=n/total,fill=ifelse(case,"cVDPV2 case","Matched test-negative control")),position=position_dodge(width=0.75))+labs(x="Total OPV SIAs",fill="",y="Proportion")+scale_y_continuous(limits=c(0,0.2))
m1$IPV<-cut(m1$doses_ipv,c(0:3),right=F,labels=c(0:2))

p<-ggarrange(p1_alt,p4,p8,p6_alt,p5,p9,common.legend = T,ncol=3,nrow=2,widths=c(1,0.5,0.5),labels="auto")
p

ggsave("figs/histograms.png",p,width=10,height=5)

# Type of IPV used in SIAs
mean(m1$IPV_SIAs[m1$doses_ipv>0]>0)
mean(m1$bOPV_IPV_SIAs[m1$doses_ipv>0]>0)
mean(m1$fIPV_SIAs[m1$doses_ipv>0]>0)
sum(m1$IPV_SIAs[m1$doses_ipv>0]>0)
sum(m1$bOPV_IPV_SIAs[m1$doses_ipv>0]>0)
sum(m1$fIPV_SIAs[m1$doses_ipv>0]>0)

mean(m1$bOPV_IPV_SIAs[m1$doses_ipv>0]>0|m1$IPV_SIAs[m1$doses_ipv>0]>0)
mean(m1$fIPV_SIAs[m1$doses_ipv>0]>0)
sum(m1$bOPV_IPV_SIAs[m1$doses_ipv>0]>0|m1$IPV_SIAs[m1$doses_ipv>0]>0)
sum(m1$fIPV_SIAs[m1$doses_ipv>0]>0)

# IPV: source RI versus source SIA
table(m1$doses_ipv,m1$case,useNA="ifany")
table(m1$doses_ipv_ri,m1$case,useNA="ifany")
table(m1$doses_ipv_sia,m1$case,useNA="ifany")
# We only know the breakdown of IPV for RI versus SIA for AFP with onset in 2022:
# 107 of 1303 test-negative controls, 33 of 363 cVDPV2 cases
table(!is.na(m1$doses_ipv_ri)&!is.na(m1$doses_ipv_sia),m1$case)
table(paste0("RI: ",m1$doses_ipv_ri," SIA: ",m1$doses_ipv_sia)[!is.na(m1$doses_ipv_ri)&!is.na(m1$doses_ipv_sia)],m1$case[!is.na(m1$doses_ipv_ri)&!is.na(m1$doses_ipv_sia)])
# For controls, 72 of 107 reported zero doses IPV, 31 reported 1 dose from RI, 2 reported 1 dose from SIA, and 2 reported 1 dose from each RI and SIA
# For cases, 28 of 32 reported zero doses IPV, 1 reported 1 dose from RI, 1 reported 2 doses from RI, and 3 reported 1 dose from each RI and SIA

m1$total_IPV_SIAs<-m1$IPV_SIAs+m1$fIPV_SIAs+m1$bOPV_IPV_SIAs

m1$doses_ipv_ri_assign<-m1$doses_ipv_ri
m1$doses_ipv_sia_assign<-m1$doses_ipv_sia
# If zero IPV reported, zero doses
m1$doses_ipv_sia_assign[is.na(m1$doses_ipv_sia)&m1$doses_ipv==0]<-0
m1$doses_ipv_ri_assign[is.na(m1$doses_ipv_ri)&m1$doses_ipv==0]<-0
# If no IPV SIAs, all reported doses come from RI
m1$doses_ipv_sia_assign[is.na(m1$doses_ipv_sia_assign)&is.na(m1$doses_ipv_sia)&m1$total_IPV_SIAs==0]<-0
m1$doses_ipv_ri_assign[is.na(m1$doses_ipv_ri_assign)&is.na(m1$doses_ipv_ri)&m1$total_IPV_SIAs==0]<-
  m1$doses_ipv[is.na(m1$doses_ipv_ri_assign)&is.na(m1$doses_ipv_ri)&m1$total_IPV_SIAs==0]
# If RI card and reporting 1 or more doses, at least 1 dose from RI
m1$doses_ipv_ri_assign[which(is.na(m1$doses_ipv_ri_assign)&is.na(m1$doses_ipv_ri)&m1$doses_ipv>1)]<-1

m1$doses_ipv_total<-ifelse(is.na(m1$doses_ipv),m1$doses_ipv_ri+m1$doses_ipv_sia,m1$doses_ipv)
m1$doses_ipv_unknown_assign<-m1$doses_ipv_total-ifelse(is.na(m1$doses_ipv_ri_assign),0,m1$doses_ipv_ri_assign)-ifelse(is.na(m1$doses_ipv_sia_assign),0,m1$doses_ipv_sia_assign)

m1$doses_ipv_ri_guess<-m1$doses_ipv_ri
m1$doses_ipv_sia_guess<-m1$doses_ipv_sia
table(m1$doses_ipv_ri_guess,useNA="ifany")
table(m1$doses_ipv_sia_guess,useNA="ifany")
# If zero IPV reported, zero doses
m1$doses_ipv_sia_guess[is.na(m1$doses_ipv_sia)&m1$doses_ipv==0]<-0
m1$doses_ipv_ri_guess[is.na(m1$doses_ipv_ri)&m1$doses_ipv==0]<-0
table(m1$doses_ipv_ri_guess,useNA="ifany")
table(m1$doses_ipv_sia_guess,useNA="ifany")
# If no IPV SIAs, all reported doses come from RI
m1$doses_ipv_sia_guess[is.na(m1$doses_ipv_sia_guess)&is.na(m1$doses_ipv_sia)&m1$total_IPV_SIAs==0]<-0
m1$doses_ipv_ri_guess[is.na(m1$doses_ipv_ri_guess)&is.na(m1$doses_ipv_ri)&m1$total_IPV_SIAs==0]<-
  m1$doses_ipv[is.na(m1$doses_ipv_ri_guess)&is.na(m1$doses_ipv_ri)&m1$total_IPV_SIAs==0]
table(m1$doses_ipv_ri_guess,useNA="ifany")
table(m1$doses_ipv_sia_guess,useNA="ifany")
# If RI card and reporting 1 or more doses and born before 2021, 1 dose from RI
m1$doses_ipv_ri_guess[which(is.na(m1$doses_ipv_ri_guess)&is.na(m1$doses_ipv_ri)&m1$doses_ipv>1&m1$dob<as.Date("2021-01-01"))]<-1
table(m1$doses_ipv_ri_guess,useNA="ifany")
# If RI card and reporting 1 or more doses and born after 2021, 1 or 2 dose from RI
m1$doses_ipv_ri_guess[which(is.na(m1$doses_ipv_ri_guess)&is.na(m1$doses_ipv_ri)&m1$doses_ipv>1&m1$dob>=as.Date("2021-01-01"))]<-1
table(m1$doses_ipv_ri_guess,useNA="ifany")


m1$ipv_source<-"Unknown"
# If zero doses
m1$ipv_source[which(m1$doses_ipv==0)]<-"None"
# If reporting separately
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv_ri>0&m1$doses_ipv_sia==0)]<-"Known RI"
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv_sia>0&m1$doses_ipv_ri==0)]<-"Known SIA"
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv_sia>0&m1$doses_ipv_ri>0)]<-"Known RI and SIA"
# If report 1 dose IPV and source of info is card
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$source=="Card"&m1$doses_ipv==1)]<-"Known RI"

m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv_ri_guess>0&m1$doses_ipv_sia_guess==0)]<-"Guess RI"
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv_sia_guess>0&m1$doses_ipv_ri_guess==0)]<-"Guess SIA"
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv_sia_guess>0&m1$doses_ipv_ri_guess>0)]<-"Guess RI and SIA"
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$total_IPV_SIAs<m1$doses_ipv)]<-"Guess RI and SIA"
m1$ipv_source[which(m1$ipv_source=="Unknown"&m1$doses_ipv>1&m1$dob<as.Date("2021-01-01"))]<-"Guess RI and SIA"

m1$ipv_source[which(m1$ipv_source=="Unknown")]<-"Guess RI or SIA"
m1$ipv_source_fac<-factor(m1$ipv_source,levels=rev(c("Known RI","Guess RI","Known RI and SIA","Guess RI and SIA","Guess RI or SIA","Known SIA","None")),ordered=T)
m1$case_print<-ifelse(m1$case,"cVDPV2 case","Test-negative control")
summarise(group_by(m1,ipv_source_fac,case_print),n=n())%>%group_by(case_print)%>%mutate(total=sum(n))%>%
  ggplot()+geom_col(aes(x=case_print,y=n/total,fill=ipv_source_fac))+labs(x="",y="Proportion",fill="Source")+scale_y_continuous(labels=scales::percent)
ggsave("figs/ipv_source.png",width=5,height=5)

summarise(group_by(m1%>%filter(doses_ipv>0),case_print),
          low_ri=sum(grepl("RI",ipv_source)&grepl("Known",ipv_source)),
          high_ri=sum(grepl("RI",ipv_source)),
          low_sia=sum(grepl("SIA",ipv_source)&grepl("Known",ipv_source)),
          high_sia=sum(grepl("SIA",ipv_source)),
          n=n())
table(grepl("SIA",m1$ipv_source[m1$doses_ipv>0]),m1$case_print[m1$doses_ipv>0])
table(grepl("SIA",m1$ipv_source[m1$doses_ipv>0])&grepl("Known",m1$ipv_source[m1$doses_ipv>0]),m1$case_print[m1$doses_ipv>0])
table(grepl("RI",m1$ipv_source[m1$doses_ipv>0]),m1$case_print[m1$doses_ipv>0])
table(grepl("RI",m1$ipv_source[m1$doses_ipv>0])&grepl("Known",m1$ipv_source[m1$doses_ipv>0]),m1$case_print[m1$doses_ipv>0])

table(m1$source[m1$doses_ipv>0],m1$case_print[m1$doses_ipv>0],useNA="ifany")

table(m1$IPV_SIAs[m1$doses_ipv>0&grepl("SIA",m1$ipv_source)]>0|
        m1$bOPV_IPV_SIAs[m1$doses_ipv>0&grepl("SIA",m1$ipv_source)]>0,m1$case_print[m1$doses_ipv>0&grepl("SIA",m1$ipv_source)])
table(m1$fIPV_SIAs[m1$doses_ipv>0&grepl("SIA",m1$ipv_source)]>0,m1$case_print[m1$doses_ipv>0&grepl("SIA",m1$ipv_source)])

x<-summarise(group_by(m1,ipv_source_fac),n=n())%>%mutate(total=sum(n),p=paste0(round(100*n/total),"%"))
y<-summarise(group_by(m1,ipv_source_fac,case_print),n=n())%>%group_by(case_print)%>%mutate(total=sum(n),p=paste0(round(100*n/total),"%"))%>%pivot_wider(names_from=case_print,values_from=c(n,total,p))
View(merge(x,y,by="ipv_source_fac")[c(6,7,3,1,4,2,5),c(1,5,9,6,10,2,4)])
summarise(group_by(m1,ipv_source,case_print),n=n())%>%filter(ipv_source!="None")%>%group_by(case_print)%>%mutate(total=sum(n))%>%
  ggplot()+geom_col(aes(x=case_print,y=n/total,fill=ipv_source))

table(m1$ipv_source,m1$case_print)

unique(m1[m1$ipv_source=="Unknown",c("IPV_SIAs","doses_ipv")])
table(paste0("RI: ",m1$doses_ipv_ri_guess," SIA: ",m1$doses_ipv_sia_guess),m1$case_print)

m1$case_print<-ifelse(m1$case,"cVDPV2 case","Matched test-negative control")
p<-ggplot()+geom_jitter(data=m1,aes(x=OPV_SIAs,y=doses_opv_sia),alpha=0.2,width=0.25,height=0.25)+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_smooth(data=m1,aes(x=OPV_SIAs,y=doses_opv_sia),method="loess")+
  # geom_text(data=summarise(group_by(as_tibble(m1),case_print),label=paste0(round(100*mean(doses_opv_sia>OPV_SIAs)),"%")),aes(x=12,y=23,label=label))+
  facet_wrap(~case_print)+labs(x="OPV SIAs exposed",y="Doses OPV from SIA reported")
ggsave("figs/overreporting.png",p,width=7,height=3.5)

m1$case_print<-ifelse(m1$case,"cVDPV2 case","Matched test-negative control")
m1$age_cat<-ifelse(m1$age<=24,"<=24 months",">24 months")
p<-ggplot()+geom_jitter(data=m1,aes(x=OPV_SIAs,y=doses_opv_sia,color=age_cat),alpha=0.2,width=0.25,height=0.25)+geom_abline(slope=1,intercept=0)+
  geom_smooth(data=m1,aes(x=OPV_SIAs,y=doses_opv_sia,color=age_cat),method="loess")+
  # geom_text(data=summarise(group_by(as_tibble(m1),case_print,age_cat),label=paste0(round(100*mean(doses_opv_sia>OPV_SIAs)),"%")),aes(x=6+ifelse(age_cat=="<=24 months",0,12),y=23,label=label,color=age_cat))+
  facet_wrap(~case_print)+labs(x="OPV SIAs exposed",y="Doses OPV from SIA reported",color="Age")
ggsave("figs/overreporting_by_age.png",p,width=7,height=3.5)


