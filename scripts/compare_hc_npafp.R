# Map cases matched to healthy controls versus non polio AFP
define.offset=28
use.date.of.last.opv=T

filename<-paste0("data/matched_HC_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv")
hc<-read.csv(filename)

hc$stratum<-as.numeric(as.factor(hc$epid_root))
mod<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=hc)
x<-summary(mod)$coef
rownames(x[1:2,])
2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))

define.age=12
define.onset=30
match.on.RI=F
define.distance=50000

filename<-paste0("data/matched_distance_",define.distance,
                 "_age_",define.age,
                 "_onset_",define.onset,
                 ifelse(match.on.RI,"_RI_",""),
                 "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv")
m1<-read.csv(filename)

hc$control<-"Healthy community"
m1$control<-"Test-negative"

summarise(group_by(hc,case),any_OPV2=sum(nOPV2_ons>0|mOPV2_ons>0),n=n(),p=any_OPV2/n)
summarise(group_by(m1,case),any_OPV2=sum(nOPV2_ons>0|mOPV2_ons>0),n=n(),p=any_OPV2/n)

hc<-filter(hc,case)
m1<-filter(m1,case)

df<-merge(m1,hc,by="EPID",all=T)
df$control<-ifelse(is.na(df$control.x),df$control.y,df$control.x)
df$control[!is.na(df$control.x)&!is.na(df$control.y)]<-"Both"

df$y<-ifelse(is.na(df$CENTER_LAT.x),df$CENTER_LAT.y,df$CENTER_LAT.x)
df$x<-ifelse(is.na(df$CENTER_LON.x),df$CENTER_LON.y,df$CENTER_LON.x)
df$adm1<-ifelse(is.na(df$adm1.x),df$adm1.y,df$adm1.x)


load("data/nigeria_maps.Rdata")
df$month<-floor_date(as.Date(ifelse(is.na(df$date.x),df$date.y,df$date.x)),"month")
df$control[df$control=="Test-negative"&as.Date(df$month)<=as.Date("2021-03-01")]<-"Test-negative (pre-2021)"

df$control<-factor(df$control,levels=rev(c("Test-negative (pre-2021)","Test-negative","Both","Healthy community")),
                   labels=rev(c("Test-negative (pre-2021)","Test-negative (post-2021)","Both","Healthy community")),ordered=T)
p1<-ggplot()+geom_sf(data=adm0data,fill="white",linetype=0)+
  geom_jitter(data=df,aes(x=x,y=y,color=control),width=0.1,height=0.1,size=1,shape=1)+scale_color_hue()+
  # scale_color_viridis_d(begin=0.05,end=0.95)+
  geom_sf(data=adm1data,fill="transparent",size=0.5)+labs(color="Control type")+theme_void()+guides(color="none")
p2<-summarise(group_by(as_tibble(df),month,control),n=n())%>%ggplot()+geom_col(aes(x=month,y=n,fill=control))+
  labs(x="",y="Cases",fill="Control type")+scale_fill_hue()+
  # scale_fill_viridis_d(begin=0.05,end=0.95)+
  theme(legend.position=c(0.3,0.8))
p3<-summarise(group_by(as_tibble(df),adm1,control),n=n())%>%ggplot()+geom_col(aes(x=gsub("Fct","FCT",str_to_title(adm1)),y=n,fill=control))+
  labs(x="",y="Cases",fill="Control type")+scale_fill_hue()+coord_flip()+guides(fill="none")
p<-ggarrange(p2,p1,p3,labels="auto",ncol=3,nrow=1,common.legend=T)
ggsave("figs/compare_hc_npafp.png",p,width=11,height=4)


define.offset=28
use.date.of.last.opv=T

filename<-paste0("data/matched_HC_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv")
hc<-read.csv(filename)
df<-unique(hc[!hc$case,c("date","investigation_date","epid_root")])
# Mean/median time of case investigation to survey (days)
mean(as.numeric(as.Date(df$investigation_date)-as.Date(df$date)))
table(floor(as.numeric(as.Date(df$investigation_date)-as.Date(df$date))/365*12))
median(as.numeric(as.Date(df$investigation_date)-as.Date(df$date)))
summary(as.numeric(as.Date(df$investigation_date)-as.Date(df$date)))

hc$age<-floor(as.numeric(as.Date(hc$investigation_date)-as.Date(hc$dob))/365*12)
hc$age[!hc$case]<-floor(as.numeric(as.Date(hc$date)-as.Date(hc$dob))/365*12)[!hc$case]


define.age=12
define.onset=30
match.on.RI=F
define.distance=50000

filename<-paste0("data/matched_distance_",define.distance,
                 "_age_",define.age,
                 "_onset_",define.onset,
                 ifelse(match.on.RI,"_RI_",""),
                 "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv")
m1<-read.csv(filename)


get_epid_root<-function(epid){
  first_half<-substring(epid,first=1,last=nchar(epid)-5)
  second_half<-substring(epid,first=nchar(epid)-4,last=nchar(epid))
  second_half<-strsplit(second_half,"C")[[1]][1]
  return(paste0(first_half,second_half))
}
m1$epid<-gsub("-","",m1$case_epid,fixed=T)
m1$epid_root<-sapply(m1$epid,get_epid_root)

a<-unique(m1$epid_root[m1$case])
b<-unique(hc$epid_root[hc$case])
c<-intersect(a,b)

m1<-m1[m1$stratum%in%m1$stratum[m1$epid_root%in%c],]
hc<-hc[hc$epid_root%in%c,]

hc$control<-"Healthy community"
m1$control<-"Test-negative"

hc$dob_orig<-as.character(hc$dob_orig)
df<-bind_rows(hc[!hc$case,],m1)

df$bOPV_SIAs_ons<-df$OPV_SIAs_ons-df$mOPV2_SIAs_ons-df$nOPV2_SIAs_ons

df$d<-as.numeric(as.Date(ifelse(df$class=="Community control",df$investigation_date,df$date)))

a<-merge(df[df$class=="Control",c("epid_root","age","d")],df[df$class=="Case",c("epid_root","age","d")],by="epid_root",all.x=T)
b<-merge(df[df$class=="Community control",c("epid_root","age","d")],df[df$class=="Case",c("epid_root","age","d")],by="epid_root",all.x=T)

rbind(cbind(a,class="TN"),cbind(b,class="HC"))%>%mutate(age=age.x-age.y,d=d.x-d.y)%>%pivot_longer(cols=c("age","d"))%>%group_by(class,name)%>%
  summarise(median=median(abs(value),na.rm=T),lwrIQR=quantile(abs(value),p=0.25,na.rm=T),uprIQR=quantile(abs(value),p=0.75,na.rm=T))

table(df$class)
table1<-pivot_longer(df,cols=c("doses_opv_ri","doses_ipv012","doses_opv_sia","mOPV2_ons","nOPV2_ons","bOPV_SIAs_ons","mOPV2_SIAs_ons","nOPV2_SIAs_ons","age"))%>%group_by(class,name)%>%summarise(mean=round(mean(value,na.rm=T),2),sd=round(sd(value,na.rm=T),2),median=median(value,na.rm=T))%>%pivot_wider(names_from="class",values_from=c("mean","sd","median"))
View(cbind(table1[c(5,4,3,8,6,9,7,2,1),1],
           round(table1[c(5,4,3,8,6,9,7,2,1),paste0(rep(c("median_","mean_","sd_"),times=3),rep(c("Case","Control","Community control"),each=3))],2)))

df$age_cat<-as.character(cut(df$age,c(0,12,24,36,60),right=F))
df$doses_opv_ri<-as.character(ifelse(df$doses_opv_ri>3,3,df$doses_opv_ri))
df$doses_ipv012<-as.character(df$doses_ipv012)
table2<-pivot_longer(df,cols=c("doses_opv_ri","doses_ipv012","sex","age_cat"))%>%group_by(class,name,value)%>%summarise(n=n())%>%group_by(class,name)%>%
  mutate(total=sum(n),p=paste0("(",round(100*n/total),"%)"))%>%pivot_wider(names_from="class",values_from=c("p","n","total"))
View(table2[c(13,14,1:4,8:12,5:7),c("name","value",paste0(rep(c("n_","p_"),times=3),rep(c("Case","Control","Community control"),each=2)))])

df%>%filter(class!="Case")%>%ggplot()+geom_histogram(aes(x=doses_opv_sia),binwidth=1)+facet_wrap(~class)


hc$stratum<-as.numeric(as.factor(hc$epid_root))
m1$stratum<-as.numeric(as.factor(m1$epid_root))
mod_hc<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=hc)
mod_np<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)

x<-summary(mod_hc)$coef
y<-summary(mod_np)$coef
rownames(x[1:3,])
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))
2*(1-pnorm(abs(x[3,1]-y[2,1])/sqrt(x[3,3]^2 + y[3,3]^2 )))

fn<-function(mod){
  ci=exp(confint(mod))
  df<-cbind.data.frame(or=exp(mod$coefficients),lwr=ci[,1],upr=ci[,2],p=summary(mod)$coef[,5])
  df$eff<-1-df$or
  df$eff_lwr<-1-df$upr
  df$eff_upr<-1-df$lwr
  df$or_print<-round(df$or,2)
  df$ci<-paste0("(",round(df$lwr,2),", ",round(df$upr,2),")")
  df$p<-ifelse(signif(df$p,2)<0.0001,"<0.0001",signif(df$p,2))
  df$eff_print<-paste0(round(100*df$eff),"")
  df$eff_ci<-paste0("(",round(100*df$eff_lwr),", ",round(100*df$eff_upr),")")
  return(df)
}

View(cbind(fn(mod_np)[,c(8,9,4,10,11)],fn(mod_hc)[,c(8,9,4,10,11)]))


cor.test(hc$doses_opv_sia[!hc$case],hc$OPV_SIAs[!hc$case])
cor.test(m1$doses_opv_sia[!m1$case],m1$OPV_SIAs[!m1$case])

table((hc$doses_opv_sia>hc$OPV_SIAs)[!hc$case])
table((m1$doses_opv_sia>m1$OPV_SIAs)[!m1$case])

