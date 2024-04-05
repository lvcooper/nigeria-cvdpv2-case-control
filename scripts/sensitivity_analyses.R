source("scripts/matching.R")
source("scripts/sias.R")
source("scripts/regression.R")

print_table<-function(mod){
  points<-cbind.data.frame(med=summary(mod)$conf[,1],lwr=summary(mod)$conf[,3],upr=summary(mod)$conf[,4],var=rownames(summary(mod)$conf),pval=signif(summary(mod)$coef[,5],2))
  points$eff<-paste0(round(100*(1-points$med)),"")
  points$eff_lwr<-round(100*(1-points$upr))
  points$eff_upr<-round(100*(1-points$lwr))
  points$eff_ci<-paste0("(",points$eff_lwr,", ",points$eff_upr,")")
  points$or<-round(points$med,2)
  points$or_ci<-paste0("(",round(points$lwr,2),", ",round(points$upr,2),")")
  rownames(points)<-1:nrow(points)
  return(points[,c("var","or","or_ci","pval","eff","eff_ci")])
}


# Double-check method for calculating p-value for difference in two odds ratios
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
mod<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
summary(mod)
x<-summary(mod)$coef
rownames(x[1:2,])
2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[2,3]^2 + x[1,3]^2 )))
m1$OPV2_ons<-m1$mOPV2_ons+m1$nOPV2_ons
mod<-clogit(case~OPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
summary(mod)

# Compare to main
results<-regression(define.distance=50000,define.age=12,define.onset=30,define.offset=28,use.date.of.last.opv=T,IPV=T)
mod<-results$model
x<-summary(mod)$coef
rownames(x[1:2,])
2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[2,3]^2 + x[1,3]^2 )))
y<-summary(mod)$coef

# Sensitivity - 1 control per case
matching(define.distance=50000,define.age=12,define.onset=30,match.limit=1)
sias(define.distance=50000,define.age=12,define.onset=30,match.limit=1,define.offset=28,use.date.of.last.opv=T)
read.csv("data/flow_match_distance_50000_age_12_onset_30_limit_1.csv")
results<-regression(define.distance=50000,define.age=12,define.onset=30,match.limit=1,define.offset=28,use.date.of.last.opv=T,IPV=T)
mod1<-results$model
x<-summary(mod1)$coef
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))
print_table(mod1)

# Sensitivity - match on RI
matching(define.distance=50000,define.age=12,define.onset=30,match.on.RI=T)
sias(define.distance=50000,define.age=12,define.onset=30,match.on.RI=T,define.offset=28,use.date.of.last.opv=T)
read.csv("data/flow_match_distance_50000_age_12_onset_30_RI.csv")
results<-regression(define.distance=50000,define.age=12,define.onset=30,define.offset=28,match.on.RI=T,use.date.of.last.opv=T,IPV=F)
mod2<-results$model
x<-summary(mod2)$coef
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))
print_table(mod2)

# Sensitivity - ignore reported date of last OPV
sias(define.distance=50000,define.age=12,define.onset=30,define.offset=28,use.date.of.last.opv=F)
results<-regression(define.distance=50000,define.age=12,define.onset=30,define.offset=28,use.date.of.last.opv=F,IPV=T)
mod3<-results$model
x<-summary(mod3)$coef
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))
print_table(mod3)

# Sensitivity - max doses reported
sias(define.distance=50000,define.age=12,define.onset=30,match.on.RI=F,define.offset=28,use.date.of.last.opv=T,define.ceiling_OPV=T)
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
mean(m1$mOPV2_ons>m1$mOPV2_SIAs_ons)
mean(m1$nOPV2_ons>m1$nOPV2_SIAs_ons)
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28_ceiling.csv")
mean(m1$mOPV2_ons>m1$mOPV2_SIAs_ons)
mean(m1$nOPV2_ons>m1$nOPV2_SIAs_ons)

results<-regression(define.distance=50000,define.age=12,define.onset=30,define.offset=28,match.on.RI=F,use.date.of.last.opv=T,IPV=T,define.ceiling_OPV=T)
mod4<-results$model
x<-summary(mod4)$coef
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))
print_table(mod4)

# Sensitivity - Age under 3 years 
matching(define.distance=50000,define.age=12,define.onset=30,under=36)
sias(define.distance=50000,define.age=12,define.onset=30,under=36,define.offset=28,use.date.of.last.opv=T)
read.csv("data/flow_match_distance_50000_age_12_onset_30_under_36.csv")
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_under_36_offset_28.csv")

mod5<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
x<-summary(mod5)$coef
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))
print_table(mod5)

# Sensitivity - Cases made positive by contact
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
m1<-m1[m1$stratum%in%m1$stratum[m1$case==1&m1$pos_cont==F],]
m1<-m1[m1$stratum%in%m1$stratum[m1$case==0],]
mod6<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
x<-summary(mod6)$coef
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))

m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
m1$fac<-m1$stratum%in%m1$stratum[m1$case==1&m1$pos_cont==F]
mod6A<-clogit(case~mOPV2_ons:fac+nOPV2_ons:fac+as.factor(doses_ipv012):fac+strata(stratum),data=m1)
x<-summary(mod6A)$coef
2*(1-pnorm(abs(x[2,1]-x[1,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))
2*(1-pnorm(abs(x[3,1]-x[4,1])/sqrt(x[3,3]^2 + x[4,3]^2 )))
2*(1-pnorm(abs(x[5,1]-x[6,1])/sqrt(x[5,3]^2 + x[6,3]^2 )))
x<-print_table(mod6A)
x$fac<-ifelse(grepl("TRUE",x$var),"Positive","Contact")
View(cbind(x[x$fac=="Positive",c("var","or","or_ci","pval","eff","eff_ci")],x[x$fac=="Contact",c("or","or_ci","pval","eff","eff_ci")]))

# Sensitivity - leave out Sokoto, Zamfara
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
m1<-m1[m1$stratum%in%m1$stratum[m1$case==1&!m1$adm1%in%c("sokoto",'zamfara')],]
m1<-m1[m1$stratum%in%m1$stratum[m1$case==0],]
mod7<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
x<-summary(mod7)$coef
# difference in mOPV2 effectiveness between main model and model excluding Sokoto, Zamfara:
2*(1-pnorm(abs(x[1,1]-y[1,1])/sqrt(x[1,3]^2 + y[1,3]^2 )))
# difference in nOPV2 effectiveness between main model and model excluding Sokoto, Zamfara:
2*(1-pnorm(abs(x[2,1]-y[2,1])/sqrt(x[2,3]^2 + y[2,3]^2 )))

# Sensitivity - separate out Sokoto, Zamfara
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
m1$fac<-m1$stratum%in%m1$stratum[m1$adm1%in%c("sokoto",'zamfara')]
mod7A<-clogit(case~mOPV2_ons:fac+nOPV2_ons:fac+as.factor(doses_ipv012)+strata(stratum),data=m1)
x<-summary(mod7A)$coef
# Difference between mOPV2 in Sokoto and mOPV2 elsewhere
2*(1-pnorm(abs(x[3,1]-x[4,1])/sqrt(x[3,3]^2 + x[4,3]^2 )))
# Difference between nOPV2 in Sokoto and nOPV2 elsewhere
2*(1-pnorm(abs(x[5,1]-x[6,1])/sqrt(x[5,3]^2 + x[6,3]^2 )))
# Difference between nOPV2 in Sokoto and mOPV2 in Sokoto
2*(1-pnorm(abs(x[4,1]-x[6,1])/sqrt(x[4,3]^2 + x[6,3]^2 )))
print_table(mod7A)

# Sensitivity: IPV history by card
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
m1$doses_ipv_source<-paste0(m1$doses_ipv012,ifelse(!is.na(m1$source)&m1$source=="Card"&m1$doses_ipv012>0,"Card",""))
mod9<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv_source)+strata(stratum),data=m1)
View(print_table(mod9))
x<-summary(mod9)$coef
rownames(x)[3:4]
2*(1-pnorm(abs(x[3,1]-x[4,1])/sqrt(x[3,3]^2 + x[4,3]^2 )))

mod9<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv_source)+strata(stratum),data=m1)

m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
mod<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)
x<-summary(mod)$coef
2*(1-pnorm(abs(x[2,1]-x[1,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))

df<-do.call("rbind.data.frame",lapply(list(mod2,mod1,mod3,mod4,mod5,mod6,mod7,mod),function(i){cbind.data.frame(eff=1-summary(i)$conf[,1],lwr=1-summary(i)$conf[,3],upr=1-summary(i)$conf[,4],var=rownames(summary(i)$conf))}))
names(df)<-c("eff","lwr","upr","var")
unlist(lapply(list(mod2,mod1,mod3,mod4,mod5,mod6,mod7,mod),function(mod){
  x<-summary(mod)$coef
2*(1-pnorm(abs(x[2,1]-x[1,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))}))

df$type<-c(rep("Match on RI history",2),rep(c("Maximum 1 control per case","Ignore date last OPV SIA","Maximum doses from OPV SIA","Exclude 36 mos. or older","Exclude cases positive by contact","Exclude Sokoto, Zamfara","Main"),each=4))
df$type<-paste0(c(rep("(a) ",2),rep(paste0("(",letters[2:7],") "),each=4),rep("",4)),df$type)
df$type<-factor(df$type,levels=rev(sort(unique(df$type))),ordered=T)
df$var<-factor(df$var,levels=c("mOPV2_ons","nOPV2_ons","as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2+ doses)"),ordered=T)
ggplot()+geom_hline(data=df%>%filter(type=="Main"),aes(yintercept=eff))+
  geom_hline(data=df%>%filter(type=="Main"),aes(yintercept=lwr),linetype=2)+
  geom_hline(data=df%>%filter(type=="Main"),aes(yintercept=upr),linetype=2)+
  geom_pointrange(data=df%>%filter(type!="Main"),aes(x=type,y=eff,ymin=lwr,ymax=upr),fatten=1)+labs(x="Sensitivity analysis",y="Effectiveness")+
  facet_wrap(~var,scales="free_x")+coord_flip()+scale_y_continuous(labels=scales::percent)
ggsave("figs/visualise_regression_SA.png",width=8,height=4)

# Gender
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
summarise(group_by(m1%>%filter(!is.na(doses_opv_ri)),sex,class),mean(doses_opv_ri>2))
summarise(group_by(m1,sex,class,doses_ipv012),n=n())
summarise(group_by(m1,sex,class),mean(doses_opv_sia),sd(doses_opv_sia))

m1$age_cat<-as.character(cut(m1$age,c(0,12,24,36,60),right=F))
View(pivot_longer(m1,cols=c(doses_opv_sia,mOPV2_ons,nOPV2_ons))%>%group_by(class,name,age_cat,sex)%>%summarise(mean=round(mean(value,na.rm=T),1),sd=round(sd(value,na.rm=T),1),print=paste0(mean," (",sd,")"))%>%
       dplyr::select(class,name,sex,print,age_cat)%>%pivot_wider(names_from=c(class,sex),values_from=print)%>%arrange(name,age_cat))

m1$doses_opv_ri<-as.character(m1$doses_opv_ri)
m1$doses_ipv012<-as.character(m1$doses_ipv012)

tab<-pivot_longer(m1,cols=c(doses_opv_ri,doses_ipv012))%>%group_by(class,sex,age_cat,name,value)%>%summarise(n=n())%>%group_by(class,sex,name,age_cat)%>%
  mutate(total=sum(n),p=paste0("(",round(100*n/total),"%)"),print=paste0(n,"/",total," ",p))%>%dplyr::select(class,name,value,sex,print,age_cat)%>%pivot_wider(names_from=c(class,sex),values_from=print)
View(tab%>%arrange(age_cat,name,value))

tab<-pivot_longer(m1,cols=c(doses_opv_ri,doses_ipv012))%>%group_by(class,sex,name,value)%>%summarise(n=n())%>%group_by(class,sex,name)%>%
  mutate(total=sum(n),p=paste0("(",round(100*n/total),"%)"),print=paste0(n,"/",total," ",p))%>%dplyr::select(class,name,value,sex,print)%>%pivot_wider(names_from=c(class,sex),values_from=print)
View(tab%>%arrange(name,value))

summary(glm((doses_opv_ri>2)~sex+class,data=m1,family="binomial"))
summary(glm((doses_ipv>0)~sex+class,data=m1,family="binomial"))
summary(glm.nb(doses_opv_sia~sex+class,data=m1))

summary(clogit(case~sex+strata(stratum),data=m1))
summary(clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+sex+strata(stratum),data=m1))

# Sensitivity - whole numbers, non-linear
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")
grps=6
s<-c(0:grps,100)
labs<-paste0(0:(grps-1)," - ",1:grps)
labs<-c(labs,paste0(grps,"+"))
m1$mOPV2_cat<-as.factor(ifelse(round(m1$mOPV2_ons)>(grps-1),paste0(grps,"+"),round(m1$mOPV2_ons)))
m1$nOPV2_cat<-as.factor(ifelse(round(m1$nOPV2_ons)>(grps-1),paste0(grps,"+"),round(m1$nOPV2_ons)))

cutoff=7

mod1<-clogit(case~mOPV2_cat+nOPV2_cat+as.factor(doses_ipv012)+strata(stratum),data=m1)

mod2<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1)

sub1<-m1[!m1$stratum%in%m1$stratum[m1$nOPV2_ons>=cutoff&m1$mOPV2_ons>=cutoff],]
sub1<-m1
sub1$mOPV2_ons[sub1$mOPV2_ons>4]<-4
sub1$nOPV2_ons[sub1$nOPV2_ons>4]<-4
mod3<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=sub1)

points<-cbind.data.frame(med=summary(mod1)$conf[,1],lwr=summary(mod1)$conf[,3],upr=summary(mod1)$conf[,4],var=rownames(summary(mod1)$conf))%>%filter(grepl("OPV2",var))
points$doses<-rep(1:(nrow(points)/2),2)

lines<-rbind(cbind.data.frame(med=summary(mod2)$conf[,1],lwr=summary(mod2)$conf[,3],upr=summary(mod2)$conf[,4],var=rownames(summary(mod2)$conf),type="All data"),
             cbind.data.frame(med=summary(mod3)$conf[,1],lwr=summary(mod3)$conf[,3],upr=summary(mod3)$conf[,4],var=rownames(summary(mod3)$conf),type=paste0("OPV2 < ",cutoff)))%>%filter(grepl("OPV2",var))
lines<-merge(lines,cbind.data.frame(doses=1:grps))
lines$med<-lines$med^lines$doses
lines$upr<-lines$upr^lines$doses
lines$lwr<-lines$lwr^lines$doses
lines<-lines%>%filter(!((doses>(cutoff-1)&type==paste0("OPV2 < ",cutoff)&grepl("mOPV2",var))))
lines<-lines%>%filter(!((doses>(cutoff-1)&type==paste0("OPV2 < ",cutoff)&grepl("nOPV2",var))))
lines<-lines%>%filter(type=="All data")

ggplot()+
  geom_ribbon(data=lines,aes(x=doses,ymin=lwr,ymax=upr),alpha=0.5)+
  geom_pointrange(data=points,aes(x=doses,y=med,ymin=lwr,ymax=upr),fatten=1)+scale_y_log10()+
  geom_line(data=lines,aes(x=doses,y=med))+scale_x_continuous(breaks=1:grps,labels=levels(m1$mOPV2_cat)[2:(grps+1)])+
  scale_color_manual(values=scales::hue_pal()(4)[c(1,3,4)])+scale_fill_manual(values=scales::hue_pal()(4)[c(1,3,4)])+
  facet_wrap(~ifelse(grepl("nOPV2",var),"nOPV2","mOPV2"),scales="free_x")+labs(y="Odds ratio",x="Modelled doses",color="",fill="")
ggsave("figs/non-linear_OR.png",width=7,height=3)

# Sensitivity - matching distance
dists<-c(1,10000*c(1:4,6:10))
# Match to controls
for(i in dists){
  matching(define.distance=i,define.age=12,define.onset=30)
  print(i)
}

x<-rbind(tail(read.csv("data/flow_match_distance_1_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_10000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_20000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_30000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_40000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_50000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_60000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_70000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_80000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_90000_age_12_onset_30.csv")[11,],1),
         tail(read.csv("data/flow_match_distance_1e+05_age_12_onset_30.csv")[11,],1))
x$dist<-10*(0:10)
x$control_per_case<-x$Control/x$Case*80

library(ggplot2)
library(tidyr)
x<-pivot_longer(x,c("control_per_case","Case"))
x$name<-factor(x$name,levels=c("control_per_case","Case"),labels=c("Controls per case","Cases matched"))
p<-ggplot()+geom_point(data=x,aes(x=dist,y=value,color=name,shape=name))+scale_y_continuous(sec.axis=sec_axis(trans=~./80,name="Controls per case"))+
  labs(x="Maximum distance (km)",y="Cases matched",color="",shape="")+theme(legend.position="top")
ggsave("figs/matched_by_distance.png",p,width=3,height=3)
p

# Assign SIAs
for(i in dists){
  sias(define.distance=i,define.age=12,define.onset=30,define.offset=28)
  print(i)
}

# Do conditional regression
results<-lapply(c(1,(10000*(1:10))),function(i){regression(define.distance=i,define.age=12,define.onset=30,define.offset=28,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)

lims<-cbind.data.frame(x=12,y=c(min(r$lwr[r$var!="as.factor(doses_ipv012)2"]),max(r$upr[r$var!="as.factor(doses_ipv012)2"])))

if(any(grepl("ipv",r$var))){
  r$dist<-rep(c(0,(10*(1:10))),each=4)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"))
  p1<-ggplot()+geom_pointrange(data=r,aes(x=dist,y=med,ymin=lwr,ymax=upr),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum distance (km)",y="Effectiveness")+facet_wrap(~var,scales="free_y")
  
}else{
  r$dist<-rep(c(0,(10*(1:10))),each=2)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons"),labels=c("mOPV2","nOPV2"))
  p1<-ggplot()+geom_pointrange(data=r,aes(x=dist,y=med,ymin=lwr,ymax=upr,color=var),position=position_dodge(width=5),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum distance (km)",y="Effectiveness (per dose)",color="Vaccine")
  
}

p1<-p1+geom_point(data=lims,aes(x=x,y=y),alpha=0)
p1
ggsave("figs/distance_sa.png",p1,width=6,height=6)


# Sensitivity - matching age
ages<-seq(6,18,by=2)
ages<-ages[ages!=12]

# Match to controls
for(i in ages){
  matching(define.distance=50000,define.age=i,define.onset=30)
  print(i)
}

x<-sapply(seq(6,18,2),function(i){tail(read.csv(paste0("data/flow_match_distance_50000_age_",i,"_onset_30.csv"))[11,],1)})
x<-as.data.frame(t(as.matrix(x)))
x$Case<-unlist(x$Case)
x$Control<-unlist(x$Control)
x$age<-seq(6,18,2)
fac<-120
x$control_per_case<-x$Control/x$Case*fac

x<-pivot_longer(x,c("control_per_case","Case"))
x$name<-factor(x$name,levels=c("control_per_case","Case"),labels=c("Controls per case","Cases matched"))
p<-ggplot()+geom_point(data=x,aes(x=age,y=value,color=name,shape=name))+scale_y_continuous(sec.axis=sec_axis(trans=~./fac,name="Controls per case"))+
  labs(x="Maximum age (months)",y="Cases matched",color="",shape="")+theme(legend.position="top")
ggsave("figs/matched_by_age.png",p,width=3,height=3)
p

# Assign SIAs
for(i in ages){
  sias(define.distance=50000,define.age=i,define.onset=30,define.offset=28)
  print(i)
}

# Regression
results<-lapply(seq(6,18,2),function(i){regression(define.distance=50000,define.age=i,define.onset=30,define.offset=28,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)

lims<-cbind.data.frame(x=12,y=c(min(r$lwr[r$var!="as.factor(doses_ipv012)2"]),max(r$upr[r$var!="as.factor(doses_ipv012)2"])))

if(any(grepl("ipv",r$var))){
  r$age<-rep(seq(6,18,2),each=4)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"))
  p2<-ggplot()+geom_pointrange(data=r,aes(x=age,y=med,ymin=lwr,ymax=upr),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum age distance (months)",y="Effectiveness")+facet_wrap(~var,scales="free_y")
  
}else{
  r$age<-rep(seq(6,18,2),each=2)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons"),labels=c("mOPV2","nOPV2"))
  p2<-ggplot()+geom_pointrange(data=r,aes(x=age,y=med,ymin=lwr,ymax=upr,color=var),position=position_dodge(width=1),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum age distance (months)",y="Effectiveness (per dose)",color="Vaccine")
  
}

p2<-p2+geom_point(data=lims,aes(x=x,y=y),alpha=0)
p2
ggsave("figs/age_sa.png",p2,width=6,height=6)


# Sensitivity - onset
ons<-c(15,45,60,90,120)

# Matching
for(i in ons){
  matching(define.distance=50000,define.age=12,define.onset=i)
  print(i)
}

x<-sapply(c(15,30,45,60,90,120),function(i){tail(read.csv(paste0("data/flow_match_distance_50000_age_12_onset_",i,".csv"))[11,],1)})
x<-as.data.frame(t(as.matrix(x)))
x$Case<-unlist(x$Case)
x$Control<-unlist(x$Control)
x$onset<-c(15,30,45,60,90,120)
fac<-80
x$control_per_case<-x$Control/x$Case*fac

x<-pivot_longer(x,c("control_per_case","Case"))
x$name<-factor(x$name,levels=c("control_per_case","Case"),labels=c("Controls per case","Cases matched"))
p<-ggplot()+geom_point(data=x,aes(x=onset,y=value,color=name,shape=name))+scale_y_continuous(sec.axis=sec_axis(trans=~./fac,name="Controls per case"))+
  labs(x="Maximum onset distance (days)",y="Cases matched",color="",shape="")+theme(legend.position="top")
ggsave("figs/matched_by_onset.png",p,width=3,height=3)
p


# Assign SIAs
for(i in ons){
  sias(define.distance=50000,define.age=12,define.onset=i,define.offset=28)
  print(i)
}

# Regression
ons<-c(15,30,45,60,90,120)
results<-lapply(ons,function(i){regression(define.distance=50000,define.age=12,define.onset=i,define.offset=28,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)

lims<-cbind.data.frame(x=28,y=c(min(r$lwr[r$var!="as.factor(doses_ipv012)2"]),max(r$upr[r$var!="as.factor(doses_ipv012)2"])))
if(any(grepl("ipv",r$var))){
  r$onset<-rep(ons,each=4)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"))
  p3<-ggplot()+geom_point(data=lims,aes(x=x,y=y),alpha=0)+geom_pointrange(data=r,aes(x=onset,y=med,ymin=lwr,ymax=upr),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum onset distance (days)",y="Effectiveness")+facet_wrap(~var,scales="free_y")
  
}else{
  r$onset<-rep(ons,each=2)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons"),labels=c("mOPV2","nOPV2"))
  p3<-ggplot()+geom_point(data=lims,aes(x=x,y=y),alpha=0)+geom_pointrange(data=r,aes(x=onset,y=med,ymin=lwr,ymax=upr,color=var),position=position_dodge(width=3),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum onset distance (days)",y="Effectiveness (per dose)",color="Vaccine")
  
}

p3
ggsave("figs/onset_sa.png",p3,width=6,height=6)


# Sensitivity - offset
# Assign SIAs
offs<-seq(0,8*7,by=7)
offs<-offs[offs!=28]
for(i in offs){
  sias(define.distance=50000,define.age=12,define.onset=30,define.offset=i)
  print(i)
}

# Regression
results<-lapply(seq(0,8*7,by=7),function(i){regression(define.distance=50000,define.age=12,define.onset=30,define.offset=i,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)

if(any(grepl("ipv",r$var))){
  r$offset<-rep(seq(0,8*7,by=7),each=4)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2","nOPV2","IPV (1 dose)","IPV (2 doses)"))
  p4<-ggplot()+geom_pointrange(data=r[grepl("OPV",r$var),],aes(x=offset,y=med,ymin=lwr,ymax=upr,color=var),position=position_dodge(width=5),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Offset (days)",y="Effectiveness (per dose)",color="Vaccine")
  
}else{
  r$offset<-rep(seq(0,8*7,by=7),each=2)
  r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2","nOPV2","IPV (1 dose)","IPV (2 doses)"))
  p4<-ggplot()+geom_pointrange(data=r[grepl("OPV",r$var),],aes(x=offset,y=med,ymin=lwr,ymax=upr,color=var),position=position_dodge(width=5),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Offset (days)",y="Effectiveness (per dose)",color="Vaccine")
  
}


p4
ggsave("figs/offset_sa.png",p4,width=5,height=3)



# Sensitivity - plot for all matching parameters
lims<-cbind.data.frame(x=12,var=rep(c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"),each=2),y=c(0,.6,0,.6,0,.6,-2.5,1))

results<-lapply(c(1,(10000*(1:10))),function(i){regression(define.distance=i,define.age=12,define.onset=30,define.offset=28,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)


r$dist<-rep(c(0,(10*(1:10))),each=4)
r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"))
p1<-ggplot()+geom_hline(data=r%>%filter(dist==50),aes(yintercept=med),linetype=1)+
  geom_hline(data=r%>%filter(dist==50),aes(yintercept=lwr),linetype=2)+
  geom_hline(data=r%>%filter(dist==50),aes(yintercept=upr),linetype=2)+
  geom_pointrange(data=r,aes(x=dist,y=med,ymin=lwr,ymax=upr),fatten=1)+
  scale_y_continuous(labels=scales::percent)+
  labs(x="Maximum distance (km)",y="Effectiveness")+
  facet_wrap(~var,scales="free_y",ncol=1)


p1<-p1+geom_point(data=lims,aes(x=x,y=y),alpha=0)
p1



# Regression
results<-lapply(seq(6,18,2),function(i){regression(define.distance=50000,define.age=i,define.onset=30,define.offset=28,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)

r$age<-rep(seq(6,18,2),each=4)
r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"))
p2<-ggplot()+geom_hline(data=r%>%filter(age==12),aes(yintercept=med),linetype=1)+
  geom_hline(data=r%>%filter(age==12),aes(yintercept=lwr),linetype=2)+
  geom_hline(data=r%>%filter(age==12),aes(yintercept=upr),linetype=2)+
  geom_pointrange(data=r,aes(x=age,y=med,ymin=lwr,ymax=upr),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum age distance (months)",y="")+facet_wrap(~var,scales="free_y",ncol=1)


p2<-p2+geom_point(data=lims,aes(x=x,y=y),alpha=0)
p2


# Regression
ons<-c(15,30,45,60,90,120)
results<-lapply(ons,function(i){regression(define.distance=50000,define.age=12,define.onset=i,define.offset=28,IPV=T)})
r<-NULL
for(i in 1:length(results)){
  r[[i]]<-summary(results[[i]]$model)$conf.int
  r[[i]]<-cbind(r[[i]],var=row.names(r[[i]]))
}
r<-do.call("rbind",r)
r<-as.data.frame(r)
for(j in 1:ncol(r)){r[,j]<-unlist(r[,j])}
r$lwr<-1-as.numeric(r$`upper .95`)
r$upr<-1-as.numeric(r$`lower .95`)
r$med<-1-as.numeric(r$`exp(coef)`)

r$onset<-rep(ons,each=4)
r$var<-factor(r$var,c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2 (per dose)","nOPV2 (per dose)","IPV (1 dose)","IPV (2 doses)"))
p3<-ggplot()+geom_hline(data=r%>%filter(onset==30),aes(yintercept=med),linetype=1)+
  geom_hline(data=r%>%filter(onset==30),aes(yintercept=lwr),linetype=2)+
  geom_hline(data=r%>%filter(onset==30),aes(yintercept=upr),linetype=2)+geom_point(data=lims,aes(x=x,y=y),alpha=0)+geom_pointrange(data=r,aes(x=onset,y=med,ymin=lwr,ymax=upr),fatten=1)+scale_y_continuous(labels=scales::percent)+labs(x="Maximum onset distance (days)",y="")+facet_wrap(~var,scales="free_y",ncol=1)

p3


p<-ggarrange(p1,p2,p3,nrow=1,ncol=3,common.legend=T)
p
ggsave("figs/distance_age_onset_sa.png",p,width=9,height=7)


