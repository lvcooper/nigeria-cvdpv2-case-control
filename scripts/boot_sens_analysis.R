source("scripts/sias_boot.R")
# This should take < 5 minutes
sias_boot(define.distance=50000,define.age=12,define.onset=30,define.offset=28,match.on.RI=F,use.date.of.last.opv=T,boot=1000)

source("scripts/regression_boot.R")
x<-regression_boot(define.distance=50000,define.age=12,define.onset=30,define.offset=28,match.on.RI=F,use.date.of.last.opv=T,boot=1000,IPV=T)
write.csv(x,"data/regression_boot_no_cutoff.csv")
x<-read.csv("data/regression_boot_no_cutoff.csv")
sub<-x[x$var=="nOPV2",]
x$boot_rank<-as.numeric(factor(x$boot, levels=sub$boot[order(sub$OR)],ordered=T))
ggplot()+geom_pointrange(data=x%>%filter(boot_rank%%25==0&grepl("OPV",var)),aes(x=boot_rank,y=eff,ymin=eff_lwr,ymax=eff_upr,color=var),position=position_dodge(width=10),fatten=1)

mutate(group_by(x,var),select=(eff_lwr==min(eff_lwr)))%>%filter(select)
mutate(group_by(x,var),select=(eff_upr==max(eff_upr)))%>%filter(select)

# Slight inverse correlation between nOPV2 effectiveness and mOPV2 effectiveness estimates
ggplot(data=pivot_wider(x,id_cols=boot_rank,names_from=var,values_from=eff),aes(x=mOPV2,y=nOPV2))+geom_point()+geom_smooth(method="lm")

source("scripts/regression.R")
y<-regression(define.distance=50000,define.age=12,define.onset=30,define.offset=28,use.date.of.last.opv=T,IPV=T)
tab<-function(mod){
  x<-cbind(summary(mod)$coef[,2],
           summary(mod)$conf[,3:4])
  x<-cbind(x,
           100*(1-x))
  x<-as.data.frame(x)
  names(x)<-c("OR","OR_lwr","OR_upr","eff","eff_upr","eff_lwr")
  x$var<-factor(rownames(x),c( "mOPV2_ons" ,"nOPV2_ons", "as.factor(doses_ipv012)1","as.factor(doses_ipv012)2"),labels=c("mOPV2","nOPV2","IPV (1 dose)","IPV (2 doses)"))
  return(x)
}
bind_rows(tab(y$model),
          summarise_all(group_by(x,var),median))%>%ggplot()+geom_pointrange(aes(x=var,y=eff,ymin=eff_lwr,ymax=eff_upr,color=as.factor(boot)),position=position_dodge(width=0.5),fatten=1)

ggplot(data=pivot_longer(x,starts_with("eff")))+geom_histogram(aes(x=value,fill=var))+
  geom_vline(data=pivot_longer(tab(y$model),starts_with("eff")),aes(xintercept=value))+facet_wrap(name~var,scales="free")


ggplot(data=pivot_longer(x%>%filter(grepl("OPV",var)),starts_with("eff"),names_transform=function(x){ifelse(x=="eff_lwr","Lower",ifelse(x=="eff_upr","Upper","Central"))}))+
  geom_histogram(aes(x=value,color=name,fill=name),alpha=0.5,binwidth=1,position=position_dodge(width=0.1))+labs(x="Effectiveness (%)",y="Frequency",color="Estimate",fill="Estimate")+
  geom_vline(data=pivot_longer(tab(y$model)%>%filter(grepl("OPV",var)),starts_with("eff"),names_transform=function(x){ifelse(x=="eff_lwr","Lower",ifelse(x=="eff_upr","Upper","Central"))}),
             aes(xintercept=value,color=name),linetype=2)+facet_wrap(~var,scales="free_x")

ggplot(data=pivot_longer(x,starts_with("eff"),names_transform=function(x){ifelse(x=="eff_lwr","Lower",ifelse(x=="eff_upr","Upper","Central"))}))+
  geom_histogram(aes(x=value,color=name,fill=name),alpha=0.5,binwidth=1,position=position_dodge(width=0.1))+labs(x="Effectiveness (%)",y="Frequency",color="Estimate",fill="Estimate")+
  geom_vline(data=pivot_longer(tab(y$model),starts_with("eff"),names_transform=function(x){ifelse(x=="eff_lwr","Lower",ifelse(x=="eff_upr","Upper","Central"))}),aes(xintercept=value,color=name),linetype=2)+facet_wrap(~var,scales="free_x")
ggsave("figs/nigeria_case_control_bootstrap.png",width=10/1.5,height=7/1.5)
summarise_all(ungroup(summarise(group_by(x,boot),mopv2_greater_nopv2=eff_lwr[var=="mOPV2"]>eff_upr[var=="nOPV2"],
                                nopv2_greater_mopv2=eff_lwr[var=="nOPV2"]>eff_upr[var=="mOPV2"])),sum)
# If you bootstrap over possible doses, 0 out of 1000 finds a significant difference between nOPV2 and mOPV2




# Non-linear
x<-regression_boot_nonlinear(define.distance=50000,define.age=12,define.onset=30,define.offset=28,match.on.RI=F,use.date.of.last.opv=T,boot=1000,grps=6,IPV=T)
write.csv(x,"data/regression_boot_nonlinear.csv")
x<-read.csv("data/regression_boot_nonlinear.csv")
summarise_all(group_by(x,vac,doses),median)%>%filter(grepl("OPV",vac))%>%ggplot()+geom_pointrange(aes(x=doses,y=eff,ymin=ifelse(eff_lwr<0,0,eff_lwr),ymax=eff_upr),fatten=1)+facet_wrap(~vac)

y<-regression_boot(define.distance=50000,define.age=12,define.onset=30,define.offset=28,match.on.RI=F,use.date.of.last.opv=T,boot=1000,cutoff=4,IPV=T)
write.csv(y,"data/regression_boot_cutoff_4.csv")
y<-read.csv("data/regression_boot_cutoff_4.csv")

z<-read.csv("data/regression_boot_no_cutoff.csv")
points<-summarise_all(group_by(x,vac,doses),median)%>%filter(grepl("OPV",vac))
z<-bind_rows(cbind(z,type="All"),cbind(y,type="<4 OPV2"))
lines<-summarise(group_by(z,var,type),or=log(median(OR)),lwr=log(median(OR_lwr)),upr=log(median(OR_upr)))%>%filter(grepl("OPV",var))
lines<-merge(lines,cbind.data.frame(doses=1:max(points$doses)))
lines$OR<-exp(lines$or*lines$doses)
lines$OR_lwr<-exp(lines$lwr*lines$doses)
lines$OR_upr<-exp(lines$upr*lines$doses)
lines$vac<-lines$var
lines<-filter(lines, type=="All"|doses<4)
ggplot()+geom_ribbon(data=lines,aes(x=doses,ymin=OR_lwr,ymax=OR_upr,fill=type),alpha=0.3)+
  geom_pointrange(data=points,aes(x=doses,y=OR,ymin=OR_lwr,ymax=OR_upr),fatten=1)+
  geom_line(data=lines,aes(x=doses,y=OR,color=type))+facet_wrap(~vac)+scale_y_log10()+
  labs(x="Doses",y="Odds ratio",color="",fill="")

summarise_all(group_by(y,var),median)
summarise_all(group_by(z,var),median)

define.distance=50000;define.age=12;define.onset=30;match.on.RI=F;define.offset=28;use.date.of.last.opv=T;boot=1000
filename<-paste0("data/matched_distance_",define.distance,
                 "_age_",define.age,
                 "_onset_",define.onset,
                 ifelse(match.on.RI,"_RI_",""),
                 "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),"_boot_",boot,".csv")
m1<-read.csv(filename)
summarise(group_by(m1,case,boot),m=sum(mOPV2_ons<4&nOPV2_ons<4),n=n())%>%group_by(case)%>%summarise(m=mean(m),n=mean(n))
table(m1$mOPV2_ons)
table(m1$nOPV2_ons)
