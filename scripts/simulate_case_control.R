# Number of times to repeat simulation
# You might want to set this lower for testing
N<-500

# Get reasonable parameter estimates
m1<-read.csv("data/matched_distance_50000_age_12_onset_30_offset_28.csv")

mean(m1$mOPV2_SIAs_ons)
mean(m1$nOPV2_SIAs_ons)
mean(m1$OPV_SIAs_ons-m1$mOPV2_SIAs_ons-m1$nOPV2_SIAs_ons)

median(table(m1$stratum))

sum(m1$case)

mean(m1$doses_opv_sia/m1$OPV_SIAs,na.rm=T)

# Functions to simulate case control analysis
simulate_case_control1<-function(N_case,n_per_strata,coverage,recall_variation,zero_mean,eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,case_control_ratio){
  # eff_nOPV2=0.8;eff_mOPV2=0.8;recall_variation=0.2;N_case=350;coverage=0.75;mean_nOPV2=1.1;mean_mOPV2=1.5;mean_bOPV=4.9;case_control_ratio=4;zero_mean=0.2;n_per_strata=100
  
  # # Number of SIAs they were exposed to (identical within strata)
  # nOPV2_SIAs=rep(rpois(n=N_case,lambda=mean_nOPV2),each=n_per_strata)
  # mOPV2_SIAs=rep(rpois(n=N_case,lambda=mean_mOPV2),each=n_per_strata)
  # bOPV_SIAs=rep(rpois(n=N_case,lambda=mean_bOPV),each=n_per_strata)
  
  # Number of SIAs they were exposed to (different within strata)
  nOPV2_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_nOPV2)
  mOPV2_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_mOPV2)
  bOPV_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_bOPV)
  
  # Number of doses they received
  nOPV2_doses=rbinom(n=length(nOPV2_SIAs),size=nOPV2_SIAs,prob=coverage)
  mOPV2_doses=rbinom(n=length(mOPV2_SIAs),size=mOPV2_SIAs,prob=coverage)
  bOPV_doses=rbinom(n=length(bOPV_SIAs),size=bOPV_SIAs,prob=coverage)
  
  # Number of total OPV SIA doses
  OPV_doses=nOPV2_doses+mOPV2_doses+bOPV_doses
  
  # Strata for matched case-control (zero mean normally distributed)
  strata=rep(1:(N_case),each=n_per_strata)
  # Log odds from strata
  strata_value_log_OR=rep(rnorm(n=max(strata),mean=0,sd=1),each=n_per_strata)
  tab<-cbind.data.frame(nOPV2_SIAs,nOPV2_doses,
                        mOPV2_SIAs,mOPV2_doses,
                        bOPV_SIAs,bOPV_doses,
                        OPV_doses,strata,strata_value_log_OR)
  
  # Number of total doses reported
  # matplot(c(0,1:10),t(sapply(c(0.4,1:10),function(i){qlnorm(p=c(0.025,0.5,0.975),meanlog=log(i),sdlog=0.2)})),type="l")
  tab<-tab[order(tab$OPV_doses),]
  a<-as.numeric(names(table(tab$OPV_doses)))
  b<-as.numeric(table(tab$OPV_doses))
  if(recall_variation>0) tab$OPV_doses_reported=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=recall_variation))}))
  if(recall_variation==0) tab$OPV_doses_reported=tab$OPV_doses
  
  # Log odds from effectiveness
  log_OR_nOPV2=log(1-eff_nOPV2)
  
  # Log odds * doses
  tab$nOPV2_value_log_OR=log_OR_nOPV2*tab$nOPV2_doses
  
  # Log odds from effectiveness
  log_OR_mOPV2=log(1-eff_mOPV2)
  
  # Log odds * doses
  tab$mOPV2_value_log_OR=log_OR_mOPV2*tab$mOPV2_doses
  
  
  # Full regression formula (logit probability)
  tab$logit_p=tab$nOPV2_value_log_OR+tab$mOPV2_value_log_OR+tab$strata_value_log_OR
  
  # Case or control based on probability
  tab$case=rbinom(n=length(tab$logit_p),size=1,prob=inv.logit(tab$logit_p))
  
  # Select strata with at least 1 case and at least case:control ratio controls
  tab<-mutate(group_by(as_tibble(tab),strata),n_case=sum(case),n_control=sum(case==0))
  sub<-tab%>%filter(n_case>=1,n_control>=case_control_ratio)
  
  # Select out N_case cases and case_control_ratio controls
  # Cases
  sub_case<-sub[sub$case==1,]
  # Order randomly
  sub_case<-sub_case[sample(x=1:nrow(sub_case),size=nrow(sub_case),replace=F),]
  # Case count per strata
  sub_case<-mutate(group_by(as_tibble(sub_case),strata),count=1:n())
  # Select first case per strata
  sel<-which(sub_case$count==1)
  # Select randomly
  if(length(sel)>N_case) sel<-sel[sample(x=1:length(sel),size=N_case,replace=F)]
  sub_case<-sub_case[sel,]
  # Controls
  sub_control<-sub[sub$case==0&sub$strata%in%sub_case$strata,]
  # Order randomly
  sub_control<-sub_control[sample(x=1:nrow(sub_control),size=nrow(sub_control),replace=F),]
  # Control count per strata
  sub_control<-mutate(group_by(as_tibble(sub_control),strata),count=1:n())
  # Select first case_control_ratio controls per strata
  sel<-which(sub_control$count<=case_control_ratio)
  sub_control<-sub_control[sel,]
  # Combine
  data<-as.data.frame(bind_rows(sub_case,sub_control))
  
  # Assign valency of doses
  data$OPV_SIAs<-rowSums(data[,c("bOPV_SIAs","nOPV2_SIAs","mOPV2_SIAs")])
  data$nOPV2_doses_frac<-data$OPV_doses*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac<-data$OPV_doses*data$mOPV2_SIAs/data$OPV_SIAs
  data$nOPV2_doses_frac_recall<-data$OPV_doses_reported*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall<-data$OPV_doses_reported*data$mOPV2_SIAs/data$OPV_SIAs
  
  # Regression
  modA<-clogit(case~mOPV2_doses+nOPV2_doses+strata(strata),data=data)
  modB<-clogit(case~mOPV2_doses_frac+nOPV2_doses_frac+strata(strata),data=data)
  modC<-clogit(case~mOPV2_doses_frac_recall+nOPV2_doses_frac_recall+strata(strata),data=data)
  
  p<-sapply(list(modA,modB,modC),function(m){x<-summary(m)$coef;2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))})
  
  fn<-function(mod){
    m<-summary(mod)$coef[,c(1,3)]
    r<-as.data.frame(1-exp(t(sapply(1:2,function(i){m[i,1]+c(0,-1,1)*1.96*m[i,2]}))))
    names(r)<-c("eff","upr","lwr")
    return(r)
  }
  
  r<-cbind.data.frame(type=rep(c("True","Frac","Frac + recall"),each=2),vaccine=rep(c("mOPV2","nOPV2"),times=3),
                      rbind(fn(modA),fn(modB),fn(modC)),pvalue=rep(p,each=2))
  ret<-cbind.data.frame(r,N_case=N_case,n_per_strata=n_per_strata,coverage=coverage,recall_variation=recall_variation,
                        zero_mean=zero_mean,
                        case_control_ratio=case_control_ratio,
                        eff_nOPV2=eff_nOPV2,eff_mOPV2=eff_mOPV2,mean_nOPV2=mean_nOPV2,mean_mOPV2=mean_mOPV2,mean_bOPV=mean_bOPV,
                        check_case=sum(data$case),check_control=sum(data$case==0),inacc=sum(data$OPV_doses_reported!=data$OPV_doses))
  return(ret)}


simulate_case_control2<-function(N_case,n_per_strata,zero_mean,coverage,eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,case_control_ratio){
  # eff_nOPV2=0.8;eff_mOPV2=0.8;N_case=350;coverage=0.75;mean_nOPV2=1.1;mean_mOPV2=1.5;mean_bOPV=4.9;case_control_ratio=4;n_per_strata=100;zero_mean=0.2
  
  # # Number of SIAs they were exposed to (identical within strata)
  # nOPV2_SIAs=rep(rpois(n=N_case,lambda=mean_nOPV2),each=n_per_strata)
  # mOPV2_SIAs=rep(rpois(n=N_case,lambda=mean_mOPV2),each=n_per_strata)
  # bOPV_SIAs=rep(rpois(n=N_case,lambda=mean_bOPV),each=n_per_strata)
  
  # Number of SIAs they were exposed to (different within strata)
  nOPV2_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_nOPV2)
  mOPV2_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_mOPV2)
  bOPV_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_bOPV)
  
  # Number of doses they received
  nOPV2_doses=rbinom(n=length(nOPV2_SIAs),size=nOPV2_SIAs,prob=coverage)
  mOPV2_doses=rbinom(n=length(mOPV2_SIAs),size=mOPV2_SIAs,prob=coverage)
  bOPV_doses=rbinom(n=length(bOPV_SIAs),size=bOPV_SIAs,prob=coverage)
  
  # Number of total OPV SIA doses
  OPV_doses=nOPV2_doses+mOPV2_doses+bOPV_doses
  
  # Strata for matched case-control (zero mean normally distributed)
  strata=rep(1:(N_case),each=n_per_strata)
  # Log odds from strata
  strata_value_log_OR=rep(rnorm(n=max(strata),mean=0,sd=1),each=n_per_strata)
  tab<-cbind.data.frame(nOPV2_SIAs,nOPV2_doses,
                        mOPV2_SIAs,mOPV2_doses,
                        bOPV_SIAs,bOPV_doses,
                        OPV_doses,strata,strata_value_log_OR)
  
  # Number of total doses reported
  tab<-tab[order(tab$OPV_doses),]
  a<-as.numeric(names(table(tab$OPV_doses)))
  b<-as.numeric(table(tab$OPV_doses))
  tab$OPV_doses_reported_0.1=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=0.1))}))
  tab$OPV_doses_reported_0.2=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=0.2))}))
  tab$OPV_doses_reported_0.3=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=0.3))}))
  tab$OPV_doses_reported_0.4=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=0.4))}))
  tab$OPV_doses_reported_0.5=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=0.5))}))
  
  # Log odds from effectiveness
  log_OR_nOPV2=log(1-eff_nOPV2)
  
  # Log odds * doses
  tab$nOPV2_value_log_OR=log_OR_nOPV2*tab$nOPV2_doses
  
  # Log odds from effectiveness
  log_OR_mOPV2=log(1-eff_mOPV2)
  
  # Log odds * doses
  tab$mOPV2_value_log_OR=log_OR_mOPV2*tab$mOPV2_doses
  
  
  # Full regression formula (logit probability)
  tab$logit_p=tab$nOPV2_value_log_OR+tab$mOPV2_value_log_OR+tab$strata_value_log_OR
  
  # Case or control based on probability
  tab$case=rbinom(n=length(tab$logit_p),size=1,prob=inv.logit(tab$logit_p))
  
  # Select strata with at least 1 case and at least case:control ratio controls
  tab<-mutate(group_by(as_tibble(tab),strata),n_case=sum(case),n_control=sum(case==0))
  sub<-tab%>%filter(n_case>=1,n_control>=case_control_ratio)
  
  # Select out N_case cases and case_control_ratio controls
  # Cases
  sub_case<-sub[sub$case==1,]
  # Order randomly
  sub_case<-sub_case[sample(x=1:nrow(sub_case),size=nrow(sub_case),replace=F),]
  # Case count per strata
  sub_case<-mutate(group_by(as_tibble(sub_case),strata),count=1:n())
  # Select first case per strata
  sel<-which(sub_case$count==1)
  # Select randomly
  if(length(sel)>N_case) sel<-sel[sample(x=1:length(sel),size=N_case,replace=F)]
  sub_case<-sub_case[sel,]
  # Controls
  sub_control<-sub[sub$case==0&sub$strata%in%sub_case$strata,]
  # Order randomly
  sub_control<-sub_control[sample(x=1:nrow(sub_control),size=nrow(sub_control),replace=F),]
  # Control count per strata
  sub_control<-mutate(group_by(as_tibble(sub_control),strata),count=1:n())
  # Select first case_control_ratio controls per strata
  sel<-which(sub_control$count<=case_control_ratio)
  sub_control<-sub_control[sel,]
  # Combine
  data<-as.data.frame(bind_rows(sub_case,sub_control))
  
  # Assign valency of doses
  data$OPV_SIAs<-rowSums(data[,c("bOPV_SIAs","nOPV2_SIAs","mOPV2_SIAs")])
  data$nOPV2_doses_frac<-data$OPV_doses*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac<-data$OPV_doses*data$mOPV2_SIAs/data$OPV_SIAs
  data$nOPV2_doses_frac_recall_0.1<-data$OPV_doses_reported_0.1*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall_0.1<-data$OPV_doses_reported_0.1*data$mOPV2_SIAs/data$OPV_SIAs
  data$nOPV2_doses_frac_recall_0.2<-data$OPV_doses_reported_0.2*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall_0.2<-data$OPV_doses_reported_0.2*data$mOPV2_SIAs/data$OPV_SIAs
  data$nOPV2_doses_frac_recall_0.3<-data$OPV_doses_reported_0.3*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall_0.3<-data$OPV_doses_reported_0.3*data$mOPV2_SIAs/data$OPV_SIAs
  data$nOPV2_doses_frac_recall_0.4<-data$OPV_doses_reported_0.4*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall_0.4<-data$OPV_doses_reported_0.4*data$mOPV2_SIAs/data$OPV_SIAs
  data$nOPV2_doses_frac_recall_0.5<-data$OPV_doses_reported_0.5*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall_0.5<-data$OPV_doses_reported_0.5*data$mOPV2_SIAs/data$OPV_SIAs
  
  # Regression
  modA<-clogit(case~mOPV2_doses+nOPV2_doses+strata(strata),data=data)
  modB<-clogit(case~mOPV2_doses_frac+nOPV2_doses_frac+strata(strata),data=data)
  modC<-clogit(case~mOPV2_doses_frac_recall_0.1+nOPV2_doses_frac_recall_0.1+strata(strata),data=data)
  modD<-clogit(case~mOPV2_doses_frac_recall_0.2+nOPV2_doses_frac_recall_0.2+strata(strata),data=data)
  modE<-clogit(case~mOPV2_doses_frac_recall_0.3+nOPV2_doses_frac_recall_0.3+strata(strata),data=data)
  modF<-clogit(case~mOPV2_doses_frac_recall_0.4+nOPV2_doses_frac_recall_0.4+strata(strata),data=data)
  modG<-clogit(case~mOPV2_doses_frac_recall_0.5+nOPV2_doses_frac_recall_0.5+strata(strata),data=data)
  
  p<-sapply(list(modA,modB,modC,modD,modE,modF,modG),function(m){x<-summary(m)$coef;2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))})
  
  fn<-function(mod){
    m<-summary(mod)$coef[,c(1,3)]
    r<-as.data.frame(1-exp(t(sapply(1:2,function(i){m[i,1]+c(0,-1,1)*1.96*m[i,2]}))))
    names(r)<-c("eff","upr","lwr")
    return(r)
  }
  
  r<-cbind.data.frame(type=rep(c("True","Frac","Frac + recall","Frac + recall","Frac + recall","Frac + recall","Frac + recall"),each=2),
                      vaccine=rep(c("mOPV2","nOPV2"),times=7),recall_variation=rep(c(0,0,0.1,0.2,0.3,0.4,0.5),each=2),
                      rbind(fn(modA),fn(modB),fn(modC),fn(modD),fn(modE),fn(modF),fn(modG)),pvalue=rep(p,each=2))
  ret<-cbind.data.frame(r,N_case=N_case,zero_mean=zero_mean,n_per_strata=n_per_strata,coverage=coverage,case_control_ratio=case_control_ratio,
                        eff_nOPV2=eff_nOPV2,eff_mOPV2=eff_mOPV2,mean_nOPV2=mean_nOPV2,mean_mOPV2=mean_mOPV2,mean_bOPV=mean_bOPV,
                        check_case=sum(data$case),check_control=sum(data$case==0))
  return(ret)}

frm<-expand.grid(boot=1:N,eff_nOPV2=seq(0.2,0.8,0.2),eff_mOPV2=seq(0.2,0.8,0.2),n_per_strata=1000,zero_mean=0.2,
                 N_case=350,coverage=0.75,mean_nOPV2=1.5,mean_mOPV2=1.1,mean_bOPV=4.9,case_control_ratio=4)

res_list=NULL
start<-Sys.time()
for(i in 1:nrow(frm)){
  res_list[[i]]<-cbind.data.frame(boot=frm$boot[i],
                                  simulate_case_control2(N_case=frm$N_case[i],coverage=frm$coverage[i],zero_mean=frm$zero_mean[i],
                                                         case_control_ratio=frm$case_control_ratio[i],n_per_strata=frm$n_per_strata[i],
                                                         eff_nOPV2=frm$eff_nOPV2[i],eff_mOPV2=frm$eff_mOPV2[i],
                                                         mean_nOPV2=frm$mean_nOPV2[i],mean_mOPV2=frm$mean_mOPV2[i],mean_bOPV=frm$mean_bOPV[i]))
  print(i)
}
Sys.time()-start
res<-do.call("rbind",res_list)
write.csv(res,"data/simulate_matched_case_control_Ncase350_boot.csv")



# Visualise/analyse simulated results
res<-read.csv("data/simulate_matched_case_control_Ncase350_boot100.csv")

summary(res$check_case)

res<-summarise(group_by(res,N_case,check_case,check_control,recall_variation,n_per_strata,zero_mean,coverage,case_control_ratio,
                        eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,type,vaccine),
               eff=median(eff),lwr=median(lwr),upr=median(upr),signif_diff=sum(pvalue<0.01),n_boot=n())

res$type<-as.character(factor(res$type,levels=c("Frac","Frac + recall","True"),labels=c("Fractional","Inaccurate recall, fractional","True")))
res$type_spec<-paste0(res$type,ifelse(res$type=="Inaccurate recall, fractional",paste0(" (SD = ",res$recall_variation,")"),""))
res$type_spec<-factor(res$type_spec,levels=c("True","Fractional",sort(unique(res$type_spec[grepl("SD",res$type_spec)]))),ordered=T)
res$eff_mOPV2_char<-paste0("True mOPV2 eff.: ",round(100*res$eff_mOPV2),"%")
ggplot(data=res%>%filter(vaccine=="nOPV2"&(recall_variation==0.2|type%in%c("True","Fractional"))))+geom_abline()+
  geom_pointrange(aes(x=eff_nOPV2,y=eff,ymin=lwr,ymax=upr,color=type),position=position_dodge(width=0.05),fatten=1)+
  facet_wrap(~eff_mOPV2_char)+labs(x="True nOPV2 effectiveness",y="Estimated nOPV2 effectiveness",color="History")+
  scale_y_continuous(labels=scales::percent)+scale_x_continuous(labels=scales::percent)+theme(legend.position="top")
ggsave("figs/simulate_matched_case_control_Ncase350.png",width=5,height=5)

res$eff_nOPV2_char<-paste0(round(100*res$eff_nOPV2),"%")
res$eff_mOPV2_char<-paste0(round(100*res$eff_mOPV2),"%")
ggplot(data=res%>%filter(vaccine=="nOPV2"&type!="True"))+geom_hline(aes(yintercept=eff_nOPV2))+
  scale_y_continuous(labels=scales::percent,breaks=seq(-0.2,0.8,0.2),sec.axis=dup_axis(name="True mOPV2 effectiveness",labels=NULL,breaks=NULL))+
  scale_x_continuous(breaks=c(0.1,0.3,0.5),sec.axis=dup_axis(name="True nOPV2 effectiveness",labels=NULL,breaks=NULL))+
  geom_pointrange(aes(x=recall_variation,y=eff,ymin=lwr,ymax=upr),position=position_dodge(width=0.01),fatten=1)+facet_grid(eff_mOPV2_char~eff_nOPV2_char)+labs(x="Recall variation (SD)",y="Estimated nOPV2 effectiveness")
ggsave("figs/simulate_matched_case_control_recall_Ncase350.png",width=5,height=5)

# How high would recall bias need to be for seroconversion estimates (~50% in Ochoge [the Gambia] after 1 dose in infants) to be consistent with case-control study (~15%)?
frm<-expand.grid(boot=1:N,eff_nOPV2=0.5,n_per_strata=1000,zero_mean=0.2,recall_variation=seq(0,1,by=0.2),
                 N_case=350,coverage=0.75,mean_nOPV2=1.5,mean_mOPV2=1.1,mean_bOPV=4.9,case_control_ratio=4)
frm$eff_mOPV2<-frm$eff_nOPV2
res_list=NULL
start<-Sys.time()
for(i in 1:nrow(frm)){
  res_list[[i]]<-cbind.data.frame(boot=frm$boot[i],
                                  simulate_case_control1(N_case=frm$N_case[i],coverage=frm$coverage[i],case_control_ratio=frm$case_control_ratio[i],
                                                         recall_variation=frm$recall_variation[i],
                                                         eff_nOPV2=frm$eff_nOPV2[i],n_per_strata=frm$n_per_strata[i],zero_mean=frm$zero_mean[i],
                                                         eff_mOPV2=frm$eff_mOPV2[i],mean_nOPV2=frm$mean_nOPV2[i],
                                                         mean_mOPV2=frm$mean_mOPV2[i],mean_bOPV=frm$mean_bOPV[i]))
  print(i)
}
Sys.time()-start
res<-do.call("rbind",res_list)
res$date<-Sys.time()
write.csv(res,"data/min_recall_var.csv",row.names=F)

res<-read.csv("data/min_recall_var.csv")
summarise(group_by(res,N_case,n_per_strata,check_case,zero_mean,check_control,recall_variation,coverage,case_control_ratio,
                   eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,type,vaccine),
          under=sum(lwr<0.15),eff=median(eff),lwr=median(lwr),upr=median(upr),n_boot=n(),inacc=sum(inacc))%>%
  filter(vaccine=='nOPV2',type=="Frac + recall")%>%
  ggplot()+geom_pointrange(aes(x=recall_variation,y=eff,ymin=lwr,ymax=upr),fatten=1)+
  scale_y_continuous(labels=scales::percent,sec.axis=dup_axis(name="Proportion of cases and controls with inaccurate reporting (%)"))+
  geom_line(aes(x=recall_variation,y=inacc/(n_boot*(check_case+check_control))),linetype=2,color="red")+
  geom_hline(yintercept=0.5)+geom_hline(yintercept=0.15,linetype=2)+labs(x="Recall variation (SD of discrete log-normal)",y="Estimated vaccine effectiveness (%)")
ggsave("figs/min_recall_var.png",width=5,height=5)

# Visualize recall bias
d<-rbind.data.frame(cbind(0.25,0:10,t(sapply(c(0.1,1:10),function(i){qlnorm(p=c(0.025,0.5,0.975),meanlog=log(i),sdlog=0.25)}))),
                    cbind(0.5,0:10,t(sapply(c(0.1,1:10),function(i){qlnorm(p=c(0.025,0.5,0.975),meanlog=log(i),sdlog=0.5)}))),
                    cbind(0.75,0:10,t(sapply(c(0.1,1:10),function(i){qlnorm(p=c(0.025,0.5,0.975),meanlog=log(i),sdlog=0.75)}))),
                    cbind(1,0:10,t(sapply(c(0.1,1:10),function(i){qlnorm(p=c(0.025,0.5,0.975),meanlog=log(i),sdlog=1)}))))
names(d)<-c("recall_variation","true_doses","lwr","med","upr")
d$recall_variation<-factor(d$recall_variation,levels=seq(1,0.25,-0.25),ordered = T)
ggplot()+
  geom_line(data=d,aes(x=true_doses,y=upr,color=recall_variation))+
  geom_ribbon(data=d,aes(x=true_doses,ymin=lwr,ymax=upr,fill=recall_variation))+
  geom_line(data=d,aes(x=true_doses,y=med))+labs(x="True doses",y="Reported doses",color="SD",fill="SD")+scale_x_continuous(breaks=seq(0,10,2))
ggsave("figs/demo_recall.png",width=5,height=5)

# At what point can we detect significant differences?
# Observed effectiveness of ~15% for nOPV2
# Option 1: true eff = 0.15 and SD = 0
# Option 2: true eff = 0.5 and SD = 0.8 (as above)
# What is the minimum detectable difference for each of these?

# Simulation function
simulate_case_control3<-function(N_case,n_per_strata,coverage,recall_variation,zero_mean,eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,case_control_ratio){
  # eff_nOPV2=0.8;eff_mOPV2=0.8;recall_variation=0.2;N_case=350;coverage=0.75;mean_nOPV2=1.1;mean_mOPV2=1.5;mean_bOPV=4.9;case_control_ratio=4;zero_mean=0.2;n_per_strata=100

  # # Number of SIAs they were exposed to (identical within strata)
  # nOPV2_SIAs=rep(rpois(n=N_case,lambda=mean_nOPV2),each=n_per_strata)
  # mOPV2_SIAs=rep(rpois(n=N_case,lambda=mean_mOPV2),each=n_per_strata)
  # bOPV_SIAs=rep(rpois(n=N_case,lambda=mean_bOPV),each=n_per_strata)

  # Number of SIAs they were exposed to (different within strata)
  nOPV2_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_nOPV2)
  mOPV2_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_mOPV2)
  bOPV_SIAs=rpois(n=n_per_strata*(N_case)*(1+case_control_ratio),lambda=mean_bOPV)

  # Number of doses they received
  nOPV2_doses=rbinom(n=length(nOPV2_SIAs),size=nOPV2_SIAs,prob=coverage)
  mOPV2_doses=rbinom(n=length(mOPV2_SIAs),size=mOPV2_SIAs,prob=coverage)
  bOPV_doses=rbinom(n=length(bOPV_SIAs),size=bOPV_SIAs,prob=coverage)

  # Number of total OPV SIA doses
  OPV_doses=nOPV2_doses+mOPV2_doses+bOPV_doses

  # Strata for matched case-control (zero mean normally distributed)
  strata=rep(1:(N_case),each=n_per_strata)
  # Log odds from strata
  strata_value_log_OR=rep(rnorm(n=max(strata),mean=0,sd=1),each=n_per_strata)
  tab<-cbind.data.frame(nOPV2_SIAs,nOPV2_doses,
                        mOPV2_SIAs,mOPV2_doses,
                        bOPV_SIAs,bOPV_doses,
                        OPV_doses,strata,strata_value_log_OR)

  # Number of total doses reported
  # matplot(c(0,1:10),t(sapply(c(0.4,1:10),function(i){qlnorm(p=c(0.025,0.5,0.975),meanlog=log(i),sdlog=0.2)})),type="l")
  tab<-tab[order(tab$OPV_doses),]
  a<-as.numeric(names(table(tab$OPV_doses)))
  b<-as.numeric(table(tab$OPV_doses))
  if(recall_variation>0) tab$OPV_doses_reported=unlist(lapply(1:length(a),function(i){round(rlnorm(b[i],meanlog=log(ifelse(a[i]==0,zero_mean,a[i])),sdlog=recall_variation))}))
  if(recall_variation==0) tab$OPV_doses_reported=tab$OPV_doses

  # Log odds from effectiveness
  log_OR_nOPV2=log(1-eff_nOPV2)

  # Log odds * doses
  tab$nOPV2_value_log_OR=log_OR_nOPV2*tab$nOPV2_doses

  # Log odds from effectiveness
  log_OR_mOPV2=log(1-eff_mOPV2)

  # Log odds * doses
  tab$mOPV2_value_log_OR=log_OR_mOPV2*tab$mOPV2_doses


  # Full regression formula (logit probability)
  tab$logit_p=tab$nOPV2_value_log_OR+tab$mOPV2_value_log_OR+tab$strata_value_log_OR

  # Case or control based on probability
  tab$case=rbinom(n=length(tab$logit_p),size=1,prob=inv.logit(tab$logit_p))

  # Select strata with at least 1 case and at least case:control ratio controls
  tab<-mutate(group_by(as_tibble(tab),strata),n_case=sum(case),n_control=sum(case==0))
  sub<-tab%>%filter(n_case>=1,n_control>=case_control_ratio)

  # Select out N_case cases and case_control_ratio controls
  # Cases
  sub_case<-sub[sub$case==1,]
  # Order randomly
  sub_case<-sub_case[sample(x=1:nrow(sub_case),size=nrow(sub_case),replace=F),]
  # Case count per strata
  sub_case<-mutate(group_by(as_tibble(sub_case),strata),count=1:n())
  # Select first case per strata
  sel<-which(sub_case$count==1)
  # Select randomly
  if(length(sel)>N_case) sel<-sel[sample(x=1:length(sel),size=N_case,replace=F)]
  sub_case<-sub_case[sel,]
  # Controls
  sub_control<-sub[sub$case==0&sub$strata%in%sub_case$strata,]
  # Order randomly
  sub_control<-sub_control[sample(x=1:nrow(sub_control),size=nrow(sub_control),replace=F),]
  # Control count per strata
  sub_control<-mutate(group_by(as_tibble(sub_control),strata),count=1:n())
  # Select first case_control_ratio controls per strata
  sel<-which(sub_control$count<=case_control_ratio)
  sub_control<-sub_control[sel,]
  # Combine
  data<-as.data.frame(bind_rows(sub_case,sub_control))

  # Assign valency of doses
  data$OPV_SIAs<-rowSums(data[,c("bOPV_SIAs","nOPV2_SIAs","mOPV2_SIAs")])
  data$nOPV2_doses_frac_recall<-data$OPV_doses_reported*data$nOPV2_SIAs/data$OPV_SIAs
  data$mOPV2_doses_frac_recall<-data$OPV_doses_reported*data$mOPV2_SIAs/data$OPV_SIAs

  # Regression
  modC<-clogit(case~mOPV2_doses_frac_recall+nOPV2_doses_frac_recall+strata(strata),data=data)

  p<-sapply(list(modC),function(m){x<-summary(m)$coef;2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[1,3]^2 + x[2,3]^2 )))})

  fn<-function(mod){
    m<-summary(mod)$coef[,c(1,3)]
    r<-as.data.frame(1-exp(t(sapply(1:2,function(i){m[i,1]+c(0,-1,1)*1.96*m[i,2]}))))
    names(r)<-c("eff","upr","lwr")
    return(r)
  }

  r<-cbind.data.frame(type=rep("Frac + recall",times=2),vaccine=c("mOPV2","nOPV2"),
                      rbind(fn(modC)),pvalue=rep(p,each=2))
  ret<-cbind.data.frame(r,N_case=N_case,n_per_strata=n_per_strata,coverage=coverage,recall_variation=recall_variation,
                        zero_mean=zero_mean,
                        case_control_ratio=case_control_ratio,
                        eff_nOPV2=eff_nOPV2,eff_mOPV2=eff_mOPV2,mean_nOPV2=mean_nOPV2,mean_mOPV2=mean_mOPV2,mean_bOPV=mean_bOPV,
                        check_case=sum(data$case),check_control=sum(data$case==0),inacc=sum(data$OPV_doses_reported!=data$OPV_doses))
  return(ret)}

# Find minimum detectable difference
# Option 1: true eff = 0.15 and SD = 0
# Option 2: true eff = 0.5 and SD = 0.8 (as above)
opts<-cbind.data.frame(eff_nOPV2=c(0.15,0.5),recall_variation=c(0,0.8))

# You may want to reduce N_boot and increase mOPV2_step for testing, otherwise this will take a long time
def.N_boot=500
def.mOPV2_step=0.05

res_k<-NULL
for(k in 1:nrow(opts)){
eff_mOPV2=opts$eff_nOPV2[k]+def.mOPV2_step
p=0
j=0
res<-NULL
while(p<0.8&eff_mOPV2<=1){
  print(paste0("mOPV2: ",eff_mOPV2,", nOPV2: ",opts$eff_nOPV2[k]))
  res_list=NULL
  for(i in 1:def.N_boot){
    res_list[[i]]<-cbind.data.frame(boot=i,
                                    simulate_case_control3(N_case=350,coverage=0.75,case_control_ratio=4,
                                                          recall_variation=opts$recall_variation[k],
                                                          eff_nOPV2=opts$eff_nOPV2[k],n_per_strata=1000,zero_mean=0.2,
                                                          eff_mOPV2=eff_mOPV2,mean_nOPV2=1.5,
                                                          mean_mOPV2=1.1,mean_bOPV=4.9))
    print(paste0("Boot: ",i))
  }
  j=j+1
  res[[j]]<-do.call("rbind",res_list)
  p<-(summarise(group_by(res[[j]]%>%filter(vaccine=='nOPV2',check_case==350,check_control==(350*4)),N_case,n_per_strata,check_case,zero_mean,check_control,recall_variation,coverage,case_control_ratio,
                             eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,type,vaccine),
                    signif_diff=sum(pvalue<0.01),n_boot=n(),p=signif_diff/n_boot))$p
  print(paste0("Proportion of simulations with significant difference: ",p))
  eff_mOPV2=eff_mOPV2+def.mOPV2_step
}
res_k[[k]]<-do.call("rbind",res)
}
res<-do.call("rbind",res_k)
write.csv(res,"data/demo_diff.csv",row.names=F)
res<-read.csv("data/demo_diff.csv")
all(res$check_case==350)
all(res$check_control==(350*4))
res<-res%>%filter(vaccine=='nOPV2')
res$eff_nOPV2_char<-paste0("True nOPV2 effectiveness: ",round(100*res$eff_nOPV2),"%")
summarise(group_by(res,N_case,n_per_strata,check_case,zero_mean,check_control,recall_variation,coverage,case_control_ratio,
                   eff_nOPV2,eff_mOPV2,mean_nOPV2,mean_mOPV2,mean_bOPV,type,vaccine,eff_nOPV2_char),
          signif_diff1=sum(pvalue<0.01),signif_diff2=sum(pvalue<0.05),n_boot=n(),p1=signif_diff1/n_boot,p2=signif_diff2/n_boot)%>%pivot_longer(cols=c(p1,p2))%>%
  ggplot()+geom_point(aes(x=eff_mOPV2,y=value,color=name))+geom_line(aes(x=eff_mOPV2,y=value,color=name))+geom_hline(yintercept=0.8,linetype=2)+
  scale_y_continuous(labels=scales::percent)+scale_x_continuous(labels=scales::percent)+facet_wrap(~eff_nOPV2_char,scales="free_x")+
  labs(x="True mOPV2 effectiveness",y="Percentage of simulations with significant difference",color="P-value")+scale_color_discrete(labels=c("<0.01","<0.05"))
ggsave("figs/min_det_diff.png",width=6,height=4)
