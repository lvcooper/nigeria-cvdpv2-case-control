regression_boot<-function(define.distance,define.age,define.onset,match.on.RI=F,define.offset,use.date.of.last.opv=T,boot,IPV,cutoff=NA){
  
  filename<-paste0("data/matched_distance_",define.distance,
                   "_age_",define.age,
                   "_onset_",define.onset,
                   ifelse(match.on.RI,"_RI_",""),
                   "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),"_boot_",boot,".csv")
  m1<-read.csv(filename)
  if(is.na(cutoff)) cutoff<-Inf
  library(survival)
  
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
  
  reg<-function(b){
    if(IPV==T) mod<-clogit(case~mOPV2_ons+nOPV2_ons+as.factor(doses_ipv012)+strata(stratum),data=m1[m1$boot==b&m1$mOPV2_ons<cutoff&m1$nOPV2_ons<cutoff,])
    if(IPV==F) mod<-clogit(case~mOPV2_ons+nOPV2_ons+strata(stratum),data=m1[m1$boot==b&m1$mOPV2_ons<cutoff&m1$nOPV2_ons<cutoff,])
    return(cbind(tab(mod),boot=b))
  }
  
  r<-NULL
  for(i in 1:boot){r[[i]]<-reg(i)
  print(i)}
  
x<-do.call("rbind",r)

return(x)

}


regression_boot_nonlinear<-function(define.distance,define.age,define.onset,match.on.RI=F,define.offset,use.date.of.last.opv=T,boot,IPV,grps){
  # define.distance=50000;define.age=12;define.onset=30;match.on.RI=F;define.offset=28;use.date.of.last.opv=T;boot=1000;IPV=T;grps=6
  filename<-paste0("data/matched_distance_",define.distance,
                   "_age_",define.age,
                   "_onset_",define.onset,
                   ifelse(match.on.RI,"_RI_",""),
                   "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),"_boot_",boot,".csv")
  m1<-read.csv(filename)
  
  library(survival)
  
  tab<-function(mod){
    x<-cbind(summary(mod)$coef[,2],
             summary(mod)$conf[,3:4])
    x<-cbind(x,
             100*(1-x))
    x<-as.data.frame(x)
    names(x)<-c("OR","OR_lwr","OR_upr","eff","eff_upr","eff_lwr")
    x$doses<-c(rep(1:grps,2),1:2)
    x$vac<-rep(c("mOPV2","nOPV2","IPV"),each=grps)[1:nrow(x)]
    return(x)
  }
  
  s<-c(0:grps,100)
  m1$mOPV2_cat<-cut(m1$mOPV2_ons,s,right=F)
  m1$nOPV2_cat<-cut(m1$nOPV2_ons,s,right=F)
  
  reg<-function(b){
    if(IPV==T) mod<-clogit(case~mOPV2_cat+nOPV2_cat+as.factor(doses_ipv012)+strata(stratum),data=m1[m1$boot==b,])
    if(IPV==F) mod<-clogit(case~mOPV2_cat+nOPV2_cat+strata(stratum),data=m1[m1$boot==b,])
    return(cbind(tab(mod),boot=b))
  }
  
  r<-NULL
  for(i in 1:boot){r[[i]]<-reg(i)
  print(i)}
  
  x<-do.call("rbind",r)
  
  return(x)
  
}