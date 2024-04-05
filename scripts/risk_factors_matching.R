load("data/nigeria_maps.Rdata")

define.distance<-50000
define.age<-12
define.onset<-30
match.on.RI<-F
if(file.exists("data/all.csv")){a<-read.csv("data/all.csv")} else{return(warning("Need to first run read_data.R"))}
if(file.exists("data/eligible.csv")){d<-read.csv("data/eligible.csv")} else{return(warning("Need to first run read_data.R"))}
filename<-paste0("data/matched_distance_",define.distance,"_age_",define.age,"_onset_",define.onset,ifelse(match.on.RI,"_RI",""),".csv")
if(file.exists(filename)){m1<-read.csv(filename)} else{return(warning("Need to first run matching.R"))}

df<-a[which(a$class=="Case"),]
df$eligible<-df$EPID%in%d$EPID
df$matched<-df$EPID%in%m1$EPID
df<-merge(df,read.csv("data/zone_key.csv"))
df$age<-cut(as.numeric(as.Date(df$date)-as.Date(df$dob))/365,c(-1,2,3,4,5,100),right=F)
df$year<-as.factor(df$year)
df$doses_opv_sia<-cut(df$doses_opv_sia,c(-1,0,10,20,40))
df$doses_opv_ri<-as.factor(df$doses_opv_ri>2)
df$doses_ipv<-as.factor(df$doses_ipv>0)
df$north<-ifelse(grepl("N",df$Zone),df$Zone,"S")

mod<-glm(eligible~doses_opv_sia,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~doses_opv_ri,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~doses_ipv,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~age,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~sex,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~pos_cont,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~year,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~0+north,data=df,family="binomial")
summary(mod)


# Cases 3-5 years of age are less likely to be eligible
# Cases with onset in 2021 are more likely to be eligible
mod<-glm(eligible~age+doses_opv_sia+year+north,data=df,family="binomial")
summary(mod)
mod<-glm(eligible~age+year,data=df,family="binomial")
summary(mod)

# Eligible cases over 4 years or from 2018 and 2019 are less likely to be matched
mod<-glm(matched~doses_opv_sia+age+doses_opv_ri+doses_ipv+sex+pos_cont+year+north,data=df[df$eligible,],family="binomial")
summary(mod)
mod<-glm(matched~age+year,data=df[df$eligible,],family="binomial")
summary(mod)

mod<-glm(matched~doses_opv_sia+age+doses_opv_ri+doses_ipv+sex+pos_cont+year+north,data=df,family="binomial")
summary(mod)
mod<-glm(matched~age+year,data=df,family="binomial")
summary(mod)

p<-cbind.data.frame(med=boot::inv.logit(mod$coefficients),
                    lwr=boot::inv.logit(mod$coefficients+1.96*summary(mod)$coef[,2]),
                    upr=boot::inv.logit(mod$coefficients-1.96*summary(mod)$coef[,2]),
                    adm1=gsub("as.factor(year)","",names(mod$coefficients),fixed=T))
ggplot()+geom_pointrange(data=p,aes(x=adm1,y=med,ymin=lwr,ymax=upr))+labs(x="Year",y="Matching rate")



a$case<-grepl("Case",a$class)
a$reason_exclude<-ifelse(a$EPID%in%m1$EPID,"Matched","No match")
a$reason_exclude[is.na(a$doses_opv_sia)]<-"No OPV SIA doses"
a$reason_exclude[a$dob<=as.Date("2016-05-01")]<-"Born before switch"
a$reason_exclude<-factor(a$reason_exclude,levels=c("Born before switch","No OPV SIA doses","No match","Matched"),ordered=T)
df<-summarise(group_by(a,case,reason_exclude,adm1),n=n())%>%filter(case)
df$adm1<-factor(df$adm1,levels=names(sort(table(a$adm1[a$case]))),ordered=T)
ggplot()+geom_col(data=df,aes(x=adm1,y=n,fill=reason_exclude))+facet_wrap(~case,scales="free")+coord_flip()
df%>%filter(adm1%in%c("jigawa","kebbi","borno","katsina"))
df<-summarise(group_by(a%>%filter(case),adm1),matched=sum(reason_exclude=="Matched"),total=sum(grepl("match",reason_exclude,ignore.case=T)),
              p=matched/total,lwr=ifelse(total>0,binom.test(matched,total)$conf[1],NA),upr=ifelse(total>0,binom.test(matched,total)$conf[2],NA))
df$adm1<-factor(df$adm1,levels=df$adm1[order(df$upr)],ordered=T)
ggplot(data=df)+geom_pointrange(aes(x=adm1,y=p,ymin=lwr,ymax=upr))+coord_flip()
df%>%filter(adm1%in%c("taraba","katsina","borno","kano","sokoto"))

# Katsina cluster
sort(table(a$adm2[a$case]))
a%>%filter(case&adm1%in%names(tail(sort(table(a$adm1[a$case])),6)))%>%ggplot()+geom_histogram(aes(x=as.Date(date)),binwidth=30)+facet_wrap(~adm1)
# Jigawa, Dutse, 25 cases between May and October 2021 (12 by contact)
a%>%filter(adm2=="dutse"&case)
# Dukku, Gombe, 11 cases (3 by contact)
a%>%filter(adm2=="dukku"&case)
a%>%filter(adm2=="baure"&case)
a%>%filter(adm2=="nganzai"&case)

a%>%filter(case)%>%ggplot()+geom_histogram(aes(x=as.Date(date)),binwidth=30)
