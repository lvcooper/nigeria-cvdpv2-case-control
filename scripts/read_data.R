# Read and combine AFP surveilllance data from different years
# Quality control check data
path_to_data<-"data/"

v<-NULL
year<-17:22
for(i in 1:length(year)){
v[[i]]<-read.csv(paste0(path_to_data,"20",year[i]," AFP CV Data.csv"),colClasses = "character")
v[[i]]$year<-as.numeric(paste0("20",year[i]))
}

sel<-c("year","date_onset","date_case_verified","dob","sex","lga_onset","lgas","state_onset","State","siaopv","riopv","ipv_ri_sia","ipv_ri","ipv_sia","last_opv_date","epid_num","age_months","age_years","SourceInfo_RI")
w<-v
for(i in 1:length(year)){
  w[[i]]<-v[[i]][,which(names(v[[i]])%in%sel)]
}
w<-do.call("bind_rows",w)
w$lga_onset[w$year==2017]<-w$lgas[w$year==2017]
w$state_onset[w$year==2017]<-w$State[w$year==2017]
w<-w[,!names(w)%in%c("lgas","State")]

# Fix duplicate EPIDs
w<-summarise(group_by(w,date_case_verified,date_onset,epid_num,dob,age_months,age_years,sex,riopv,siaopv,ipv_ri_sia,last_opv_date,lga_onset,state_onset,ipv_ri,ipv_sia,SourceInfo_RI),year=max(year))
duplicates<-mutate(group_by(w,epid_num),n=n())%>%filter(n>1)
w<-mutate(group_by(w,epid_num),n=n())%>%filter(n==1)

y<-summarise(group_by(duplicates,date_case_verified,date_onset,epid_num,dob,age_months,age_years,sex,riopv,siaopv,ipv_ri_sia,last_opv_date,lga_onset,state_onset),
          ipv_ri=unique(ipv_ri[!is.na(ipv_ri)]),
          ipv_sia=unique(ipv_sia[!is.na(ipv_sia)]),
          SourceInfo_RI=unique(SourceInfo_RI[!is.na(SourceInfo_RI)]),
          year=max(year))
mutate(group_by(y,epid_num),n=n())%>%filter(n>1)
duplicates<-duplicates[!duplicates$epid_num%in%y$epid_num,]
duplicates

w<-rbind(w,y)

# Lab result from POLIS
z<-read.csv(paste0(path_to_data,"Nigeria POLIS AFP 2017 2022.csv"))
z$Year<-format(as.Date(z$Date.Onset,"%d/%m/%Y"),"%Y")
z$EPID<-as.character(z$EPID)

# Lab result from Nigeria AFP DB
q<-read.csv(paste0(path_to_data,"2017_2022_AFP information_230223.csv"))
q$EPID<-as.character(q$EPID)

table(q$Year[q$Finalclassification=="cVDPV2"])
table(q$Finalclassification=="cVDPV2")
table(w$epid_num%in%z$EPID)
table(w$epid_num%in%q$EPID)
table(w$epid_num%in%c(z$EPID,q$EPID))

z$class_polis<-NA
z$class_polis[grepl("cVDPV2",z$Virus.Type.s.)&z$Classification=="VDPV"]<-"Case"
z$class_polis[!grepl("VDPV",z$Virus.Type.s.)&z$Classification=="Discarded"]<-"Control"

q$class_nga<-NA
q$class_nga[q$Virus.Type.s.=="cVDPV2"&q$Finalclassification=="cVDPV2"]<-"Case"
q$class_nga[q$Virus.Type.s.%in%c("2-Negative","3-NPENT")&q$Finalclassification=="Discarded"]<-"Control"

w<-left_join(w,z[,c("EPID","class_polis")],by=c("epid_num"="EPID"))
w<-left_join(w,q[,c("EPID","class_nga")],by=c("epid_num"="EPID"))

# Cross-check number of AFP POLIS vs CV vs AFP NGA
table(w$year)
table(q$Year)
table(z$Year)

# Cross-check number of cVDPV2 cases
tapply(z$class_polis=="Case",z$Year,sum,na.rm=T)
tapply(q$class_nga=="Case",q$Year,sum,na.rm=T)
tapply(ifelse(is.na(w$class_nga),w$class_polis,w$class_nga)=="Case",w$year,sum,na.rm=T)

# Cases positive by contact
r<-read.csv(paste0(path_to_data,"20230530 IndexFlaggedPositiveByContact in comment.csv"))
w$pos_cont<-w$epid_num%in%r$EPID

# Needed variables
# date_onset (date of onset of paralysis)
# date_case_verified (investigation date)
# dob (date of birth)
# sex (sex)
# lgas (LGA - ADM2)
# State (state - ADM1)
# siaopv (OPV doses from SIAs)
# riopv (OPV doses from RI)
# ipv_ri_sia (total IPV doses from both RI and SIA)
# last_opv_date (last date of OPV)
# epid (EPID - need for linking to lab result)

# Completeness by year - need composite age
tapply(is.na(w$date_onset),w$year,mean)
tapply(is.na(w$date_case_verified),w$year,mean)
tapply(is.na(w$dob),w$year,mean)
tapply(is.na(w$age_months),w$year,mean)
tapply(is.na(w$age_years),w$year,mean)
tapply(is.na(w$sex),w$year,mean)
tapply(is.na(w$riopv),w$year,mean)
tapply(is.na(w$siaopv),w$year,mean)
tapply(is.na(w$ipv_ri_sia),w$year,mean)
tapply(is.na(w$last_opv_date),w$year,mean)
tapply(is.na(w$lga_onset),w$year,mean)
tapply(is.na(w$state_onset),w$year,mean)
tapply(is.na(w$ipv_ri),w$year,mean)
tapply(is.na(w$ipv_sia),w$year,mean)
tapply(is.na(w$SourceInfo_RI),w$year,mean)

v<-w
rm(w)

# Dates
v$date<-ifelse(v$year==2018,as.Date(v$date_onset,format="%d/%m/%Y"),as.Date(as.numeric(v$date_onset),origin=as.Date("1900-01-01")))
v$date<-as.Date(as.numeric(v$date),origin=as.Date("1970-01-01"))

v$investigation_date<-ifelse(v$year==2018,as.Date(v$date_case_verified,format="%d/%m/%Y"),as.Date(as.numeric(v$date_case_verified),origin=as.Date("1900-01-01")))
v$investigation_date<-as.Date(as.numeric(v$investigation_date),origin=as.Date("1970-01-01"))

v$dob_orig<-v$dob
v$dob<-ifelse(v$year==2018,as.Date(v$dob_orig,format="%d/%m/%Y"),as.Date(as.numeric(v$dob_orig),origin=as.Date("1900-01-01")))
v$dob<-as.Date(as.numeric(v$dob),origin=as.Date("1970-01-01"))
v$age_orig<-365*ifelse(is.na(as.numeric(v$age_years)),0,as.numeric(v$age_years))+365/12*ifelse(is.na(as.numeric(v$age_months)),0,as.numeric(v$age_months))
v$age_orig[is.na(v$age_months)&is.na(v$age_years)]<-NA
v$age_orig[!is.na(v$age_months)&!is.na(v$age_years)]

v$dob[is.na(v$dob)]<-(v$investigation_date-
                                  365*ifelse(is.na(as.numeric(v$age_years)),0,as.numeric(v$age_years))-
                                  365/12*ifelse(is.na(as.numeric(v$age_months)),0,as.numeric(v$age_months)))[is.na(v$dob)]

v$last_opv_date_orig<-v$last_opv_date
v$last_opv_date<-ifelse(v$year==2017,as.Date(v$last_opv_date,format="%d-%b-%y"),as.Date(as.numeric(v$last_opv_date),origin=as.Date("1900-01-01")))
v$last_opv_date<-as.Date(as.numeric(v$last_opv_date),origin=as.Date("1970-01-01"))

v$doses_ipv_ri<-as.numeric(v$ipv_ri)
v$doses_ipv_sia<-as.numeric(v$ipv_sia)
v$doses_ipv<-as.numeric(v$ipv_ri_sia)

v$doses_ipv_ri[which(v$doses_ipv_ri==99)]<-NA
v$doses_ipv_sia[which(v$doses_ipv_sia==99)]<-NA
v$doses_ipv[which(v$doses_ipv==99)]<-NA

v$doses_opv_sia<-as.numeric(v$siaopv)
v$doses_opv_ri<-as.numeric(v$riopv)
v$doses_opv_sia[which(v$doses_opv_sia==99)]<-NA
v$doses_opv_ri[which(v$doses_opv_ri==99)]<-NA

table(v$riopv,v$doses_opv_ri,useNA="ifany")
table(v$siaopv,v$doses_opv_sia,useNA="ifany")
table(v$ipv_ri,v$doses_ipv_ri,useNA="ifany")
table(v$ipv_sia,v$doses_ipv_sia,useNA="ifany")
table(v$ipv_ri_sia,v$doses_ipv,useNA="ifany")

clean<-function(x){gsub("/","",gsub("-","",gsub(" ","",gsub(",","",tolower(x)),fixed=T),fixed=T),fixed=T)}

v$adm1<-clean(v$state_onset)

v$adm2<-clean(v$lga_onset)

v$sex[v$sex=="FEMALE"]<-"F"
v$sex[v$sex=="MALE"]<-"M"

table(paste0("NGA: ",v$class_nga),paste0("POLIS: ",v$class_polis))
v$class<-ifelse(is.na(v$class_nga),v$class_polis,v$class_nga)

v$EPID<-v$epid_num

v$source<-ifelse(v$SourceInfo_RI%in%c("0","epi_registry","n/a"),NA,v$SourceInfo_RI)

cols<-c("EPID","date","investigation_date","year","dob","dob_orig","age_orig","sex","adm1","adm2","class","pos_cont","doses_opv_sia","doses_opv_ri","doses_ipv","doses_ipv_ri","doses_ipv_sia","last_opv_date","source")


d<-v[!is.na(v$class),cols]

# Check date quality 
# Onset date should come before investigation date
d$diff1<-as.numeric(as.Date(d$investigation_date)-as.Date(d$date))
# ggplot()+geom_histogram(data=d,aes(x=diff1),binwidth=1)
table(d$diff1<=0,useNA="ifany")
table(d$year[which(d$diff1<=0)])

# table(format(as.Date(d$investigation_date[which(d$diff1<=0)]),"%d"))
# 
# d$investigation_date_flipped<-d$investigation_date
# d$investigation_date_flipped[which(d$diff1<=0)]<-as.Date(paste0(format(as.Date(d$investigation_date),"%m"),"-",
#                                              format(as.Date(d$investigation_date),"%d"),"-",
#                                              format(as.Date(d$investigation_date),"%Y")),format="%d-%m-%Y")[which(d$diff1<=0)]
# d$investigation_date_fixed<-d$investigation_date
# d$investigation_date_fixed[which(as.Date(d$investigation_date)<=as.Date(d$date)&as.Date(d$investigation_date_flipped)>as.Date(d$date))]<-
#   d$investigation_date_flipped[which(as.Date(d$investigation_date)<=as.Date(d$date)&as.Date(d$investigation_date_flipped)>as.Date(d$date))]
# 
# d$diff1<-as.numeric(as.Date(d$investigation_date_fixed)-as.Date(d$date))
# ggplot()+geom_histogram(data=d,aes(x=diff1),binwidth=1)
# table(d$diff1<=0,useNA="ifany")
# table(d$year[which(d$diff1<=0)])

# Date last OPV should come before investigation date
d$diff2<-as.numeric(as.Date(d$investigation_date)-as.Date(d$last_opv_date))
table(d$diff2<=0,useNA="ifany")
table(d$year[which(d$diff2<=0)])

# Date of birth should come before onset date
d$diff3<-as.numeric(as.Date(d$date)-as.Date(d$dob))
table(d$diff3<=0,useNA="ifany")
table(d$year[which(d$diff3<=0)])

# Date of birth should come before investigation date
d$diff4<-as.numeric(as.Date(d$investigation_date)-as.Date(d$dob))
table(d$diff4<=0,useNA="ifany")
table(d$year[which(d$diff4<=0)])

d$dates_correct<-d$diff1>0&d$diff3>0&d$diff4>0
# d$dates_correct[which(!is.na(d$last_opv_date)&d$diff2<=0)]<-F

# # Check age patterns
# ggplot()+geom_histogram(data=d,aes(x=floor(as.numeric(investigation_date-dob)/365*12)),binwidth=1)+
#   geom_vline(data=cbind.data.frame(x=seq(12,60,12)),aes(xintercept=x),alpha=0.5)+
#   facet_wrap(~year,scales="free")+scale_x_continuous(limits=c(-1,61))
# 
# 
# ggplot()+geom_histogram(data=d,aes(x=wday(as.Date(dob))),binwidth=1)
# ggplot()+geom_histogram(data=d,aes(x=as.numeric(format(as.Date(dob),"%d"))),binwidth=1)
# ggplot()+geom_histogram(data=d,aes(x=as.numeric(format(as.Date(dob),"%m"))),binwidth=1)
# 
# ggplot()+geom_histogram(data=d,aes(x=wday(as.Date(date))),binwidth=1)
# ggplot()+geom_histogram(data=d,aes(x=as.numeric(format(as.Date(date),"%d"))),binwidth=1)
# ggplot()+geom_histogram(data=d,aes(x=as.numeric(format(as.Date(date),"%m"))),binwidth=1)
# 
# ggplot()+geom_histogram(data=d,aes(x=wday(as.Date(investigation_date))),binwidth=1)
# ggplot()+geom_histogram(data=d,aes(x=as.numeric(format(as.Date(investigation_date),"%d"))),binwidth=1)
# ggplot()+geom_histogram(data=d,aes(x=as.numeric(format(as.Date(investigation_date),"%m"))),binwidth=1)
# 
# 
# table(table(d$EPID))
# 
# table(d$doses_opv_ri,useNA="ifany")

table(d$doses_opv_sia,useNA="ifany")

# Match to GUID
load("data/nigeria_maps.Rdata")
key<-unique(st_drop_geometry(adm2data)[,c("adm2","adm1","GUID","CENTER_LAT","CENTER_LON")])


d$adm1[d$adm1=="fctabuja" & !is.na(d$adm1)]<-"fct"
d$adm2[d$adm2=="ayedade"&d$adm1=="osun" & !is.na(d$adm2) & !is.na(d$adm1)]<-"aiyedade" 
d$adm2[d$adm2=="ayedire"&d$adm1=="osun" & !is.na(d$adm2)& !is.na(d$adm1)]<-"aiyedire" 
d$adm2[d$adm2=="girei"&d$adm1=="adamawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"girie" 
d$adm2[d$adm2=="lamurde"&d$adm1=="adamawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"larmurde" 
d$adm2[d$adm2=="toungo"&d$adm1=="adamawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"teungo" 
d$adm2[d$adm2=="buruku"&d$adm1=="benue" & !is.na(d$adm2)& !is.na(d$adm1)]<-"bukuru" 
d$adm2[d$adm2=="burutu"&d$adm1=="benue" & !is.na(d$adm2)& !is.na(d$adm1)]<-"bukuru" 
d$adm2[d$adm2=="onuimo"&d$adm1=="imo" & !is.na(d$adm2)& !is.na(d$adm1)]<-"unuimo" 
d$adm2[d$adm2=="birniwa"&d$adm1=="jigawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"biriniwa" 
d$adm2[d$adm2=="suletankarkar"&d$adm1=="jigawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"suletankakar" 
d$adm2[d$adm2=="kabau"&d$adm1=="kaduna" & !is.na(d$adm2)& !is.na(d$adm1)]<-"kuban" 
d$adm2[d$adm2=="ungogo"&d$adm1=="kano" & !is.na(d$adm2)& !is.na(d$adm1)]<-"ungongo" 
d$adm2[d$adm2=="munya"&d$adm1=="niger" & !is.na(d$adm2)& !is.na(d$adm1)]<-"muya" 
d$adm2[d$adm2=="yewanorth"&d$adm1=="ogun" & !is.na(d$adm2)& !is.na(d$adm1)]<-"egbadonorth" 
d$adm2[d$adm2=="yewasouth"&d$adm1=="ogun" & !is.na(d$adm2)& !is.na(d$adm1)]<-"egbadosouth" 
d$adm2[d$adm2=="akokoakokosoutheast"&d$adm1=="ondo" & !is.na(d$adm2)& !is.na(d$adm1)]<-"akokosoutheast" 
d$adm2[d$adm2=="wamakko"&d$adm1=="sokoto" & !is.na(d$adm2)& !is.na(d$adm1)]<-"wamako" 
d$adm2[d$adm2=="bade"&d$adm1=="yobe" & !is.na(d$adm2)& !is.na(d$adm1)]<-"barde" 
d$adm2[d$adm2=="bursari"&d$adm1=="yobe" & !is.na(d$adm2)& !is.na(d$adm1)]<-"borsari" 
d$adm2[d$adm2=="tarmuwa"&d$adm1=="yobe" & !is.na(d$adm2)& !is.na(d$adm1)]<-"tarmua"
d$adm2[d$adm2=="kubau"&d$adm1=="kaduna" & !is.na(d$adm2)& !is.na(d$adm1)]<-"kuban" 
d$adm2[d$adm2=="kirikasamma"&d$adm1=="jigawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"kirikasama" 
d$adm2[d$adm2=="malammadori"&d$adm1=="jigawa" & !is.na(d$adm2)& !is.na(d$adm1)]<-"malammaduri" 
d$adm2[d$adm2=="jamaare"&d$adm1=="bauchi" & !is.na(d$adm2)& !is.na(d$adm1)]<-"jama'are" 
d$adm2[d$adm2=="jemaa"&d$adm1=="kaduna" & !is.na(d$adm2)& !is.na(d$adm1)]<-"jema'a" 
d$adm2[d$adm2=="quaanpan"&d$adm1=="plateau" & !is.na(d$adm2)& !is.na(d$adm1)]<-"qua'anpan"  
d$adm2[d$adm2=="maiadua"&d$adm1=="katsina" & !is.na(d$adm2)& !is.na(d$adm1)]<-"mai'adua" 


#check which ones are not in shapefile
unique(d$adm1[which(!d$adm1 %in% unique(key$adm1))])
unique(d$adm2[which(!d$adm2 %in% unique(key$adm2))]) # identifies "kubau"

unique(d[which(!d$adm2 %in% unique(key$adm2)),c("adm1","adm2")])


d<-left_join(d,key)

# Exclusion flow chart
total<-table(d$class)

d<-d[,!names(d)%in%c("diff1","diff2","diff3","diff4")]
write.csv(d,paste0(path_to_data,"all.csv"),row.names=F)

d$case<-d$class=="Case"

d$doses_ipv012<-d$doses_ipv
d$doses_ipv012[d$doses_ipv>2]<-2

d$age<-floor(as.numeric(as.Date(d$investigation_date)-as.Date(d$dob))/365*12)

table(paste0(ifelse(is.na(d$dob_orig)&is.na(d$age_orig),"Miss age",""),
             ifelse(is.na(d$GUID),"Miss GUID",""),
             ifelse(is.na(d$date),"Miss onset",""),
             ifelse(is.na(d$investigation_date),"Miss inv",""),
             ifelse(ifelse(d$dates_correct==T|is.na(d$dates_correct),T,F),"","Dates wrong")),d$class)

d<-d[!is.na(d$investigation_date),]
flow<-rbind(total,inv_date=table(d$class))

d<-d[!is.na(d$date),]
flow<-rbind(flow,ons_date=table(d$class))

d<-d[!is.na(d$dob),]
flow<-rbind(flow,known_age=table(d$class))

d<-d[!is.na(d$GUID),]
flow<-rbind(flow,guid=table(d$class))

d<-d[d$dates_correct,]
flow<-rbind(flow,dates_correct=table(d$class))

table(paste0(ifelse(is.na(d$doses_opv_sia),"Miss OPV",""),ifelse(is.na(d$doses_ipv),"Miss IPV","")),d$class)

d<-d[!is.na(d$doses_opv_sia),]
flow<-rbind(flow,doses_opv_sia=table(d$class))

d<-d[!is.na(d$doses_ipv),]
flow<-rbind(flow,doses_ipv=table(d$class))

table(paste0(ifelse(d$dob>as.Date("2016-05-01"),"","Before"),ifelse(d$age<60,"","Over 5")),d$class)

d<-d[d$dob>as.Date("2016-05-01"),]
flow<-rbind(flow,dob_switch=table(d$class))

d<-d[d$age<60,]
flow<-rbind(flow,under_5=table(d$class))

write.csv(d,paste0(path_to_data,"eligible.csv"),row.names=F)
write.csv(flow,paste0(path_to_data,"flow.csv"))
rm(list=ls())