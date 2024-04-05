rm(list=ls())
load("data/nigeria_maps.Rdata")
# Read in raw survey data
d<-read.csv("data/control_survey_data/combined.csv")
# Get EPID (remove formatting)
d$case.epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",d$Case.EPID),fixed=T),fixed=T),fixed=T)
# Is the EPID a contact?
d$contact<-grepl("C",substring(d$case.epid,first=nchar(d$case.epid)-3,last=nchar(d$case.epid)),ignore.case=T)
# EPID root
# Remove contact suffix if required
get_epid_root<-function(epid){
  first_half<-substring(epid,first=1,last=nchar(epid)-5)
  second_half<-substring(epid,first=nchar(epid)-4,last=nchar(epid))
  second_half<-strsplit(second_half,"C")[[1]][1]
  return(paste0(first_half,second_half))
}
d$epid_root<-sapply(d$case.epid,get_epid_root)

# We have 
nrow(d) 
# controls, from 
length(unique(d$File.Name)) 
# surveys with 
length(unique(d$case.epid))
# unique EPIDs. Of these, 
sum(unique(d[,c("case.epid","contact")])$contact) 
# have EPIDs that indicate that they are contacts. We have 
length(unique(d$epid_root)) 
# unique EPID roots (i.e. removing the contact suffix).

# Clean healthy control case EPIDs, cross-checking with national database
# National CV AFP database
a<-read.csv("data/all.csv")
a$case<-a$class=="Case"
a$epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",a$EPID),fixed=T),fixed=T),fixed=T)
national<-a$epid[a$case]

d$case.epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",d$Case.EPID),fixed=T),fixed=T),fixed=T)
d$contact<-grepl("C",substring(d$case.epid,first=nchar(d$case.epid)-2,last=nchar(d$case.epid)),ignore.case=T)
d$epid_root<-sapply(d$case.epid,get_epid_root)
# 47 case EPIDs which are not in the national database (41 not at all, 6 not classed as cases)
length(unique(d$epid_root[!d$epid_root%in%national]))
length(unique(d$epid_root[!d$epid_root%in%a$epid]))
# Total of 271
length(unique(d$epid_root))

# Too short?
unique(d$epid_root[nchar(d$epid_root)<14])
d$case.epid<-gsub("NIEKTSKAT004","NIEKTSKAT21004",d$case.epid)
d$case.epid<-gsub("NEBASZAK21012","NIEBASZAK21012",d$case.epid)
d$case.epid<-gsub("NIEJISDUT060","NIEJISDUT21060",d$case.epid)
d$case.epid<-gsub("NIEBOSMAG2120","NIEBOSMAG21020",d$case.epid)
d$case.epid<-gsub("NIEKBSGWN015","NIEKBSGWN21015",d$case.epid)
d$case.epid<-gsub("NIEKBSGWN021","NIEKBSGWN21021",d$case.epid)
d$epid_root<-sapply(d$case.epid,get_epid_root)
# Fixed 5 case EPIDs which are not in the national database
length(unique(d$epid_root[!d$epid_root%in%national]))

# Too long?
unique(d$epid_root[nchar(d$epid_root)>14])
d$case.epid<-gsub("BOSMAFA","BOSMAF",d$case.epid)
d$case.epid<-gsub("AESI","",d$case.epid)
d$case.epid<-gsub("NIEGMSFKY2021021","NIEGMSFKY21021",d$case.epid)
d$case.epid<-gsub("NIEGMSKWM2021024","NIEGMSKWM21024",d$case.epid)
d$case.epid<-gsub("NIEGMSFKY2021023","NIEGMSFKY21023",d$case.epid)
d$case.epid<-gsub("NIEJISDUT210154","NIEJISDUT21154",d$case.epid)
d$case.epid<-gsub("NIEOYSKSH021008","NIEOYSKSH21008",d$case.epid)
d$case.epid<-gsub("NIEPLSWASE21022","NIEPLSWAS21022",d$case.epid)
d$epid_root<-sapply(d$case.epid,get_epid_root)
# Fixed 12 case EPIDs which are not in the national database
length(unique(d$epid_root[!d$epid_root%in%national]))

# Something else
unique(d$epid_root[!d$epid_root%in%national])
d$case.epid<-gsub("NIEBAUDBM","NIEBASDBM",d$case.epid)
d$case.epid<-gsub("NIEJISJHN21065","NIEJISDUT21065",d$case.epid)
d$case.epid<-gsub("NIEJISJHN21071","NIEJISDUT21071",d$case.epid)
d$case.epid<-gsub("NISJISDUT21130","NIEJISDUT21130",d$case.epid)
d$case.epid<-gsub("NIEBOSGUB02114","NIEBOSGUB21014",d$case.epid)
d$case.epid<-gsub("NIESOSGRY21OO6","NIESOSGRY21006",d$case.epid)
d$case.epid<-gsub("NIEKBSMYM21020","NIEKBSSNA21020",d$case.epid)
d$case.epid<-gsub("NIESOSRBH21003","NIESOSRBA21003",d$case.epid)
d$case.epid<-gsub("NIEZASMRR22017","NIEZASMRR22010",d$case.epid)
d$epid_root<-sapply(d$case.epid,get_epid_root)
# Fixed 9 case EPIDs which are not in the national database
length(unique(d$epid_root[!d$epid_root%in%national]))

# Check against POLIS
b<-read.csv("data/Nigeria POLIS AFP 2017 2022.csv",sep=";")
b<-b[b$Place.Admin.0=="NIGERIA"&b$Surveillance.Type=="AFP",]
b$case<-grepl("cVDPV2",b$Virus.Type.s.)
b$epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",b$EPID),fixed=T),fixed=T),fixed=T)
polis<-b$epid
# Only 2023 remaining
length(unique(d$epid_root[!d$epid_root%in%polis]))

# Merge to POLIS case date
b$case_onset_polis<-as.Date(b$Date.Onset,"%d/%m/%Y")
d<-merge(d,b[,c("epid","case_onset_polis")],all.x=T,by.x="epid_root",by.y="epid")
# Merge to national AFP CV DB case date
a$case_onset_nat<-as.Date(a$date)
d<-merge(d,a[,c("epid","case_onset_nat")],all.x=T,by.x="epid_root",by.y="epid")

d$case_onset<-as.Date(ifelse(is.na(d$case_onset_nat),d$case_onset_polis,d$case_onset_nat),origin=as.Date("1970-01-01"))

# Clean formatting of dates
fix_dates<-function(char_vector,min_date,max_date){
  # char_vector<-"2020-03-01";min_date=as.Date("2022-05-12");max_date=Sys.Date()
  char<-gsub(".","",gsub("th","",gsub("\\","-",gsub("/","-",gsub(" ","-",gsub(",","-",str_squish(char_vector)))),fixed=T)),fixed=T)
  formats<-c("%d-%m-%Y","%d-%m-%y","%m-%d-%Y","%Y-%m-%d","%Y-%d-%m")
  dates<-sapply(1:length(formats),function(f){as.Date(char,format=formats[f])})
  dates<-dates[!is.na(dates)]
  dates<-as.Date(dates,origin=as.Date("1970-01-01"))
  if(length(dates)==1) return(dates)
  dates_sub<-dates[dates<max_date&dates>min_date]
  if(length(dates_sub)==0) return(NA) 
  if(length(dates_sub)>0) return(sample(dates_sub,1))
}

# Survey dates between April 2021 and present, should be after case onset
k<-unique(d[,c("Survey.Date","case_onset","case_onset_polis")])
k$survey_date<-as.Date(unlist(sapply(1:nrow(k),function(i){fix_dates(k$Survey.Date[i],min_date=as.Date(ifelse(is.na(k$case_onset[i]),'2021-04-01',k$case_onset[i]),origin=as.Date("1970-01-01")),max_date=Sys.Date())})),origin=as.Date("1970-01-01"))
k[is.na(k$survey_date),]
k$survey_date[k$Survey.Date=="22/4/2022"]<-as.Date("2022-04-22")
k$survey_date[k$Survey.Date=="14th August 2021"]<-as.Date("2021-08-14")
k$survey_date[k$Survey.Date=="10TH NOVEMBER 2021"]<-as.Date("2021-11-10")
k$survey_date[k$Survey.Date=="31/8/21"]<-as.Date("2021-08-31")
k$survey_date[k$Survey.Date=="4th August 2021"]<-as.Date("2021-08-04")
k$survey_date[k$Survey.Date=="23/10/23"]<-as.Date("2023-10-23")
k$survey_date[k$Survey.Date=="14th October 2021"]<-as.Date("2021-10-14")
k$survey_date[k$Survey.Date=="2021-07-30"]<-as.Date("2021-07-30")
k$survey_date[k$Survey.Date=="17-18/11/2021"]<-as.Date("2021-11-18")
k$survey_date[k$Survey.Date=="23/9/2019"]<-as.Date("2019-09-23")
d<-merge(d,k,by=c("Survey.Date","case_onset","case_onset_polis"))


# Clean age
# Trim whitespace case epid
d$Case.EPID<-trimws(d$Case.EPID)
v<-gsub("!","",gsub("m","",gsub("mon","",gsub("month","",gsub("months","",gsub("mnt","",gsub("mnth","",gsub("mnths","",gsub(" ","",tolower(d$Age.month))))))))))
v[is.na(as.numeric(v))]
d$Age.month[which(is.na(as.numeric(v)))]
d$age_month<-as.numeric(v)
d$age_month[grep("wk",v)]<-floor(as.numeric(gsub("wk","",v[grep("wk",v)]))/4)
d$age_month[grep("days",v)]<-floor(as.numeric(gsub("days","",v[grep("days",v)]))/30)
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==1]<-24
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==2]<-24
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==3]<-27
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==4]<-36
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==5]<-48
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==6]<-36
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==7]<-27
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==8]<-48
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==9]<-39
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==10]<-24
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==11]<-40
d$age_month[d$Case.EPID=="NIE-SOS-BDN-21-002"&d$ID==12]<-48
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==1]<-21
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==2]<-36
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==3]<-32
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==4]<-18
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==5]<-27
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==6]<-19
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==7]<-22
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==8]<-29
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==10]<-37
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==11]<-15
d$age_month[d$Case.EPID=="NIE-KNS-BBJ-21-005C1"&d$ID==12]<-16
unique(d[is.na(d$age_month),c("Age.month","Case.EPID","ID","Date.of.Birth","survey_date")])

d$dob<-d$survey_date-d$age_month*365/12

# Date of last SIA
d$date_last_opv_sia<-as.Date(sapply(1:nrow(d),function(i){fix_dates(d$Date.of.last.OPV.dose.received.through.SIA[i],
                                                                    min_date=as.Date(ifelse(is.na(d$dob[i]),'2000-01-01',d$dob[i]-30),origin=as.Date("1970-01-01")),
                                                                    max_date=as.Date(ifelse(is.na(d$survey_date[i]),Sys.Date(),d$survey_date[i]),origin=as.Date("1970-01-01")))}),origin=as.Date("1970-01-01"))

d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="D1(31st Oct 2020)"]<-as.Date("2020-10-31")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="D1(18th SEPT 2021)"]<-as.Date("2021-09-18")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="D1(7th Nov 2019)"]<-as.Date("2019-11-07")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="April 14th 2021"]<-as.Date("2021-04-14")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="10-13 Mar 2021"]<-as.Date("2021-03-10")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="25-05.2021"]<-as.Date("2021-05-25")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="31/06/2021"]<-as.Date("2021-06-30")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="31/6/20"]<-as.Date("2020-06-30")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="21-20-2021"&d$Case.EPID=="NIE-FCT-BWR-21-013"]<-as.Date("2021-10-21")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="19/20/2021"&d$Case.EPID=="NIE/KTS/KAT/21/002C3"]<-as.Date("2021-09-19")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="21/92021"&d$Case.EPID=="NIE-JIS-BUJ-21-025C1"]<-as.Date("2021-09-21")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="81/10/2020"&d$Case.EPID=="NIE-TRS-TTM-21-005"]<-as.Date("2020-10-18")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="6-9/12/2021"]<-as.Date("2021-12-06")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="17/10/2021"]<-as.Date("2021-10-17")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="19/9/2021"]<-as.Date("2021-09-19")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="11 April 2021"]<-as.Date("2021-04-11")
d$date_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="13 April 2021"]<-as.Date("2021-04-13")

d$date_last_opv_sia[is.na(d$date_last_opv_sia)]<-
  as.Date(sapply(which(is.na(d$date_last_opv_sia)),function(i){fix_dates(d$Date.of.last.OPV.dose.received.through.SIA[i],
                                                                         min_date=as.Date(ifelse(is.na(d$dob[i]),'2000-01-01',d$dob[i]-30),origin=as.Date("1970-01-01")),
                                                                         max_date=Sys.Date())}),origin=as.Date("1970-01-01"))

d$date_last_opv_sia[is.na(d$date_last_opv_sia)]<-
  as.Date(sapply(which(is.na(d$date_last_opv_sia)),function(i){fix_dates(d$Date.of.last.OPV.dose.received.through.SIA[i],
                                                 min_date=as.Date('2000-01-01'),
                                                max_date=as.Date(ifelse(is.na(d$survey_date[i]),Sys.Date(),d$survey_date[i]),origin=as.Date("1970-01-01")))}),origin=as.Date("1970-01-01"))

d$date_last_opv_sia[is.na(d$date_last_opv_sia)]<-
  as.Date(sapply(which(is.na(d$date_last_opv_sia)),function(i){fix_dates(d$Date.of.last.OPV.dose.received.through.SIA[i],
                                                                         min_date=as.Date('2000-01-01'),
                                                                         max_date=Sys.Date())}),origin=as.Date("1970-01-01"))

unique(d[is.na(d$date_last_opv_sia),c("Date.of.last.OPV.dose.received.through.SIA")])

# Month of last OPV SIA
d$month_last_opv_sia<-format(d$date_last_opv_sia,"%m")
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("jan",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"01"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("feb",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"02"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("mar",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"03"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("apr",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"04"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("may",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"05"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("jun",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"06"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("jul",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"07"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("aug",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"08"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("sep",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"09"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("oct",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"10"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("nov",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"11"
d$month_last_opv_sia[is.na(d$month_last_opv_sia)&grepl("dec",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"12"

# Year of last OPV SIA
d$year_last_opv_sia<-format(d$date_last_opv_sia,"%Y")
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2003",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2003"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2004",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2004"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2005",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2005"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2006",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2006"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2007",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2007"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2008",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2008"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2009",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2009"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2010",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2010"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2011",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2011"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2012",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2012"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2013",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2013"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2014",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2014"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2015",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2015"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2016",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2016"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2017",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2017"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2018",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2018"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2019",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2019"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2020",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2020"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2021",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2021"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("2022",d$Date.of.last.OPV.dose.received.through.SIA,ignore.case=T)]<-"2022"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("-20",d$Date.of.last.OPV.dose.received.through.SIA,fixed=T)]<-"2020"
d$year_last_opv_sia[is.na(d$year_last_opv_sia)&grepl("-21",d$Date.of.last.OPV.dose.received.through.SIA,fixed=T)]<-"2021"
d$year_last_opv_sia[d$Date.of.last.OPV.dose.received.through.SIA=="NOV, 18"]<-"2018"
unique(d[is.na(d$date_last_opv_sia),c("Date.of.last.OPV.dose.received.through.SIA","month_last_opv_sia","year_last_opv_sia")])
d$date_last_opv_sia[is.na(d$date_last_opv_sia)]<-as.Date(paste0(d$year_last_opv_sia[is.na(d$date_last_opv_sia)],"-",d$month_last_opv_sia[is.na(d$date_last_opv_sia)],"-15"),format="%Y-%m-%d")

# Clean LGA/state
d$adm1<-trimws(toupper(d$State))
d$adm1[d$adm1=="GMOBA"]<-"GOMBE"
d$adm1[d$adm1=="JIAGWA"]<-"JIGAWA"
d$adm1[!d$adm1%in%adm1data$ADM1_NAME]

d$case_adm1<-trimws(toupper(d$Case.ADM1))
d$case_adm1[!d$case_adm1%in%adm1data$ADM1_NAME]

d$adm2<-trimws(toupper(d$LGA))
d$adm2[d$adm2=="DAMBAM"]<-"DAMBAN"
d$adm2[d$adm2=="MMC"]<-"MAIDUGURI"
d$adm2[d$adm2=="AMAC"]<-"MUNICIPAL AREA COUNCIL"
d$adm2[d$adm2=="KUBAU"]<-"KUBAN"
d$adm2[d$adm2=="UNGOGO"]<-"UNGONGO"
d$adm2[d$adm2=="BIRNIN MAGAJI"]<-"BIRNIN MAGAJI/KIYAW"
d$adm2[d$adm2=="YENAGOA"]<-"YENEGOA"
d$adm2[d$adm2=="YAMALTU DEBA"]<-"YAMALTU/DEBA"
d$adm2[d$adm2=="NSR"]<-"NASSARAWA"
d$adm2[d$adm2=="DANGE SHUNI"]<-"DANGE-SHUNI"
d$adm2[d$adm2=="WAMAKKO"]<-"WAMAKO"
d$adm2[d$adm2=="ARDOKOLA"]<-"ARDO-KOLA"
d$adm2[d$adm2=="BURSARI"]<-"BORSARI"
d$adm2[d$adm2%in%c("KOKO BESSE","KOKO-BESSE","KOKO BESSE")]<-"KOKO/BESSE"
d$adm2[d$adm2=="IDO/OSI"]<-"IDO-OSI"
d$adm2[d$adm2=="TALATAMAFARA"]<-"TALATA MAFARA"
d$adm2[d$adm2=="LAGOS MAILAND"]<-"LAGOS MAINLAND"
unique(d[!d$adm2%in%adm2data$ADM2_NAME,c("adm1","adm2")])

d$case_adm2<-trimws(toupper(d$Case.ADM2))
d$case_adm2[d$case_adm2=="DAMBAM"]<-"DAMBAN"
d$case_adm2[d$case_adm2=="MMC"]<-"MAIDUGURI"
d$case_adm2[d$case_adm2=="YENAGOA"]<-"YENEGOA"
d$case_adm2[d$case_adm2=="DANGE SHUNI"]<-"DANGE-SHUNI"
d$case_adm2[d$case_adm2=="WAMAKKO"]<-"WAMAKO"
d$case_adm2[d$case_adm2=="ARDOKOLA"]<-"ARDO-KOLA"
d$case_adm2[d$case_adm2=="BURSARI"]<-"BORSARI"
d$case_adm2[d$case_adm2%in%c("KOKO BESSE","KOKO-BESSE","KOKO BESSE","KOKO - BESSE")]<-"KOKO/BESSE"
d$case_adm2[d$case_adm2=="IDO/OSI"]<-"IDO-OSI"
d$case_adm2[d$case_adm2=="MUNICIPAL"]<-"MUNICIPAL AREA COUNCIL"
d$case_adm2[d$case_adm2=="KUBAU"]<-"KUBAN"
d$case_adm2[d$case_adm2=="UNGOGO"]<-"UNGONGO"
d$case_adm2[d$case_adm2=="NNEWI NORTH LGA"]<-"NNEWI NORTH"
d$case_adm2[d$case_adm2=="BIRNIN MAGAJI"]<-"BIRNIN MAGAJI/KIYAW"
unique(d[!d$case_adm2%in%adm2data$ADM2_NAME,c("case_adm1","case_adm2")])

d[d$case_adm1!=d$adm1&d$case_adm1!="",c("adm1","case_adm1","Case.EPID")]

d$adm1[d$case.epid=="NIEKNSSML21009"&d$adm1==""]<-"KANO"
d$adm2[d$case.epid=="NIEKNSSML21009"&d$adm2==""]<-"SUMAILA"
d$adm1[d$case.epid=="NIELASKJA21005"&d$adm1%in%c("BORNO","IMO","KWARA")]<-"LAGOS"
d$adm1[d$case.epid=="NIESOSRBA21003C1"&is.na(d$adm1)]<-"LAGOS"
d<-merge(d,st_drop_geometry(adm2data)[,c("ADM1_NAME","ADM2_NAME","GUID","CENTER_LAT","CENTER_LON")],by.x=c("adm1","adm2"),by.y=c("ADM1_NAME","ADM2_NAME"),all.x=T)
d[is.na(d$CENTER_LAT),c("case_adm1","case_adm2","adm1","adm2","case.epid")]

# Clean source data
d$source<-trimws(tolower(d$Source.of.RI.vaccination.information))
d$source[d$source%in%c("c","ri card","imm. card","carad","cart","1")]<-"card"
d$source[d$source%in%c("r","recal","history","re-call","recall-new born referal","no card","0")]<-"recall"
d$source[d$source%in%c("nil","nill","none","n","","-","card/recall","call","recard","99","recall/card")]<-NA
table(d$source,useNA="ifany")

# Clean dose data
unique(d$Total.OPV.doses.received.through.RI[is.na(as.numeric(d$Total.OPV.doses.received.through.RI))])
unique(d[is.na(as.numeric(d$Total.OPV.doses.received.through.SIA)),c("Total.OPV.doses.received.through.SIA")])
unique(d[is.na(as.numeric(d$Total.IPV.doses.received.through.RI)),c("Total.IPV.doses.received.through.RI")])
# Where a date is reported for IPV RI, assume this is the date of IPV1
d[d$Total.IPV.doses.received.through.RI%in%c("14/07/2020","2020-10-06","29/06/2021"),c("Total.IPV.doses.received.through.RI","dob")]
d$Total.IPV.doses.received.through.RI[d$Total.IPV.doses.received.through.RI%in%c("14/07/2020","2020-10-06","29/06/2021")]<-1
unique(d[is.na(as.numeric(d$Total.IPV.doses.received.through.SIA)),c("Total.IPV.doses.received.through.SIA")])

unique(d[is.na(as.numeric(d$Total.OPV.doses.received.through.RI)),c("Total.OPV.doses.received.through.RI","File.Name")])
unique(d[is.na(as.numeric(d$Total.OPV.doses.received.through.SIA)),c("Total.OPV.doses.received.through.SIA","File.Name")])
unique(d[is.na(as.numeric(d$Total.IPV.doses.received.through.RI)),c("Total.IPV.doses.received.through.RI","File.Name")])
unique(d[is.na(as.numeric(d$Total.IPV.doses.received.through.SIA)),c("Total.IPV.doses.received.through.SIA","File.Name")])
unique(d[is.na(as.numeric(d$Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case)),c("Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case","File.Name")])

x<-unique(d[is.na(as.numeric(d$Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case)),c("Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case")])


clean_dose<-function(x){
  x<-tolower(gsub(" ","",gsub("s","",x)))
  x[x%in%c("notapplicable",  "notappplicable", "n", "na","n/a","nil",",",99,999)]<-NA
  return(as.numeric(x))}


d$doses_opv_sia_overall<-clean_dose(d$Total.OPV.doses.received.through.SIA)
d$doses_opv_ri<-clean_dose(d$Total.OPV.doses.received.through.RI)
d$doses_ipv_sia<-clean_dose(d$Total.IPV.doses.received.through.SIA)
d$doses_ipv_ri<-clean_dose(d$Total.IPV.doses.received.through.RI)
d$doses_opv_sia_since<-clean_dose(d$Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case)

# Clean sex data
d$sex<-NA
d$sex[d$Sex%in%c("M","m"," M","M ","Male","MALE","M "," M")]<-"M"
d$sex[d$Sex%in%c("F","f"," F","F ","Female","FEMALE","FALSE"," F" , "F ","Female ")]<-"F"
d$sex[grepl("M",d$Sex)&is.na(d$sex)]<-"M"
d$sex[grepl("F",d$Sex)&is.na(d$sex)]<-"F"

 

select<-c("case.epid","survey_date","date_last_opv_sia","age_month","sex","dob","adm1","adm2","doses_opv_sia_overall","source","doses_opv_ri","doses_ipv_sia","doses_ipv_ri","doses_opv_sia_since","GUID","CENTER_LAT","CENTER_LON","File.Name")
write.csv(d[,select],"data/healthy_control_data.csv",row.names=F)
