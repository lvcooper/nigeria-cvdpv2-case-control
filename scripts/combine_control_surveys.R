# Read and combine Nigeria healthy community control survey data

subdir<-"data/control_survey_data/"
f<-list.files(subdir)
f<-f[f!="combined.csv"]

d_list<-NULL
for(i in 1:length(f)){
  if(grepl("csv",f[i])){d<-read.csv(paste0(subdir,f[i]),skip=10,stringsAsFactors=F,fileEncoding="latin1")}
  if(grepl("xls",f[i])&!f[i]=="NIE-KNS-GSW-21-005C1.xlsx"){d<-read.xls(paste0(subdir,f[i]),skip=10,stringsAsFactors=F,fileEncoding="latin1")}
  if(f[i]=="NIE-KNS-GSW-21-005C1.xlsx"){d<-read.xls(paste0(subdir,f[i]),skip=9,stringsAsFactors=F,fileEncoding="latin1")}
  d<-d[,names(d)[!grepl("name",names(d),ignore.case=T)]]
  d<-d[,-(1:2)]
  names(d)<-c("Date.of.Birth","Age.month","Sex","State","LGA","Ward","Total.OPV.doses.received.through.SIA",
              "Date.of.last.OPV.dose.received.through.SIA","Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case",
              "Total.IPV.doses.received.through.SIA","Date.of.last.IPV.dose.received.through.SIA","Total.OPV.doses.received.through.RI","Total.IPV.doses.received.through.RI","Source.of.RI.vaccination.information")
  d$Date.of.last.IPV.dose.received.through.SIA<-as.character(d$Date.of.last.IPV.dose.received.through.SIA)
  d$Age.month<-as.character(d$Age.month)
  d$Total.IPV.doses.received.through.SIA<-as.character(d$Total.IPV.doses.received.through.SIA)
  d$Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case<-as.character(d$Number.of.OPV.doses.received.through.SIA.since.date.of.onset.of.case)
  d$ID<-1:nrow(d)
  
  if(grepl("csv",f[i])){survey_date<-read.csv(paste0(subdir,f[i]),skip=1)[1,15]}
  if(grepl("xls",f[i])){survey_date<-read.xls(paste0(subdir,f[i]),skip=1)[1,15]}
  
  if(grepl("csv",f[i])){case<-read.csv(paste0(subdir,f[i]),skip=4)[1:2,]}
  if(grepl("xls",f[i])){case<-read.xls(paste0(subdir,f[i]),skip=4)[1:2,]}
  
  case_epid<-case[1,3]
  case_dob<-case[1,16]
  case_age_months<-case[1,18]
  case_adm1<-case[2,3]
  case_adm2<-case[2,7]
  case_ward<-case[2,13]
  case_onset<-case[2,18]
  
  if(f[i]=="NIE-KBS-WRR-21-013.xlsx"){
    case_epid<-case[1,3]
    case_dob<-case[1,12]
    case_age_months<-case[1,14]
    case_adm1<-case[2,3]
    case_adm2<-case[2,6]
    case_ward<-case[2,9]
    case_onset<-case[2,14]
  }
  d_list[[i]]<-cbind.data.frame(d,File.Name=f[i],Survey.Date=survey_date,Case.EPID=case_epid,Case.ADM1=case_adm1,Case.ADM2=case_adm2,Case.Ward=case_ward,Case.Onset=case_onset,Case.DOB=case_dob,Case.Age.Months=case_age_months)
  print(i)
}

# Do all sheets have the same number of columns?
length(unique(sapply(d_list,ncol)))==1
# If not, which one is the odd one out?
x<-which(sapply(d_list,ncol)!=names(tail(sort(table(sapply(d_list,ncol))),1)))
x
f[which(sapply(d_list,ncol)!=names(tail(sort(table(sapply(d_list,ncol))),1)))]
if(length(x)>0) d_list[[which(sapply(d_list,ncol)!=names(tail(sort(table(sapply(d_list,ncol))),1)))]]

for(i in 1:length(d_list)){d_list[[i]]$Total.OPV.doses.received.through.RI<-as.character(d_list[[i]]$Total.OPV.doses.received.through.RI)}
for(i in 1:length(d_list)){d_list[[i]]$Total.OPV.doses.received.through.SIA<-as.character(d_list[[i]]$Total.OPV.doses.received.through.SIA)}
for(i in 1:length(d_list)){d_list[[i]]$Total.IPV.doses.received.through.RI<-as.character(d_list[[i]]$Total.IPV.doses.received.through.RI)}
for(i in 1:length(d_list)){d_list[[i]]$Date.of.last.OPV.dose.received.through.SIA<-as.character(d_list[[i]]$Date.of.last.OPV.dose.received.through.SIA)}
for(i in 1:length(d_list)){d_list[[i]]$Source.of.RI.vaccination.information<-as.character(d_list[[i]]$Source.of.RI.vaccination.information)}
for(i in 1:length(d_list)){d_list[[i]]$Sex<-as.character(d_list[[i]]$Sex)}
for(i in 1:length(d_list)){d_list[[i]]$Date.of.Birth<-as.character(d_list[[i]]$Date.of.Birth)}

d<-do.call("bind_rows",d_list)

# Look for missing case EPIDs
d$case.epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",d$Case.EPID),fixed=T),fixed=T),fixed=T)
case_data<-unique(d[,c("case.epid","Case.EPID","Case.ADM1", "Case.ADM2", "Case.Ward", "Case.Onset", "Case.DOB", "Case.Age.Months","File.Name")])
case_data[case_data$case.epid=="",]

f<-"NIE-FCT-ABC-21-021.xlsx"
d$Case.EPID[d$File.Name==f]<-"NIE-FCT-ABC-21-012"

f<-"NIE-SOS-RBH-21-003.xls"
d$Case.EPID[d$File.Name==f]<-"NIE-SOS-RBH-21-003"

f<-"NIE-SOS-RBH-21-003C1.xls"
d$Case.EPID[d$File.Name==f]<-"NIE-SOS-RBH-21-003C1"

f<-"NIE-ZAS-MRD-22-018.xlsx"
d$Case.EPID[d$File.Name==f]<-"NIE-ZAS-MRD-22-018"

# Look for missing survey date
case_data<-unique(d[,c("Survey.Date","File.Name")])
case_data[case_data$Survey.Date=="",]
d$Survey.Date[d$File.Name=="NIE-BAU-DBM-21-007.xlsx"]<-"24/08/2021"
d$Survey.Date[d$File.Name=="NIE-KBS-DRD-23-007.xlsx"]<-"14/05/2023"
d$Survey.Date[d$File.Name=="NIE-KBS-MYM-21-016.xlsx"]<-"23/12/2021"
d$Survey.Date[d$File.Name=="NIE-KBS-WRR-21-013.xlsx"]<-"25/12/2021"
d$Survey.Date[d$File.Name=="NIE-KDS-ANC-21-004.xlsx"]<-"07/01/2022"
d$Survey.Date[d$File.Name=="NIE-KNS-AJG-21-003C2.xlsx"]<-"12/10/2021"
d$Survey.Date[d$File.Name=="NIE-KNS-GZW-21-003.xlsx"]<-"14/10/2021"
d$Survey.Date[d$File.Name=="NIE-NAS-NSW-22-003.xlsx"]<-"25-03-2022"
d$Survey.Date[d$File.Name=="NIE-SOS-GAD-21-015.xls"]<-"11/02/2022"
d$Survey.Date[d$File.Name=="NIE-SOS-RBH-21-003.xls"]<-"07/01/2022"
d$Survey.Date[d$File.Name=="NIE-SOS-RBH-21-003C1.xls"]<-"10/8/2021"

# No survey date: NIE-JIS-GGW-21-017C1.xlsx, NIE-JIS-GGW-21-018.xlsx, NIE-JIS-KGM -21-018.xlsx, NIE-KBS-BES-21-018C2-4.xlsx, NIE-KBS-BES-21-018C2.xlsx


# Look for duplicates in case EPIDs
d$case.epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",d$Case.EPID),fixed=T),fixed=T),fixed=T)
case_data<-unique(d[,c("case.epid","Case.EPID","Case.ADM1", "Case.ADM2", "Case.Ward", "Case.Onset", "Case.DOB", "Case.Age.Months","File.Name")])
s<-summarise(group_by(as_tibble(case_data),case.epid),n=n())
s[s$n>1,]

d_orig<-d

a<-unique(s$case.epid[s$n==2])
a
diff<-sapply(a,function(epid){
  sub<-case_data[case_data$case.epid==epid,]
  sub_a<-d[d$File.Name==sub$File.Name[1],]
  sub_b<-d[d$File.Name==sub$File.Name[2],]
  names<-intersect(names(sub_a),names(sub_b))
  names<-names[!names%in%c("File.Name","Case.EPID")]
  x<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    b<-as.character(sub_b[,n])
    a[is.na(a)|a==""]<-"XXXX"
    b[is.na(b)|b==""]<-"XXXX"
    return(all(a==b))
  })
  cols<-names(x[is.na(x)|!x])
  ret<-NULL
  if(length(cols)>0) {
    cols<-c(names(x[is.na(x)|!x]),"File.Name","Case.EPID")
    if(any(cols=="Age.month")) cols<-c(cols,"Case.Age.Months")
    ret<-cbind(sub_a[,cols],sub_b[,cols],epid=epid)
  }
  return(ret)
})

# Duplicate files
epids<-names(diff[which(sapply(diff,length)==0)])
files<-sapply(epids,function(e){unique(d$File.Name[d$case.epid==e])[1]})
dropped<-files
d<-d[!d$File.Name%in%files,]

# NIE-AESI-BAS-JMA-21-013 C3.xlsx, NIE-BAS-JMA-21-013C3-4.xlsx
# 2 controls differ in age 24 vs 12 months, 32 vs 12 months, case age = 36 months
diff[which(sapply(diff,length)>0)[1]]

# NIE-GMS-KWM-21-019-C3.xlsx   NIE-GMS-KWM-21-019C3-4.xlsx
# Date format
diff[which(sapply(diff,length)>0)[2]]

# NIE-JIS-DUT-21-058C3-4.xlsx  NIE-JIS-DUT-21-058C3.xlsx
# 7 controls differ in age
diff[which(sapply(diff,length)>0)[3]]
sum(diff[which(sapply(diff,length)>0)[3]][[1]][,2]!=diff[which(sapply(diff,length)>0)[3]][[1]][,7])
sum(abs(36-as.numeric(diff[which(sapply(diff,length)>0)[3]][[1]][,2]))>12)
sum(abs(36-as.numeric(diff[which(sapply(diff,length)>0)[3]][[1]][,7]))>12)

# NIE-JIS-GWW-21-002C3-4.xlsx  NIE-JIS-GWW-21-002C3.xlsx
# 8 controls differ in age
diff[which(sapply(diff,length)>0)[4]]
sum(diff[which(sapply(diff,length)>0)[4]][[1]][,2]!=diff[which(sapply(diff,length)>0)[4]][[1]][,7])
sum(abs(120-as.numeric(diff[which(sapply(diff,length)>0)[4]][[1]][,2]))>12)
sum(abs(120-as.numeric(diff[which(sapply(diff,length)>0)[4]][[1]][,7]))>12)

# NIE-KBS-BES-21-023-C1 2.xlsx  NIE-KBS-BES-21-023-C1.xlsx
# Date format
diff[which(sapply(diff,length)>0)[5]]

# NIE-KBS-MHT-21-016 2.xlsx, NIE-KBS-MHT-21-016.xlsx
# mislabeled case EPID?
diff[which(sapply(diff,length)>0)[6]]

# NIE-KBS-SNA-21-011C-4.xlsx   NIE-KBS-SNA-21-011C2.xlsx
# Date format
diff[which(sapply(diff,length)>0)[7]]

# NIE-KTS-BRE-21-003-4.xlsx  NIE-KTS-BRE-21-005-4.xlsx
# mislabeled case EPID?
diff[which(sapply(diff,length)>0)[8]]
d$Case.EPID[d$File.Name=="NIE-KTS-BRE-21-005-4.xlsx"]<-"NIE-KTS-BRE-21-005"
d$case.epid<-gsub("-","",gsub("_","",gsub("/","-",gsub(" ","",d$Case.EPID),fixed=T),fixed=T),fixed=T)

# NIE-TRS-GAS-21-012-C3 (2).xlsx  NIE-TRS-GAS-21-012-C3.xlsx
# 3 controls differ in age
diff[which(sapply(diff,length)>0)[9]]
sum(diff[which(sapply(diff,length)>0)[9]][[1]][,1]!=diff[which(sapply(diff,length)>0)[9]][[1]][,5])
sum(abs(30-as.numeric(diff[which(sapply(diff,length)>0)[9]][[1]][,1]))>12)
sum(abs(30-as.numeric(diff[which(sapply(diff,length)>0)[9]][[1]][,5]))>12)

# Drop duplicate survey files (16 Dec 21 email communication from Dr. Walla)
d<-d[!d$File.Name%in%c("NIE-AESI-BAS-JMA-21-013 C3.xlsx", # Different ages
                       "NIE-GMS-KWM-21-019-C3.xlsx", # Date format
                       "NIE-JIS-DUT-21-058C3.xlsx", # Different ages
                       "NIE-JIS-GWW-21-002C3.xlsx", # Different ages
                       "NIE-KBS-BES-21-023-C1.xlsx", # Date format
                       "NIE-TRS-GAS-21-012-C3.xlsx", # Different ages
                       "NIE-KBS-SNA-21-011C2.xlsx"),]  # Date format

dropped<-c(dropped,c("NIE-AESI-BAS-JMA-21-013 C3.xlsx", 
                     "NIE-GMS-KWM-21-019-C3.xlsx", 
                     "NIE-JIS-DUT-21-058C3.xlsx", 
                     "NIE-JIS-GWW-21-002C3.xlsx", 
                     "NIE-KBS-BES-21-023-C1.xlsx",
                     "NIE-TRS-GAS-21-012-C3.xlsx", 
                     "NIE-KBS-SNA-21-011C2.xlsx"))
# Triplicate files
a<-unique(s$case.epid[s$n==3])
a
drop<-sapply(a,function(epid){
  sub<-case_data[case_data$case.epid==epid,]
  sub_a<-d[d$File.Name==sub$File.Name[1],]
  sub_b<-d[d$File.Name==sub$File.Name[2],]
  sub_c<-d[d$File.Name==sub$File.Name[3],]
  names<-intersect(intersect(names(sub_a),names(sub_b)),names(sub_c))
  names<-names[!names%in%c("File.Name","Case.EPID")]
  x<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    b<-as.character(sub_b[,n])
    a[is.na(a)|a==""]<-"XXXX"
    b[is.na(b)|b==""]<-"XXXX"
    return(all(a==b))
  })
  y<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    c<-as.character(sub_c[,n])
    a[is.na(a)|a==""]<-"XXXX"
    c[is.na(c)|c==""]<-"XXXX"
    return(all(a==c))
  })
  z<-sapply(names,function(n){
    b<-as.character(sub_b[,n])
    c<-as.character(sub_c[,n])
    b[is.na(b)|b==""]<-"XXXX"
    c[is.na(c)|c==""]<-"XXXX"
    return(all(b==c))
  })
  # A&B identical
  ab<-all(x)
  # A&C identical
  ac<-all(y)
  # B&C identical
  bc<-all(z)
  if(all(c(ab,ac,bc))) drop<-sub$File.Name[1:2]
  if(ab&!ac&!bc) drop<-sub$File.Name[1]
  if(!ab&ac&!bc) drop<-sub$File.Name[1]
  if(!ab&!ac&bc) drop<-sub$File.Name[2]
  return(drop)})
# Duplicate files
dropped<-c(dropped,unlist(drop))
d<-d[!d$File.Name%in%unlist(drop),]

# Quadruplicate files
a<-unique(s$case.epid[s$n==4])
a

drop<-sapply(a,function(epid){
  sub<-case_data[case_data$case.epid==epid,]
  sub_a<-d[d$File.Name==sub$File.Name[1],]
  sub_b<-d[d$File.Name==sub$File.Name[2],]
  sub_c<-d[d$File.Name==sub$File.Name[3],]
  sub_d<-d[d$File.Name==sub$File.Name[4],]
  names<-intersect(intersect(intersect(names(sub_a),names(sub_b)),names(sub_c)),names(sub_d))
  names<-names[!names%in%c("File.Name","Case.EPID")]
  ab<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    b<-as.character(sub_b[,n])
    a[is.na(a)|a==""]<-"XXXX"
    b[is.na(b)|b==""]<-"XXXX"
    return(all(a==b))
  })
  ac<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    c<-as.character(sub_c[,n])
    a[is.na(a)|a==""]<-"XXXX"
    c[is.na(c)|c==""]<-"XXXX"
    return(all(a==c))
  })
  ad<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    d<-as.character(sub_d[,n])
    a[is.na(a)|a==""]<-"XXXX"
    d[is.na(d)|d==""]<-"XXXX"
    return(all(a==d))
  })
  bc<-sapply(names,function(n){
    b<-as.character(sub_b[,n])
    c<-as.character(sub_c[,n])
    b[is.na(b)|b==""]<-"XXXX"
    c[is.na(c)|c==""]<-"XXXX"
    return(all(b==c))
  })
  bd<-sapply(names,function(n){
    b<-as.character(sub_b[,n])
    d<-as.character(sub_d[,n])
    b[is.na(b)|b==""]<-"XXXX"
    d[is.na(d)|d==""]<-"XXXX"
    return(all(b==d))
  })
  cd<-sapply(names,function(n){
    d<-as.character(sub_d[,n])
    c<-as.character(sub_c[,n])
    d[is.na(d)|d==""]<-"XXXX"
    c[is.na(c)|c==""]<-"XXXX"
    return(all(d==c))
  })
  # A&B identical
  ab<-all(ab)
  # A&C identical
  ac<-all(ac)
  # A&D identical
  ad<-all(ad)
  # B&C identical
  bc<-all(bc)
  # B&D identical
  bd<-all(bd)
  # C&D identical
  cd<-all(cd)
  if(all(c(ab,ac,ad,bc,bd,cd))) drop<-sub$File.Name[1:3]
  if(ab&!ac& ad&!bc& bd&!cd) drop<-sub$File.Name[1:2]
  if(ab&!ac&!ad&!bc&!bd&!cd) drop<-sub$File.Name[1]
  return(drop)})
# Duplicate files
dropped<-c(dropped,unlist(drop))

d<-d[!d$File.Name%in%unlist(drop),]


# Repeat
case_data<-unique(d[,c("case.epid","Case.EPID","Case.ADM1", "Case.ADM2", "Case.Ward", "Case.Onset", "Case.DOB", "Case.Age.Months","File.Name")])
s<-summarise(group_by(as_tibble(case_data),case.epid),n=n())
s[s$n>1,]

a<-unique(s$case.epid[s$n==2])
a

diff<-sapply(a,function(epid){
  sub<-case_data[case_data$case.epid==epid,]
  sub_a<-d[d$File.Name==sub$File.Name[1],]
  sub_b<-d[d$File.Name==sub$File.Name[2],]
  names<-intersect(names(sub_a),names(sub_b))
  names<-names[!names%in%c("File.Name","Case.EPID")]
  x<-sapply(names,function(n){
    a<-as.character(sub_a[,n])
    b<-as.character(sub_b[,n])
    a[is.na(a)|a==""]<-"XXXX"
    b[is.na(b)|b==""]<-"XXXX"
    return(all(a==b))
  })
  cols<-names(x[is.na(x)|!x])
  ret<-NULL
  if(length(cols)>0) {
    cols<-c(names(x[is.na(x)|!x]),"File.Name","Case.EPID")
    if(any(cols=="Age.month")) cols<-c(cols,"Case.Age.Months")
    ret<-cbind(sub_a[,cols],sub_b[,cols],EPID=epid)
  }
  return(ret)
})

# Duplicate files
epids<-names(diff[which(sapply(diff,length)==0)])
files<-sapply(epids,function(e){unique(d$File.Name[d$case.epid==e])[1]})
dropped<-c(dropped,files)
d<-d[!d$File.Name%in%files,]

# NIE-KBS-MHT-21-016 2.xlsx, NIE-KBS-MHT-21-016.xlsx 
# mislabeled case EPID?
diff[which(sapply(diff,length)>0)[1]]

# NIE-SOS-BDN-21-002 (2).xlsx  NIE-SOS-BDN-21-002.xlsx
# Case age 36/36?
diff[which(sapply(diff,length)>0)[2]]

# NIE-SOS-WMK-21-002C2-4.xlsx NIE-SOS-WMK-21-002C2.xlsx 
# Date format
diff[which(sapply(diff,length)>0)[3]]

d<-d[!d$File.Name%in%c("NIE-SOS-BDN-21-002.xlsx", # Case age format
                       "NIE-SOS-WMK-21-002C2.xlsx", # Date format
                       "NIE-KBS-MHT-21-016 2.xlsx"),]  # Communication from Dr. Walla 16 Dec 21

dropped<-c(dropped,c("NIE-SOS-BDN-21-002.xlsx", 
                     "NIE-SOS-WMK-21-002C2.xlsx", 
                     "NIE-KBS-MHT-21-016 2.xlsx"))

# Look for duplicates
case_data<-unique(d[,c("case.epid","Case.EPID","Case.ADM1", "Case.ADM2", "Case.Ward", "Case.Onset", "Case.DOB", "Case.Age.Months","File.Name")])
s<-summarise(group_by(as_tibble(case_data),case.epid),n=n())
s[s$n>1,]


length(unique(d$File.Name))
length(unique(d_orig$File.Name))
length(unlist(dropped))
write.csv(d,paste0(subdir,"combined.csv"),row.names=F)



