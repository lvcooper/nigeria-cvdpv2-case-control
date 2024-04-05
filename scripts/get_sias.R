set_token("token.txt")
sia<-get_polis_subactivities(country_iso3="NGA",min_date="2016-05-01")

# Subset SIAs which were done
sia<-sia[sia$status%in%c("Done"),]
sia$date<-as.Date(sia$activity_parent_planned_date_from)
unique(sia[is.na(sia$date),c("activity_parent_planned_date_from","sia_sub_activity_code")])
sia<-sia[which(sia$date>=as.Date("2016-05-01")),]

# Age group
key<-cbind.data.frame(age_group=c("0 to 5 years","2-35M","2-59M","4-59M","6-59M","18M-47M","14W-59M","0 to 9 years","14W-35M",
                                  "4-23M","2-24M","6W-59M","0-23M","4-5Y","0 to 15 years","All ages","0-14Y","0-11M","9-59M","0-71M","1-4Y","0-10Y"),
                      age_min=c(0,2,2,4,6,18,0,0,0,4,2,0,0,4*12,0,0,0,0,9,0,12,0),
                      age_max=c(5*12,35,59,59,59,47,59,9*12,35,23,24,59,23,5*12,15*12,15*12,14*12,11,59,71,47,10*12))
sia<-merge(sia,key,all.x=T)
if(any(is.na(sia$age_min))) warning(paste0("New SIA age group(s): ",paste0(unique(sia$age_group[is.na(sia$age_min)]),collapse=", ")))
write.csv(sia,"data/SIAs_Nigeria_case_control.csv")
rm(list=ls())