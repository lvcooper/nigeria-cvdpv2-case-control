source("scripts/libraries.R")
source("scripts/read_data.R")
source("scripts/matching.R")

# Data for Figure 1
x<-read.csv("data/flow.csv")
x
rev(diff(rev(x$Case)))
round(100*(1-x$Case[2:10]/x$Case[1:9]))
rev(diff(rev(x$Control)))
round(100*(1-x$Control[2:10]/x$Control[1:9]))

# Match baseline parameters (50 km, 12 months, 30 days)
matching(define.distance=50000,define.age=12,define.onset=30)
read.csv("data/flow_match_distance_50000_age_12_onset_30.csv")


# Get SIAs from POLIS API
if(!file.exists("data/SIAs_Nigeria_case_control.csv")) source("scripts/get_sias.R")
source("scripts/sias.R")
sias(define.distance=50000,define.age=12,define.onset=30,define.offset=28,use.date.of.last.opv=T)

# Figure 2
source("scripts/map_cases.R")

# Figure 3
source("scripts/results_figs.R")

# # Analysis of reported date of last OPV SIA
# temp<-read.csv(paste0("data/matched_distance_50000_age_12_onset_30_offset_28.csv"))
# 
# x<-table(ifelse(is.na(temp$last_opv_date),"Date of last OPV dose not reported","Date of last OPV dose reported"),ifelse(temp$case,"Case","Control"))
# x
# paste0(x, " (",round(100*x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x)))),"%)")
# 
# temp$date_used<-temp$reported_recent_SIA_type!="Inaccurate date"
# x<-table(ifelse(!temp$date_used[!is.na(temp$last_opv_date)],"Date of last OPV reported, but > 30 days from SIA calendar",
#                 "Date of last OPV reported within 30 days from SIA calendar"),ifelse(temp$case[!is.na(temp$last_opv_date)],"Case","Control"))
# x
# paste0(x, " (",round(100*x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x)))),"%)")
# 
# table((temp$reported_recent_SIA_code==temp$most_recent_SIA_code)[!is.na(temp$last_opv_date)&temp$date_used&temp$case],temp$most_recent_SIA_type[!is.na(temp$last_opv_date)&temp$date_used&temp$case])
# table((temp$reported_recent_SIA_code==temp$most_recent_SIA_code)[!is.na(temp$last_opv_date)&temp$date_used&!temp$case],temp$most_recent_SIA_type[!is.na(temp$last_opv_date)&temp$date_used&!temp$case])
# 
# temp$after<-(as.Date(as.numeric(temp$reported_recent_SIA_date),origin=as.Date("1970-01-01")))>(as.Date(temp$date)-28)
# # x<-table(ifelse(temp$after[!is.na(temp$last_opv_date)&temp$date_used],"Date reported after onset","Date reported before onset"),
# #          ifelse(temp$case[!is.na(temp$last_opv_date)&temp$date_used],"Case","Control"))
# # x
# # x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x)))
# 
# x<-table(paste0("Date of last OPV reported, consistent with ",
#                 temp$reported_recent_SIA_type[!is.na(temp$last_opv_date)&temp$date_used]),ifelse(temp$case[!is.na(temp$last_opv_date)&temp$date_used],"Case","Control"))
# x
# x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x)))
# 
# y<-table(paste0("Date of last OPV reported before onset, consistent with ",temp$reported_recent_SIA_type[!is.na(temp$last_opv_date)&temp$date_used&!temp$after]),
#          ifelse(temp$case[!is.na(temp$last_opv_date)&temp$date_used&!temp$after],"Case","Control"))
# y
# y/x
# 
# # x<-table(paste0("Date of last OPV reported, consistent with ",
# #                 temp$type[!is.na(temp$last_opv_date)&temp$date_used&!temp$case]),
# #          ifelse(temp$after[!is.na(temp$last_opv_date)&temp$date_used&!temp$case],"After onset","Before onset"))
# # x
# temp$cat<-NA
# temp$cat[is.na(temp$last_opv_date)&is.na(temp$cat)]<-"No date"
# temp$cat[!temp$date_used&is.na(temp$cat)]<-"Inaccurate date"
# temp$cat[temp$reported_recent_SIA_type=="mOPV2"&temp$after&is.na(temp$cat)]<-"Received mOPV2 after"
# temp$cat[temp$reported_recent_SIA_type=="mOPV2"&!temp$after&is.na(temp$cat)]<-"Received mOPV2 before"
# temp$cat[temp$reported_recent_SIA_type=="bOPV"&temp$after&is.na(temp$cat)]<-"Received bOPV after"
# temp$cat[temp$reported_recent_SIA_type=="bOPV"&!temp$after&is.na(temp$cat)]<-"Received bOPV before"
# temp$cat[temp$reported_recent_SIA_type=="nOPV2"&temp$after&is.na(temp$cat)]<-"Received nOPV2 after"
# temp$cat[temp$reported_recent_SIA_type=="nOPV2"&!temp$after&is.na(temp$cat)]<-"Received nOPV2 before"
# x<-table(temp$cat,temp$class,useNA="ifany")
# x
# round(100*x/cbind(rep(colSums(x)[1],nrow(x)),rep(colSums(x)[2],nrow(x))))

source("scripts/regression.R")
# Main analysis (Table 1)
results<-regression(define.distance=50000,define.age=12,define.onset=30,define.offset=28,use.date.of.last.opv=T,IPV=T)
summary(results$model)
View(cbind(results$table1[c(5,4,3,7,9,8,10,2,6,1),],
           round(results$table1[c(5,4,3,7,9,8,10,2,6,1),2:5],2),
           results$table1[c(5,4,3,7,9,8,10,2,6,1),6:7]))
View(results$table2[c(13,14,1:4,8:12,5:7),c(1,2,5,3,6,4)])
x<-summary(results$model)$coef
rownames(x[1:2,])
2*(1-pnorm(abs(x[1,1]-x[2,1])/sqrt(x[2,3]^2 + x[1,3]^2 )))



# Healthy community control analysis
# Combine individual community control surveys
source("scripts/combine_control_surveys.R")
# Clean data
source("scripts/clean_control_survey_data_new.R")
# Exclusion criteria and match to cases
source("scripts/match_healthy_control.R")
# SIAs exposed
source("scripts/sias_healthy_control.R")
# Case-control analysis
source("scripts/regression_healthy_control.R")
summary(results$model)


# Compare results with healthy community controls versus non-polio AFP controls
source("scripts/compare_hc_npafp.R")

# Sensitivity analyses
# Bootstrap dose valency
source("scripts/boot_sens_analysis.R")
# Other sensitivity analyses
source("scripts/sensitivity_analyses.R")


# Simulate case control analysis to test methodology, explore role of inaccurate dose recall
source("scripts/simulate_case_control.R")
