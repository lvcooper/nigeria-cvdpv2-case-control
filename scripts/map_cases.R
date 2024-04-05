load("data/nigeria_maps.Rdata")

define.distance<-50000
define.age<-12
define.onset<-30
match.on.RI<-F
if(file.exists("data/all.csv")){a<-read.csv("data/all.csv")} else{return(warning("Need to first run read_data.R"))}
if(file.exists("data/eligible.csv")){d<-read.csv("data/eligible.csv")} else{return(warning("Need to first run read_data.R"))}
filename<-paste0("data/matched_distance_",define.distance,"_age_",define.age,"_onset_",define.onset,ifelse(match.on.RI,"_RI",""),".csv")
if(file.exists(filename)){m1<-read.csv(filename)} else{return(warning("Need to first run matching.R"))}

map<-merge(adm2data,summarise(group_by(as_tibble(m1),GUID),matched=sum(case)),all.x=T)
map<-merge(map,summarise(group_by(as_tibble(a),GUID),full_data=sum(grepl("Case",class))),all.x=T)
map<-merge(map,summarise(group_by(as_tibble(d),GUID),eligible=sum(grepl("Case",class))),all.x=T)
map$matched[is.na(map$matched)]<-0
map$eligible[is.na(map$eligible)]<-0
map$full_data[is.na(map$full_data)]<-0
map$p_matched<-map$matched/map$full_data

x=0.2
p1<-ggplot()+geom_sf(data=adm0data,fill="dark gray",linetype=0)+
  geom_sf(data=map,aes(fill=as.factor(full_data)),linetype=0)+
  geom_sf(data=map[map$full_data>0,],fill="transparent",size=x)+
  geom_sf(data=adm1data,fill="transparent",size=3*x,color="black")+
  scale_fill_manual(values=c("light gray",viridis(n=length(unique(map$full_data))-1)))+labs(fill="Cases")+theme_void()

map<-merge(adm1data,summarise(group_by(as_tibble(m1),adm1),matched=sum(case)),all.x=T)
map<-merge(map,summarise(group_by(as_tibble(a),adm1),full_data=sum(grepl("Case",class))),all.x=T)
map<-merge(map,summarise(group_by(as_tibble(d),adm1),eligible=sum(grepl("Case",class))),all.x=T)
map$matched[is.na(map$matched)]<-0
map$eligible[is.na(map$eligible)]<-0
map$full_data[is.na(map$full_data)]<-0
map$p_matched<-map$matched/map$full_data
map$p_eligible<-map$eligible/map$full_data

map$CENTER_LON[map$adm1=="kebbi"]<-4
map$CENTER_LAT[map$adm1=="lagos"]<-6.2
map$CENTER_LAT[map$adm1=="delta"]<-5.4
map$CENTER_LAT[map$adm1=="kwara"]<-9.2
map$CENTER_LAT[map$adm1=="abia"]<-5.5
p2<-ggplot()+geom_sf(data=map,aes(fill=p_matched),color="black")+geom_text(data=map[map$full_data>0&map$p_matched>0.3,],aes(x=CENTER_LON,y=CENTER_LAT,label=full_data),size=3)+
geom_text(data=map[map$full_data>0&map$p_matched<=0.3,],aes(x=CENTER_LON,y=CENTER_LAT,label=full_data),color="white",size=3)+
  scale_fill_viridis_c(label=scales::percent)+labs(fill="Matched",x="",y="")+theme_void()

if(!dir.exists("figs")) dir.create("figs")
ggsave(paste0("figs/map_cases_match_",define.distance,".png"),ggarrange(p1,p2,labels="auto",widths=c(1,1.07)),width=10,height=5)
p3<-ggarrange(p1,p2,labels="auto",widths=c(1,1.07))

define.offset<-28
use.date.of.last.opv<-T
filename<-paste0("data/matched_distance_",define.distance,
                 "_age_",define.age,
                 "_onset_",define.onset,
                 ifelse(match.on.RI,"_RI",""),
                 "_offset_",define.offset,ifelse(use.date.of.last.opv,"","_no_date"),".csv")
if(file.exists(filename)){m1<-read.csv(filename)} else{return(warning("Need to first run sias.R"))}
m1$exp<-paste(ifelse(m1$nOPV2_SIAs_ons>0,"nOPV2",""),ifelse(m1$mOPV2_SIAs_ons,"mOPV2",""))
m1$exp<-factor(m1$exp,labels=c("Neither","mOPV2 SIAs","nOPV2 SIAs","Both"))
table(m1$exp,m1$case)
map<-merge(adm2data,m1[m1$case,],all.y=T)
map<-cbind(map,t(sapply(1:nrow(map),function(i){x<-st_sample(map[i,],1);return(st_bbox(x)[1:2])})))
p1<-ggplot()+geom_sf(data=adm0data,fill="white",linetype=0)+
  geom_point(data=map,aes(x=xmin,y=ymin,color=exp),size=1.5)+
  geom_sf(data=adm1data,fill="transparent",size=0.5)+labs(color="Exposed to")+theme_void()+guides(color="none")

m1$month<-floor_date(as.Date(m1$date),"month")
a$month<-floor_date(as.Date(a$date),"month")
x<-summarise(group_by(as_tibble(m1),month,case,exp),n=n())%>%filter(case)
y<-merge(summarise(group_by(as_tibble(a),month),full_data=sum(grepl("Case",class))),
      summarise(group_by(as_tibble(m1),month),matched=sum(case)),all.x=T)
y$matched[is.na(y$matched)]<-0
y$n<-y$full_data-y$matched
y<-y[!is.na(y$month)&y$n>0,]
x$exp<-as.character(x$exp)
x<-bind_rows(x,cbind(y,case=TRUE,exp="Unmatched"))
x$exp<-factor(x$exp,levels=c("Unmatched","Neither","mOPV2 SIAs","nOPV2 SIAs","Both"),labels=c("Unmatched","Exposed to neither","Exposed to mOPV2 SIAs","Exposed to nOPV2 SIAs","Exposed to nOPV2 & mOPV2 SIAs"),ordered=T)
p2<-x%>%ggplot()+geom_col(aes(x=month,y=n,fill=exp))+
  labs(x="",y="Cases",fill="Category")+theme(legend.position=c(0.3,0.7))+scale_fill_manual(values=c("dark grey",scales::hue_pal()(4)))
p<-ggarrange(p2,p1,labels="auto")
ggsave(paste0("figs/exposure_map_",define.distance,".png"),p,width=10,height=4)

p2<-x%>%filter(exp!="Unmatched")%>%ggplot()+geom_col(aes(x=month,y=n,fill=exp))+
  labs(x="",y="Cases",fill="Category")+theme_minimal()+theme(legend.position=c(0.3,0.7))+scale_fill_manual(values=c(scales::hue_pal()(4)))
p<-ggarrange(p2,p1,labels=c("c","d"),widths=c(1.3,1))
ggsave(paste0("figs/multi_panel.png"),ggarrange(p3,p,nrow=2),width=10,height=8)
