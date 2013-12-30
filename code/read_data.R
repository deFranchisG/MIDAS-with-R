rm(list=ls())
setwd("/Users/carlo/Dropbox/R/Midas/donnees")
library(ggplot2)
library(zoo)
library(reshape2)
library(zoo)
library(lubridate)
library(chron)


cac_40<-na.omit(read.csv2("base_j_2.csv",sep=";",dec=",",na.strings=c("#N/A")))[,-3]
colnames(cac_40)<-c("Date","cac_40")
cac_40$Date<-as.Date(cac_40$Date,format="%d/%m/%Y")

dow_jones<-na.omit(read.csv2("base_j_2.csv",sep=";",dec=",",na.strings=c("#N/A")))[,-2]
colnames(dow_jones)<-c("Date","dow_jones")
dow_jones$Date<-as.Date(dow_jones$Date,format="%d/%m/%Y")

change<-na.omit(read.csv2("base_j.csv",sep=";",dec=",",na.strings=c("#N/A")))[,-3]
colnames(change)<-c("Date","change")
change$Date<-as.Date(change$Date,format="%d/%m/%y")

volatilite<-na.omit(read.csv2("base_j.csv",sep=";",dec=",",na.strings=c("#N/A")))[,-2]
colnames(volatilite)<-c("Date","volatilite")
volatilite$Date<-as.Date(volatilite$Date,format="%d/%m/%y")

tmp<-na.omit(read.csv2("base_j_3.csv",sep=",",dec=",",na.strings=c("#N/A")))
colnames(tmp)<-c("Date","tmp")
tmp$Date<-as.Date(tmp$Date,format="%d/%m/%Y")

vix<-na.omit(read.csv2("base_j_5.csv",sep=",",dec=",",na.strings=c("#N/A")))
colnames(vix)<-c("Date","vix")
vix$Date<-as.Date(vix$Date,format="%d/%m/%y")

piboff<-na.omit(read.csv2("base_j_4.csv",dec="."))
colnames(piboff)<-c("Date","piboff")
piboff$Date<-as.Date(piboff$Date,format="%d/%m/%Y")



climat<-na.omit(read.csv2("base_m.csv",sep=";",dec=",",na.strings=c("#N/A")))
colnames(climat)<-c("Date","climat")
climat$Date<-gsub("/","/01/",climat$Date)
climat$Date<-as.Date(climat$Date,format="%m/%d/%y")



data<-merge(cac_40,tmp,by.x="Date",by.y="Date")
data<-merge(data,vix,by.x="Date",by.y="Date")
data<-merge(data,dow_jones,by.x="Date",by.y="Date")
data<-merge(data,climat,by.x="Date",by.y="Date",all.x=TRUE,all.y=TRUE)
data<-merge(data,piboff,by.x="Date",by.y="Date",all.x=TRUE,all.y=TRUE)
data<-merge(data,change,by.x="Date",by.y="Date",all.x=TRUE,all.y=TRUE)


base_t<-read.csv2("base_t.csv",header=TRUE,sep="\t",dec=",",allowEscapes = TRUE)
base_t[,2]<-log(base_t[,2])
base_t[(2:nrow(base_t)),2]<-base_t[(2:nrow(base_t)),2]-base_t[(1:nrow(base_t)-1),2]

base_t[,1]<-gsub("Q1","-03-01",base_t[,1])
base_t[,1]<-gsub("Q2","-06-01",base_t[,1])
base_t[,1]<-gsub("Q3","-09-01",base_t[,1])
base_t[,1]<-gsub("Q4","-12-01",base_t[,1])
base_t[,1]<-as.Date(base_t[,1],format="%Y-%m-%d")


data<-merge(base_t,data,by.x="date",by.y="Date",all.x=TRUE,all.y=TRUE)
colnames(data)[1]<-"Date"

data$dow_jones<-c(NA,diff(data$dow_jones))
data$change<-c(NA,diff(data$change))
data$cac_40<-c(NA,diff(data$cac_40))
data$spread<-data$piboff-data$tmp
data$spread<-c(NA,diff(data$spread))



save(list = ls(all=TRUE), file = "midas_data.RData")


ts.plot(data$vix)
ts.plot(data$cac_40)

ts.plot(data$dow_jones)
ts.plot(data$tmp)
ts.plot(data$piboff)
ts.plot(data$spread)


