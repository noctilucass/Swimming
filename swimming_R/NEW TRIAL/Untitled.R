library(readr)
library(tidyverse)
library(rlist)
library(data.table) 
library(tibble)
library(lubridate)
library(chron)
options(digits.secs=3)





setwd("~/Desktop/ECOLAB!/paper tesis/swimming R/NEW TRIAL/")
list_of_files <- list.files(path = ".", recursive = TRUE,
                            pattern = "\\.txt$", 
                            full.names = TRUE)

#Lapply permite importar todos los archivos .txt de nado y apilarlos en una lista
loop_full<-lapply(list_of_files, read.delim2)
names(loop_full)<-list_of_files
loop_full<-lapply(loop_full, function(x) { x["frame"] <- NULL; x })


t1_12<-list()
mean_t1_12<-list()
for (i in 1:25){ a=loop_full[[i]][,]
a$x<-as.numeric(a$x)
a$y<-as.numeric(a$y)
df<-a %>%
  group_by(X) %>%
  mutate(XX = x - lag(x)) %>%
  mutate(Y = y - lag(y)) 
df$time<-as.POSIXct(df$time, format="%M:%OS")
df$tdif<-ifelse(difftime(df$time,lag(df$time), units = "secs")<1,difftime(df$time,lag(df$time), units = "secs"),NA)
#z1 <- df %>% group_by(X) %>% count(X)
z1<-aggregate(df$tdif, by=list(Trayectoria=df$X), FUN=sum, na.rm=T)
df$iv <- (sqrt((df$XX)^2+(df$Y)^2))/0.125
z1$IV <- aggregate(df$iv, list(df$X), FUN = mean, na.rm=TRUE)[,2]
mean_t1_12[[i]]<-z1
t1_12[[i]]=df }

names(t1_12)<-list_of_files[1:25]
t1_12<-Map(cbind, t1_12, group = names(t1_12))
t1_12<- list.rbind(t1_12)

names(mean_t1_12)<-list_of_files[1:25]
mean_t1_12<-Map(cbind, mean_t1_12, group = names(mean_t1_12))
mean_t1_12<- list.rbind(mean_t1_12)
mean_t1_12<-mean_t1_12 %>% separate(group, c("A", "B","C","D","e","f","g","h"))
mean_t1_12<-mean_t1_12[,c(9,10,11,1,2,3,4)]
colnames(mean_t1_12)<-c("Temp","Pote","Larva","Trayectoria","Time","IV")

summary(t1_12$iv)

write_excel_csv(mean_t1_12,"mean_t1_12.csv" ,delim = ",")


t1_15<-list()
mean_t1_15<-list()
for (i in 26:49){ 
a=loop_full[[i]][,]
a$x<-as.numeric(a$x)
a$y<-as.numeric(a$y)
df<-a %>%
  group_by(X) %>%
  mutate(XX = x - lag(x)) %>%
  mutate(Y = y - lag(y))
df$time<-as.POSIXct(df$time, format="%M:%OS")
df$tdif<-ifelse(difftime(df$time,lag(df$time), units = "secs")<1,difftime(df$time,lag(df$time), units = "secs"),NA)
df$iv <- (sqrt((df$XX)^2+(df$Y)^2))/0.125
z1<-aggregate(df$tdif, by=list(Trayectoria=df$X), FUN=sum, na.rm=T)
z1$IV <- aggregate(df$iv, list(df$X), FUN = mean, na.rm=TRUE)[,2]
mean_t1_15[[i]]<-z1
t1_15[[i]]=df }

t1_15<-t1_15[c(26:49)]
names(t1_15)<-list_of_files[26:49]
t1_15<-Map(cbind, t1_15, group = names(t1_15))
t1_15<- list.rbind(t1_15)

mean_t1_15<-mean_t1_15[c(26:49)]
names(mean_t1_15) <- list_of_files[26:49]
mean_t1_15<-Map(cbind, mean_t1_15, group = names(mean_t1_15))
mean_t1_15<- list.rbind(mean_t1_15)

mean_t1_15<-mean_t1_15 %>% separate(group, c("A", "B","C","D","e","f","g"))
mean_t1_15<-mean_t1_15[,c(8,9,10,1,2,3)]
colnames(mean_t1_15)<-c("Temp","Pote","Larva","Trayectoria","Time","IV")

t1_15$iv<-na_if(t1_15$iv, 0)
summary(t1_15$iv)

write_excel_csv(mean_t1_15,"mean_t1_15.csv" ,delim = ",")

#############################
t2_15<-list()
mean_t2_15<-list()
for (i in 50:58){ 
  a=loop_full[[i]][,]
  a$x<-as.numeric(a$x)
  a$y<-as.numeric(a$y)
  df<-a %>%
    group_by(X) %>%
    mutate(XX = x - lag(x)) %>%
    mutate(Y = y - lag(y))
  df$iv <- (sqrt((df$XX)^2+(df$Y)^2))/0.125
  mean_t2_15[[i]]<-aggregate(df$iv, list(df$X), FUN = mean, na.rm=TRUE)
  t2_15[[i]]=df }
t2_15<-t2_15[c(50:58)]
names(t2_15)<-list_of_files[50:58]
t2_15<- list.rbind(t2_15)
names(mean_t2_15)<-list_of_files[50:58]
mean_t2_15<- list.rbind(mean_t2_15)
mean_t2_15$group<-rep("15")
t2_15$iv<-na_if(t2_15$iv, 0)
summary(t2_15$iv)
#############################

t1_17<-list()
mean_t1_17<-list()
for (i in 59:85){ a=loop_full[[i]][,]
a$x<-as.numeric(a$x)
a$y<-as.numeric(a$y)
df<-a %>%
  group_by(X) %>%
  mutate(XX = x - lag(x)) %>%
  mutate(Y = y - lag(y))
df$time<-as.POSIXct(df$time, format="%M:%OS")
df$tdif<-ifelse(difftime(df$time,lag(df$time), units = "secs")<1,difftime(df$time,lag(df$time), units = "secs"),NA)
df$iv <- (sqrt((df$XX)^2+(df$Y)^2))/0.125
z1<-aggregate(df$tdif, by=list(Trayectoria=df$X), FUN=sum, na.rm=T)
z1$IV <- aggregate(df$iv, list(df$X), FUN = mean, na.rm=TRUE)[,2]
mean_t1_17[[i]]<-z1
t1_17[[i]]=df }

t1_17<-t1_17[c(59:85)]
names(t1_17)<-list_of_files[59:85]
t1_17<- list.rbind(t1_17)


mean_t1_17<-mean_t1_17[c(59:85)]
names(mean_t1_17) <- list_of_files[59:85]
mean_t1_17<-Map(cbind, mean_t1_17, group = names(mean_t1_17))
mean_t1_17<- list.rbind(mean_t1_17)

mean_t1_17<-mean_t1_17 %>% separate(group, c("A", "B","C","D","e","f","g"))
mean_t1_17<-mean_t1_17[,c(8,9,10,1,2,3)]
colnames(mean_t1_17)<-c("Temp","Pote","Larva","Trayectoria","Time","IV")


t1_17$iv<-na_if(t1_17$iv, 0)
summary(t1_17$iv)

write_excel_csv(mean_t1_17,"mean_t1_17.csv" ,delim = ",")

ggplot(t1_12, aes((iv))) + geom_density() + xlab("Instant velocity (mm/s)")
ggplot(t1_15, aes((iv))) + geom_density() + xlab("Instant velocity (mm/s)")
ggplot(t1_17, aes((iv))) + geom_density() + xlab("Instant velocity (mm/s)")
ggplot(t2_15, aes((iv))) + geom_density() + xlab("Instant velocity (mm/s)")

full_mean<-rbind(mean_t1_12,mean_t1_15,mean_t1_17)
write_excel_csv(full_mean,"full_mean.csv" ,delim = ",")

ggplot(full_mean,aes(x=Temp, y=IV)) + geom_boxplot()

write.csv(A_12_1,"~/Desktop/ECOLAB!/paper tesis/swimming/grabación1_12/Pote A/A_12_1.csv", row.names = FALSE)



##########
#####Calculo de IV
# calculo de tiempo en una columna, agregar group_by para que sea por trayectoria
a %>% group_by(X) %>% mutate(tdif=difftime(a$time,lag(a$time), units = "secs"))
#ifelse (difftime(a$time,lag(a$time), units = "secs")>1,,NA)
ave(a$time, a$X, FUN=function(x) difftime(a$time,lag(a$time)))

sum(df$tdif,na.rm=TRUE)


