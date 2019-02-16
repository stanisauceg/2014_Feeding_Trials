library(plyr)
library(ggplot2)

scaph<-read.csv("scaph-leftover-sizes.csv",header=TRUE)
head(scaph)

macscaph<-subset(scaph,scaph$sal.sp=='AMMA')
opscaph<-subset(scaph,scaph$sal.sp=='AMOP')
  opscaph1<-subset(opscaph,opscaph$rep==1)
  opscaph2<-subset(opscaph,opscaph$rep==2)
  opscaph3<-subset(opscaph,opscaph$rep==3)
  opscaph4<-subset(opscaph,opscaph$rep==4)
  opscaph5<-subset(opscaph,opscaph$rep==5)
  
# plot histograms w/ density plots

  # all scaph
  hist(scaph$size, breaks=15, prob=TRUE)
  lines(density(scaph$size))
  
  hist(scaph$size, breaks=30, prob=TRUE)
  lines(density(scaph$size))
  abline(v=350,lty=3)
  abline(v=400,lty=3)
  min(scaph$size)
  boxplot(scaph$size)
  
 # using ggplot2
qplot(scaph$size, geom="histogram")
qplot(scaph$size, geom="density")
m<-ggplot(data=macscaph,aes(macscaph$size))
m+geom_histogram()
+n+geom_histogram()
n<-ggplot(data=opscaph,aes(opscaph$size))
n+geom_histogram()

p<-ggplot(data=scaph,aes(scaph$size,fill=scaph$sal.sp))+geom_histogram(alpha=0.2)
p+geom_vline(xintercept=375)+geom_vline(xintercept=380)

p<-ggplot(data=scaph,aes(y=..density..,scaph$size,fill=scaph$sal.sp))+geom_histogram(alpha=0.2)
p+geom_vline(xintercept=375)+
  labs(x="Scapholeberis size",y="Density")+
  scale_fill_discrete(name="Salamander species",
                      breaks=c("AMOP", "AMMA"),
                      labels=c("A. opacum","A. maculatum"))


#alternate method w/o ggplot2
p1<-hist(macscaph$size,breaks=seq(from=0,to=2150,by=50))
p2<-hist(opscaph$size,breaks=seq(from=0,to=2150,by=50))
plot(p2,col=rgb(0,0,1,1/4),xlim=c(0,2000))
plot(p1,col=rgb(0,0,0,1/4),xlim=c(0,2000),add=TRUE)
abline(v=350,lty=3)

p3<-hist(opscaph1$size,breaks=seq(from=0,to=2150,by=50))
p4<-hist(opscaph2$size,breaks=seq(from=0,to=2150,by=50))
p5<-hist(opscaph3$size,breaks=seq(from=0,to=2150,by=50))
p6<-hist(opscaph4$size,breaks=seq(from=0,to=2150,by=50))
p7<-hist(opscaph5$size,breaks=seq(from=0,to=2150,by=50))

plot(p7,col=rgb(1,1,0,1/4),xlim=c(0,2200),ylim=c(0,40))
plot(p6,col=rgb(1,0,0,1/4),add=TRUE)
plot(p5,col=rgb(0,1,1,1/4),add=TRUE)
plot(p4,col=rgb(0,0,1,1/4),add=TRUE)
plot(p3,col=rgb(0,1,0,1/4),add=TRUE)
plot(p1,col=rgb(1,0,0,1/4),add=TRUE)
abline(v=350,lty=3)

p8<-hist(scaph$size,breaks=seq(from=0,to=2150,by=50))
abline(v=350,lty=3)

# use histogram to estimate and define size cutoff
scaph.cut<-subset(scaph,scaph$size>375)
head(scaph.cut)
summary(scaph.cut)
class(scaph.cut$treat)
# new counts after cutoff
    # create empty data frame
counts<-data.frame(matrix(data=NA, nrow=(5*5+9),ncol=5))
    # re-label columns

counts<-rename(counts,c("X1"="sal.sp","X2"="rep","X3"="sal.treat",
                        "X4"="treat","X5"="scaph_end"))
    # populate with treatment info
counts[,1]<-c(rep("AMOP",25),rep("AMMA",9))
counts[,2]<-c(rep(1:6,each=5),rep(6,4))
counts[,3]<-c(rep(c("ctrl",rep("AMOP",4)),5),"ctrl",rep("AMMA.uncued",4),rep("AMMA.cued",4))
counts[,4]<-c(rep(c("ctrl","avg","T1","T2","T3"),5),"ctrl",rep(c("avg","T1","T2","T3"),2))

counts

#rename sal.treat to correspond w/ better naming convention
head(scaph.cut)
summary(scaph.cut)

#scaph.cut$sal.treat[scaph.cut$sal.treat==AsFactor("present")]<-"AMOP"
#scaph.cut$sal.treat[scaph.cut$sal.treat=="uncued"]<-"AMMA.uncued"
#scaph.cut$sal.treat[scaph.cut$sal.treat=="cued"]<-"AMMA.cued"

attach(scaph.cut)
    
# what are actual counts
table(sal.treat,treat,rep)
        # this shortcut will work if there are no actual 0's
count.vector1<-c(t(table(sal.treat,treat,rep=1)))
count.vector<-count.vector[count.vector!=0]

    # populate with counts
counts[,3]<-c(70,26,7,53,75,41,3,69,35)
        # or, using shortcut
counts.alt<-counts
counts.alt[,3]<-count.vector

# ta-da!
counts
counts.alt


# # # # # # # # # # # #
# Analyze actual data #
# # # # # # # # # # # #

# create new table with only the scaph that exceed threshold size

library(vcd)

tapply(scaph.cut$size,list(scaph.cut$treat, scaph.cut$rep),length)
c(tapply(scaph.cut$size,list(scaph.cut$treat, scaph.cut$rep),length))

counts[,5]<-c(c(tapply(scaph.cut$size,list(scaph.cut$treat, scaph.cut$rep),length)),rep("na",4))
counts[27:34,5]<-"na"

counts

AMOP.counts<-subset(counts,counts$sal.sp=="AMOP")
AMOP.counts

AMMA.counts<-subset(counts,counts$sal.sp=="AMMA")

scaph.cut.amma<-subset(scaph.cut,scaph.cut$sal.sp=="AMMA")
temp.list<-c(tapply(scaph.cut.amma$size,list(scaph.cut.amma$treat,scaph.cut.amma$sal.treat),length))
temp.list<-temp.list[!is.na(temp.list)]
AMMA.counts[,5]<-temp.list
rm(temp.list)

AMMA.counts

counts<-rbind(AMOP.counts,AMMA.counts)
counts

# create dataframe of trial results, and swap in new counts

# actually, better to do this mannually in excel

write.csv(AMMA.counts,file="AMMA-counts.csv")
write.csv(AMOP.counts,file="AMOP-counts.csv")

# once I have a .csv with the updated trial results 
# & chesson's alpha/electivity index
# I can try making the linear models
