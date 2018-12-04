### Sensitivity Analysis -------------------------------------------

names(data)

# strictly speaking, if 100% of taxon A is consumed and any of taxa B and C remain,
# preference for taxon A = 1 and preference for taxa B & C = -1
data$ostr.e[which(data$ostr.end == 0)] <- 1
data$cope.e[which(data$ostr.end == 0)] <- -1
data$clado.e[which(data$ostr.end == 0)] <- -1


par(mfrow=c(1,2))
par(mar = c(5, 4, 4, 2), oma = c(0, 0, 0, 0))

# plot overall feeding preferences ----
y=c(mean(cycl.e),mean(scaph.e),mean(harp.e),mean(ost.s.e),mean(simo.e),mean(ost.l.e))
x=c("Cycl","Scaph","Harp","Ostr","Simo","Cypr")
sem <- c(sd(cycl.e)/sqrt(20),sd(scaph.e)/sqrt(20),sd(harp.e)/sqrt(20),
         sd(ost.s.e)/sqrt(20),sd(simo.e)/sqrt(20),sd(ost.l.e)/sqrt(20))
error <- qt(.975,df=19)*c(sd(cycl.e)/sqrt(20),sd(scaph.e)/sqrt(20),sd(harp.e)/sqrt(20),
                          sd(ost.s.e)/sqrt(20),sd(simo.e)/sqrt(20),sd(ost.l.e)/sqrt(20))

plot(y,xaxt="n",ylim=c(min(y-error),max(y+error)),xlab="Prey taxa",ylab="Preference", 
     family="C",main="",pch=16,bty="L")
axis(1,at=c(1,2,3,4,5,6),labels=x,family="C")
for(i in 1:6){segments(i,y[i]-error[i],i,y[i]+error[i]
                       #,col=colors.p[i]
)}
abline(0,0,lty=2)
mtext("(a)", side = 3, line = -16, adj = -.25, cex = 1)

rm(i,sem,x,y)

# now for pooled organisms ----
y=c(mean(cope.e),mean(clado.e),mean(ostr.e))
x=c("Copepoda","Cladocera","Ostracoda")
sem <- c(sd(cope.e)/sqrt(20),sd(clado.e)/sqrt(20),sd(ostr.e)/sqrt(20))
error <- qt(.975,df=19)*sem
plot(y,xaxt="n",ylim=c(min(y-error),max(y+error)),xlab="Prey taxa",ylab="Preference", 
     main="",pch=16,bty="L",family="A")
axis(1,at=c(1,2,3),labels=x,family="A")
for(i in 1:3){segments(i,y[i]-error[i],i,y[i]+error[i]
)}
abline(0,0,lty=2)
mtext("(b)", side = 3, line = -16, adj = -.25, cex = 1)

rm(i,sem,x,y)


# per capita mortality of Cypricercus, by initial frequency
par(mfrow=c(1,1))

seq(from=1,by=4,length.out=5)
avg.mean<-mean(ost.l.eat[c(1,5,9,13,17)]/ost.l.start[c(1,5,9,13,17)])
avg.sem<-sd(ost.l.eat[c(1,5,9,13,17)]/ost.l.start[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(ost.l.eat[c(2,6,10,14,18)]/ost.l.start[c(2,6,10,14,18)])
t1.sem<-sd(ost.l.eat[c(2,6,10,14,18)]/ost.l.start[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(ost.l.eat[c(3,7,11,15,19)]/ost.l.start[c(3,7,11,15,19)])
t2.sem<-sd(ost.l.eat[c(3,7,11,15,19)]/ost.l.start[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(ost.l.eat[c(4,8,12,16,20)]/ost.l.start[c(4,8,12,16,20)])
t3.sem<-sd(ost.l.eat[c(4,8,12,16,20)]/ost.l.start[c(4,8,12,16,20)])/sqrt(5)

y=c(t3.mean,avg.mean,t1.mean,t2.mean)
x=ost.l.start[c(4,1:3)]
sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
plot(y~x,ylim=c(min(y-sem),max(y+sem)),main="",
     xlab="Cypricercus Initial Frequency",ylab="Per Capita Mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}

rm(i,sem,x,y)

# # amt of other zoop eaten by amt of Cypricercus eaten 
# plot((total.eat-ost.l.eat)~ost.l.eat, main="",xlab="Cypricercus consumed",
#      ylab="Other zooplankton consumed",ylim=c(0,450),xlim=c(0,100))
# 
# # amt of other zoop eaten by initial frequency of Cypricercus
# plot((total.eat-ost.l.eat)~ost.l.start, main="",xlab="Initial frequency of Cypricercus",
#      ylab="Other zooplankton consumed",ylim=c(0,450),xlim=c(0,200))
# 
# avg.mean<-mean(total.eat[c(1,5,9,13,17)]-ost.l.eat[c(1,5,9,13,17)])
# avg.sem<-sd(total.eat[c(1,5,9,13,17)]-ost.l.eat[c(1,5,9,13,17)])/sqrt(5)
# t1.mean<-mean(total.eat[c(2,6,10,14,18)]-ost.l.eat[c(2,6,10,14,18)])
# t1.sem<-sd(total.eat[c(2,6,10,14,18)]-ost.l.eat[c(2,6,10,14,18)])/sqrt(5)
# t2.mean<-mean(total.eat[c(3,7,11,15,19)]-ost.l.eat[c(3,7,11,15,19)])
# t2.sem<-sd(total.eat[c(3,7,11,15,19)]-ost.l.eat[c(3,7,11,15,19)])/sqrt(5)
# t3.mean<-mean(total.eat[c(4,8,12,16,20)]-ost.l.eat[c(4,8,12,16,20)])
# t3.sem<-sd(total.eat[c(4,8,12,16,20)]-ost.l.eat[c(4,8,12,16,20)])/sqrt(5)
# 
# y=c(t3.mean,avg.mean,t1.mean,t2.mean)
# x=ost.l.start[c(4,1:3)]
# sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
# plot(y~x,ylim=c(0,450),xlim=c(0,200),main="",
#      xlab="Cypricercus Initial Frequency",ylab="Other zooplankton consumed", 
#      pch=16,bty="L")
# for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}
# 
# rm(i,sem,x,y)
# 
# # cladoceran per capita mortality by initial frequency
# avg.mean<-mean(clado.eat[c(1,5,9,13,17)]/clado.start[c(1,5,9,13,17)])
# avg.sem<-sd(clado.eat[c(1,5,9,13,17)]/clado.start[c(1,5,9,13,17)])/sqrt(5)
# t1.mean<-mean(clado.eat[c(2,6,10,14,18)]/clado.start[c(2,6,10,14,18)])
# t1.sem<-sd(clado.eat[c(2,6,10,14,18)]/clado.start[c(2,6,10,14,18)])/sqrt(5)
# t2.mean<-mean(clado.eat[c(3,7,11,15,19)]/clado.start[c(3,7,11,15,19)])
# t2.sem<-sd(clado.eat[c(3,7,11,15,19)]/clado.start[c(3,7,11,15,19)])/sqrt(5)
# t3.mean<-mean(clado.eat[c(4,8,12,16,20)]/clado.start[c(4,8,12,16,20)])
# t3.sem<-sd(clado.eat[c(4,8,12,16,20)]/clado.start[c(4,8,12,16,20)])/sqrt(5)
# 
# y=c(t3.mean,avg.mean,t1.mean,t2.mean)
# x=clado.start[c(4,1:3)]
# sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
# plot(y~x,ylim=c(min(y-sem),max(y+sem)),main="",
#      xlab="Cladocera Initial Frequency",ylab="Per Capita Mortality", 
#      pch=16,bty="L")
# for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}
# 
# rm(i,sem,x,y)
# 
# # copepod per capita mortality by initial frequency
# avg.mean<-mean(cope.eat[c(1,5,9,13,17)]/cope.start[c(1,5,9,13,17)])
# avg.sem<-sd(cope.eat[c(1,5,9,13,17)]/cope.start[c(1,5,9,13,17)])/sqrt(5)
# t1.mean<-mean(cope.eat[c(2,6,10,14,18)]/cope.start[c(2,6,10,14,18)])
# t1.sem<-sd(cope.eat[c(2,6,10,14,18)]/cope.start[c(2,6,10,14,18)])/sqrt(5)
# t2.mean<-mean(cope.eat[c(3,7,11,15,19)]/cope.start[c(3,7,11,15,19)])
# t2.sem<-sd(cope.eat[c(3,7,11,15,19)]/cope.start[c(3,7,11,15,19)])/sqrt(5)
# t3.mean<-mean(cope.eat[c(4,8,12,16,20)]/cope.start[c(4,8,12,16,20)])
# t3.sem<-sd(cope.eat[c(4,8,12,16,20)]/cope.start[c(4,8,12,16,20)])/sqrt(5)
# 
# y=c(t3.mean,avg.mean,t1.mean,t2.mean)
# x=cope.start[c(4,1:3)]
# sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
# plot(y~x,ylim=c(min(y-sem),max(y+sem)),main="",
#      xlab="Copepoda Initial Frequency",ylab="Per Capita Mortality", 
#      pch=16,bty="L")
# for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}
# 
# rm(i,sem,x,y)
# 
# # ostracod per capita mortality by initial frequency
# avg.mean<-mean(ostr.eat[c(1,5,9,13,17)]/ostr.start[c(1,5,9,13,17)])
# avg.sem<-sd(ostr.eat[c(1,5,9,13,17)]/ostr.start[c(1,5,9,13,17)])/sqrt(5)
# t1.mean<-mean(ostr.eat[c(2,6,10,14,18)]/ostr.start[c(2,6,10,14,18)])
# t1.sem<-sd(ostr.eat[c(2,6,10,14,18)]/ostr.start[c(2,6,10,14,18)])/sqrt(5)
# t2.mean<-mean(ostr.eat[c(3,7,11,15,19)]/ostr.start[c(3,7,11,15,19)])
# t2.sem<-sd(ostr.eat[c(3,7,11,15,19)]/ostr.start[c(3,7,11,15,19)])/sqrt(5)
# t3.mean<-mean(ostr.eat[c(4,8,12,16,20)]/ostr.start[c(4,8,12,16,20)])
# t3.sem<-sd(ostr.eat[c(4,8,12,16,20)]/ostr.start[c(4,8,12,16,20)])/sqrt(5)
# 
# y=c(t3.mean,avg.mean,t1.mean,t2.mean)
# x=ostr.start[c(4,1:3)]
# sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
# plot(y~x,ylim=c(min(y-sem),max(y+sem)),main="",
#      xlab="Ostracoda Initial Frequency",ylab="Per Capita Mortality", 
#      pch=16,bty="L")
# for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}
# 
# rm(i,sem,x,y)



# trial.mlm1<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+weight:clado.start+weight:cope.start+clado.start:cope.start)
# 
# trial.mlm2<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+weight:clado.start+weight:cope.start)
# trial.mlm3<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+weight:cope.start+cope.start:clado.start)
# trial.mlm4<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+weight:clado.start+cope.start:clado.start)
# 
# trial.mlm5<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+weight:clado.start)
# trial.mlm6<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+weight:cope.start)
# trial.mlm7<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+cope.start+clado.start:cope.start)
# 
# trial.mlm8<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+clado.start+weight:clado.start)
# trial.mlm9<-lm(cbind(clado.e,cope.e,ostr.e)~
#                  weight+cope.start+weight:cope.start)
# trial.mlm10<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   clado.start+cope.start+clado.start:cope.start)
# 
# trial.mlm11<-lm(cbind(clado.e,cope.e,ostr.e)~weight+clado.start+cope.start)
# 
# trial.mlm12<-lm(cbind(clado.e,cope.e,ostr.e)~weight+clado.start)
# trial.mlm13<-lm(cbind(clado.e,cope.e,ostr.e)~weight+cope.start)
# trial.mlm14<-lm(cbind(clado.e,cope.e,ostr.e)~clado.start+cope.start)
# 
# trial.mlm15<-lm(cbind(clado.e,cope.e,ostr.e)~clado.start)
# trial.mlm16<-lm(cbind(clado.e,cope.e,ostr.e)~cope.start)
# trial.mlm17<-lm(cbind(clado.e,cope.e,ostr.e)~weight)
# 
# trial.mlm18<-lm(cbind(clado.e,cope.e,ostr.e)~1)
# 
# trial.mlm19<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+clado.start+ostr.start+weight:clado.start+weight:ostr.start)
# trial.mlm20<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+clado.start+ostr.start+weight:clado.start+ostr.start:clado.start)
# trial.mlm21<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+clado.start+ostr.start+weight:ostr.start+ostr.start:clado.start)
# 
# trial.mlm22<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+clado.start+ostr.start+weight:clado.start)
# trial.mlm23<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+clado.start+ostr.start+weight:ostr.start)
# trial.mlm24<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+clado.start+ostr.start+clado.start:ostr.start)
# 
# trial.mlm25<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+ostr.start+weight:ostr.start)
# trial.mlm26<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   clado.start+ostr.start+clado.start:ostr.start)
# 
# trial.mlm27<-lm(cbind(clado.e,cope.e,ostr.e)~weight+clado.start+ostr.start)
# 
# trial.mlm28<-lm(cbind(clado.e,cope.e,ostr.e)~weight+ostr.start)
# trial.mlm29<-lm(cbind(clado.e,cope.e,ostr.e)~clado.start+ostr.start)
# 
# trial.mlm30<-lm(cbind(clado.e,cope.e,ostr.e)~ostr.start)
# 
# trial.mlm31<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+cope.start+ostr.start+weight:cope.start+weight:ostr.start)
# trial.mlm32<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+cope.start+ostr.start+weight:cope.start+cope.start:ostr.start)
# trial.mlm33<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+cope.start+ostr.start+weight:ostr.start+cope.start:ostr.start)
# 
# trial.mlm34<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+cope.start+ostr.start+weight:cope.start)
# trial.mlm35<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+cope.start+ostr.start+weight:ostr.start)
# trial.mlm36<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   weight+cope.start+ostr.start+cope.start:ostr.start)
# 
# trial.mlm37<-lm(cbind(clado.e,cope.e,ostr.e)~
#                   cope.start+ostr.start+cope.start:ostr.start)
# 
# trial.mlm38<-lm(cbind(clado.e,cope.e,ostr.e)~weight+cope.start+ostr.start)
# 
# trial.mlm39<-lm(cbind(clado.e,cope.e,ostr.e)~cope.start+ostr.start)
# 
# models<-list(trial.mlm1,trial.mlm2,trial.mlm3,trial.mlm4,trial.mlm5,
#              trial.mlm6,trial.mlm7,trial.mlm8,trial.mlm9,trial.mlm10,
#              trial.mlm11, trial.mlm12, trial.mlm13, trial.mlm14,
#              trial.mlm15,trial.mlm16,trial.mlm17,trial.mlm18,trial.mlm19,
#              trial.mlm20,trial.mlm21,trial.mlm22,trial.mlm23,trial.mlm24,
#              trial.mlm25,trial.mlm26,trial.mlm27,trial.mlm28,trial.mlm29,
#              trial.mlm30,trial.mlm31,trial.mlm32,trial.mlm33,trial.mlm34,
#              trial.mlm35,trial.mlm36,trial.mlm37,trial.mlm38,trial.mlm39)



# models mlm1 & mlm2 are best supported by AIC
rm(models,trial.mlm1,trial.mlm2c,trial.mlm3,trial.mlm4,trial.mlm5,
   trial.mlm6,trial.mlm7,trial.mlm8,trial.mlm9,trial.mlm10,
   trial.mlm14,trial.mlm15,trial.mlm16,trial.mlm17)



models<-list(trial.mlm1,trial.mlm2c,trial.mlm3,trial.mlm4,
             trial.mlm7,trial.mlm8,trial.mlm9,trial.mlm10,
             trial.mlm14,trial.mlm15,
             trial.mlm16,trial.mlm17,trial.mlm5,trial.mlm6)


trial.mlm1<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+weight:clado.start+weight:cope.start+clado.start:cope.start)

trial.mlm2<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+weight:clado.start+weight:cope.start)
trial.mlm3<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+weight:cope.start+cope.start:clado.start)
trial.mlm4<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+weight:clado.start+cope.start:clado.start)

trial.mlm5<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+weight:clado.start)
trial.mlm6<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+weight:cope.start)
trial.mlm7<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+cope.start+clado.start:cope.start)

trial.mlm8<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+clado.start+weight:clado.start)
trial.mlm9<-lm(cbind(clado.a,cope.a,ostr.a)~
                 weight+cope.start+weight:cope.start)
trial.mlm10<-lm(cbind(clado.a,cope.a,ostr.a)~
                  clado.start+cope.start+clado.start:cope.start)

trial.mlm11<-lm(cbind(clado.a,cope.a,ostr.a)~weight+clado.start+cope.start)

trial.mlm12<-lm(cbind(clado.a,cope.a,ostr.a)~weight+clado.start)
trial.mlm13<-lm(cbind(clado.a,cope.a,ostr.a)~weight+cope.start)
trial.mlm14<-lm(cbind(clado.a,cope.a,ostr.a)~clado.start+cope.start)

trial.mlm15<-lm(cbind(clado.a,cope.a,ostr.a)~clado.start)
trial.mlm16<-lm(cbind(clado.a,cope.a,ostr.a)~cope.start)
trial.mlm17<-lm(cbind(clado.a,cope.a,ostr.a)~weight)

trial.mlm18<-lm(cbind(clado.a,cope.a,ostr.a)~1)

trial.mlm19<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+clado.start+ostr.start+weight:clado.start+weight:ostr.start)
trial.mlm20<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+clado.start+ostr.start+weight:clado.start+ostr.start:clado.start)
trial.mlm21<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+clado.start+ostr.start+weight:ostr.start+ostr.start:clado.start)

trial.mlm22<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+clado.start+ostr.start+weight:clado.start)
trial.mlm23<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+clado.start+ostr.start+weight:ostr.start)
trial.mlm24<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+clado.start+ostr.start+clado.start:ostr.start)

trial.mlm25<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+ostr.start+weight:ostr.start)
trial.mlm26<-lm(cbind(clado.a,cope.a,ostr.a)~
                  clado.start+ostr.start+clado.start:ostr.start)

trial.mlm27<-lm(cbind(clado.a,cope.a,ostr.a)~weight+clado.start+ostr.start)

trial.mlm28<-lm(cbind(clado.a,cope.a,ostr.a)~weight+ostr.start)
trial.mlm29<-lm(cbind(clado.a,cope.a,ostr.a)~clado.start+ostr.start)

trial.mlm30<-lm(cbind(clado.a,cope.a,ostr.a)~ostr.start)

trial.mlm31<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+cope.start+ostr.start+weight:cope.start+weight:ostr.start)
trial.mlm32<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+cope.start+ostr.start+weight:cope.start+cope.start:ostr.start)
trial.mlm33<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+cope.start+ostr.start+weight:ostr.start+cope.start:ostr.start)

trial.mlm34<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+cope.start+ostr.start+weight:cope.start)
trial.mlm35<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+cope.start+ostr.start+weight:ostr.start)
trial.mlm36<-lm(cbind(clado.a,cope.a,ostr.a)~
                  weight+cope.start+ostr.start+cope.start:ostr.start)

trial.mlm37<-lm(cbind(clado.a,cope.a,ostr.a)~
                  cope.start+ostr.start+cope.start:ostr.start)

trial.mlm38<-lm(cbind(clado.a,cope.a,ostr.a)~weight+cope.start+ostr.start)

trial.mlm39<-lm(cbind(clado.a,cope.a,ostr.a)~cope.start+ostr.start)

models<-list(trial.mlm1,trial.mlm2,trial.mlm3,trial.mlm4,trial.mlm5,
             trial.mlm6,trial.mlm7,trial.mlm8,trial.mlm9,trial.mlm10,
             trial.mlm11, trial.mlm12, trial.mlm13, trial.mlm14,
             trial.mlm15,trial.mlm16,trial.mlm17,trial.mlm18,trial.mlm19,
             trial.mlm20,trial.mlm21,trial.mlm22,trial.mlm23,trial.mlm24,
             trial.mlm25,trial.mlm26,trial.mlm27,trial.mlm28,trial.mlm29,
             trial.mlm30,trial.mlm31,trial.mlm32,trial.mlm33,trial.mlm34,
             trial.mlm35,trial.mlm36,trial.mlm37,trial.mlm38,trial.mlm39)

(unlist(lapply(models,extractAIC)))

# now return to investigating the univariate responses, in a more orderly manner ----

# cladocerans ----

clado.lm<-lm(clado.e~weight*clado.start)
summary(clado.lm)
par(mfrow=c(2,2))
avPlots(clado.lm,main="Partial Regressions")
plot(clado.lm)
# residuals show humped trend so add a quadratic term

clado.lm1<-lm(clado.e~weight*clado.start+I(clado.start^2))
clado.lm2<-lm(clado.e~weight+clado.start+I(clado.start^2))
clado.lm3<-lm(clado.e~weight*clado.start)
clado.lm4<-lm(clado.e~clado.start+I(clado.start^2))
clado.lm5<-lm(clado.e~weight+clado.start)
clado.lm6<-lm(clado.e~clado.start)
clado.lm7<-lm(clado.e~weight)
clado.lm8<-lm(clado.e~1)

clado.lm9<-lm(clado.e~weight+log(clado.start))
clado.lm10<-lm(log(clado.e+1)~clado.start*weight)
clado.lm11<-lm(log(clado.e+1)~clado.start+weight)
clado.lm12<-lm(log(clado.e+1)~clado.start)
clado.lm13<-lm(log(clado.e+1)~weight)
clado.lm14<-lm(log(clado.e+1)~1)

clado.models1<-list(clado.lm1,clado.lm2,clado.lm3,clado.lm4,
                    clado.lm5,clado.lm6,clado.lm7,clado.lm8)#,clado.lm9,clado.lm10,clado.lm11,clado.lm12,clado.lm13,clado.lm14)

clado.lm1a<-lm(clado.a~weight*clado.start+I(clado.start^2))
clado.lm2a<-lm(clado.a~weight+clado.start+I(clado.start^2))
clado.lm3a<-lm(clado.a~weight*clado.start)
clado.lm4a<-lm(clado.a~clado.start+I(clado.start^2))
clado.lm5a<-lm(clado.a~weight+clado.start)
clado.lm6a<-lm(clado.a~clado.start)
clado.lm7a<-lm(clado.a~weight)
clado.lm8a<-lm(clado.a~1)
clado.lm9a<-lm(clado.a~weight+log(clado.start))

clado.lm10a<-lm(log(clado.a)~clado.start+weight)
clado.lm11a<-lm(log(clado.a)~clado.start)
clado.lm12a<-lm(log(clado.a)~weight)
clado.lm13a<-lm(log(clado.a)~1)
clado.lm14a<-lm(log(clado.a)~clado.start*weight)

clado.models1a<-list(clado.lm1a,clado.lm2a,clado.lm3a,clado.lm4a,
                     clado.lm5a,clado.lm6a,clado.lm7a,clado.lm8a,clado.lm9a,
                     clado.lm10a,clado.lm11a,clado.lm12a,clado.lm13a,clado.lm14a)


# per Burnham & Anderson, AICc is recommended when n/K<40
# here, n=20, K is from 1 to 4

library(AICcmodavg)

(aicc<-unlist(lapply(clado.models1a,AICc)))
which(aicc<=(min(aicc)+4))
# clado.models1[which(aicc<=(min(aicc)+4))]

(aicc<-unlist(lapply(clado.models1,AICc)))
which(aicc<=(min(aicc)+2))

# clado.lm1 is w/in 2 units of minimum using AIC & 3.82 units using AICc
summary(clado.lm1)
plot(clado.lm1)
avPlots(clado.lm1)
# but residuals are still hump-shaped

# clado.lm2 is minimum using both methods
summary(clado.lm2)
plot(clado.lm2)
avPlots(clado.lm2)
# but residuals are still hump-shaped

# clado.lm4 is w/in 2 units of minimum using AICc & 2.44 units using AIC
summary(clado.lm4)
plot(clado.lm4)
avPlots(clado.lm4)
# residuals do not have any clear pattern (subjective judgement)

AIC(clado.lm4)-AIC(clado.lm2)
AICc(clado.lm4)-AICc(clado.lm2)

AICc(clado.lm1)
AICc(clado.lm2)
AICc(clado.lm4)
#clado.lm2
summary(clado.lm2)
summary(clado.lm4)
summary(clado.lm1)

summary(clado.lm2a)
plot(clado.lm2a)
avPlots(clado.lm2a)

summary(clado.lm4a)
plot(clado.lm4a)
avPlots(clado.lm4a)


rm(clado.lm,clado.lm3,clado.lm4p,clado.lm5,clado.lm6,
   clado.lm7,clado.lm8,clado.lm9,clado.models1,clado.models1a,aic,aicc)
rm(clado.lm10,clado.lm11,clado.lm12,clado.lm13,clado.lm14)
rm(clado.lm1a,clado.lm3a,clado.lm5a,clado.lm6a,clado.lm7a,
   clado.lm8a,clado.lm9a,clado.lm10a,clado.lm11a,clado.lm12a,clado.lm13a,clado.lm14a)



# copepods ----

cope.lm1<-lm(cope.e~weight*cope.start+I(cope.start^2))
cope.lm2<-lm(cope.e~weight+cope.start+I(cope.start^2))
cope.lm3<-lm(cope.e~cope.start+I(cope.start^2))
cope.lm4<-lm(cope.e~weight*cope.start)
cope.lm5<-lm(cope.e~weight+cope.start)
cope.lm6<-lm(cope.e~cope.start)
cope.lm7<-lm(cope.e~weight)
cope.lm8<-lm(cope.e~1)
cope.lm9<-lm(cope.e~weight+log(cope.start))
cope.lm10<-lm(cope.e~weight*log(cope.start)+I((log(cope.start))^2))
cope.lm11<-lm(cope.e~weight+log(cope.start)+I((log(cope.start))^2))
cope.lm12<-lm(cope.e~log(cope.start)+I((log(cope.start))^2))
cope.lm13<-lm(log(cope.e+1)~weight*cope.start)
cope.lm14<-lm(log(cope.e+1)~weight+cope.start)
cope.lm15<-lm(log(cope.e+1)~cope.start)
cope.lm16<-lm(log(cope.e+1)~weight)
cope.lm17<-lm(log(cope.e+1)~1)

cope.lm1a<-lm(cope.a~weight*cope.start+I(cope.start^2))
cope.lm2a<-lm(cope.a~weight+cope.start+I(cope.start^2))
cope.lm3a<-lm(cope.a~cope.start+I(cope.start^2))
cope.lm4a<-lm(cope.a~weight*cope.start)
cope.lm5a<-lm(cope.a~weight+cope.start)
cope.lm6a<-lm(cope.a~cope.start)
cope.lm7a<-lm(cope.a~weight)
cope.lm8a<-lm(cope.a~1)
cope.lm9a<-lm(log(cope.a)~weight*cope.start)
cope.lm10a<-lm(log(cope.a)~weight+cope.start)
cope.lm11a<-lm(log(cope.a)~cope.start)
cope.lm12a<-lm(log(cope.a)~weight)
cope.lm13a<-lm(log(cope.a)~1)

cope.models<-list(cope.lm1,cope.lm2,cope.lm3,cope.lm4,cope.lm5,cope.lm6,
                  cope.lm7,cope.lm8)#,cope.lm9,cope.lm10,cope.lm11,cope.lm12,cope.lm13,cope.lm14,cope.lm15,cope.lm16,cope.lm17)

cope.models.a<-list(cope.lm1a,cope.lm2a,cope.lm3a,cope.lm4a,cope.lm5a,cope.lm6a,
                    cope.lm7a,cope.lm8a,cope.lm9a,cope.lm10a,cope.lm11a,cope.lm12a,cope.lm13a)


(aicc<-unlist(lapply(cope.models.a,AICc)))
which(aicc<=(min(aicc)+4))
# cope.models[which(aicc<=(min(aicc)+2))]
(aicc<-unlist(lapply(cope.models,AICc)))
which(aicc<=(min(aicc)+2))

summary(cope.lm3a)

# aic 3.9514 (+1.448), aicc 10.4130 (+3.624)
summary(cope.lm1)
plot(cope.lm1)
# residuals humped?
avPlots(cope.lm1)

# aic 2.5030, aicc 6.7888
summary(cope.lm2)
plot(cope.lm2)
# weird pattern in resids
avPlots(cope.lm2)

# aic 5.7402 (+3.237), aicc 8.4068 (+1.618)
summary(cope.lm3)
plot(cope.lm3)
# resids pattern unclear
avPlots(cope.lm3)

plot(cope.lm11,main="cope.lm11")
plot(cope.lm2,main="cope.lm2")
plot(cope.lm3,main="cope.lm3")
avPlots(cope.lm11)

summary(cope.lm5)
plot(cope.lm5)
avPlots(cope.lm5)

summary(cope.lm7)
plot(cope.lm7)
avPlots(cope.lm7)

summary(cope.lm8)
plot(cope.lm8)
avPlots(cope.lm8)


AICc(cope.lm2)
AICc(cope.lm5)
AICc(cope.lm7)

rm(aicc,cope.models,cope.lm1,cope.lm3,cope.lm4,
   cope.lm6,cope.lm8,cope.lm9,cope.lm10,cope.lm11,cope.lm12)
rm(cope.lm13,cope.lm14,cope.lm15,cope.lm16,cope.lm17)

rm(cope.lm1a,cope.lm2a,cope.lm3a,cope.lm4a,cope.lm5a,cope.lm6a,cope.lm7a,cope.lm8a,
   cope.lm9a,cope.lm10a,cope.lm11a,cope.lm12a,cope.lm13a,cope.models.a)

# ostracods ----

ostr.lm1<-lm(ostr.e~weight*ostr.start+I(ostr.start^2))
ostr.lm2<-lm(ostr.e~weight+ostr.start+I(ostr.start^2))
ostr.lm3<-lm(ostr.e~ostr.start+I(ostr.start^2))
ostr.lm4<-lm(ostr.e~weight*ostr.start)
ostr.lm5<-lm(ostr.e~weight+ostr.start)
ostr.lm6<-lm(ostr.e~ostr.start)
ostr.lm7<-lm(ostr.e~weight)
ostr.lm8<-lm(ostr.e~1)
ostr.lm9<-lm(log(ostr.e+1)~weight*ostr.start)
ostr.lm10<-lm(log(ostr.e+1)~weight+ostr.start)
ostr.lm11<-lm(log(ostr.e+1)~ostr.start)
ostr.lm12<-lm(log(ostr.e+1)~weight)
ostr.lm13<-lm(log(ostr.e+1)~1)

ostr.models<-list(ostr.lm1,ostr.lm2,ostr.lm3,ostr.lm4,ostr.lm5,ostr.lm6,
                  ostr.lm7,ostr.lm8)#,ostr.lm9,ostr.lm10,ostr.lm11,ostr.lm12,ostr.lm13)

(aicc<-unlist(lapply(ostr.models,AICc)))
which(aicc<=(min(aicc)+4))
# ostr.models[which(aicc<=(min(aicc)+2))]
# 1,2,4,5; 9,10


ostr.lm1a<-lm(ostr.a~weight*ostr.start+I(ostr.start^2))
ostr.lm2a<-lm(ostr.a~weight+ostr.start+I(ostr.start^2))
ostr.lm3a<-lm(ostr.a~ostr.start+I(ostr.start^2))
ostr.lm4a<-lm(ostr.a~weight*ostr.start)
ostr.lm5a<-lm(ostr.a~weight+ostr.start)
ostr.lm6a<-lm(ostr.a~ostr.start)
ostr.lm7a<-lm(ostr.a~weight)
ostr.lm8a<-lm(ostr.a~1)
ostr.lm9a<-lm(log(ostr.a)~weight*ostr.start)
ostr.lm10a<-lm(log(ostr.a)~weight+ostr.start)
ostr.lm11a<-lm(log(ostr.a)~ostr.start)
ostr.lm12a<-lm(log(ostr.a)~weight)
ostr.lm13a<-lm(log(ostr.a)~1)

ostr.models.a<-list(ostr.lm1a,ostr.lm2a,ostr.lm3a,ostr.lm4a,ostr.lm5a,ostr.lm6a,
                    ostr.lm7a,ostr.lm8a,ostr.lm9a,ostr.lm10a,ostr.lm11a,ostr.lm12a,ostr.lm13a)

(aicc<-unlist(lapply(ostr.models.a,AICc)))
which(aicc<=(min(aicc)+4))


AICc(ostr.lm1)
AICc(ostr.lm2)
AICc(ostr.lm5)
#lm1
AICc(ostr.lm9)
AICc(ostr.lm10)
#lm10,lm9

AICc(ostr.lm1a)
AICc(ostr.lm2a)
AICc(ostr.lm4a)
AICc(ostr.lm5a)
# lm2a & lm5a

rm(aicc,ostr.models,ostr.lm3,ostr.lm4,
   ostr.lm6,ostr.lm7,ostr.lm8,ostr.lm11,ostr.lm12,ostr.lm13)
rm(ostr.models.a,ostr.lm3a,ostr.lm6a,ostr.lm7a,ostr.lm8a,
   ostr.lm9a,ostr.lm10a,ostr.lm11a,ostr.lm12a,ostr.lm13a)

summary(ostr.lm1)
plot(ostr.lm1)
avPlots(ostr.lm1)


summary(ostr.lm2)
plot(ostr.lm2)
avPlots(ostr.lm2)

summary(ostr.lm4)
plot(ostr.lm4)
avPlots(ostr.lm4)

summary(ostr.lm5)
plot(ostr.lm5)
avPlots(ostr.lm5)


summary(ostr.lm9)
plot(ostr.lm9)
avPlots(ostr.lm9)
# resids strongly hump-shaped

summary(ostr.lm10)
plot(ostr.lm10)
avPlots(ostr.lm10)
# resids strongly hump-shaped

summary(ostr.lm1a)
plot(ostr.lm1a)
avPlots(ostr.lm1a)


summary(ostr.lm2a)
plot(ostr.lm2a)
avPlots(ostr.lm2a)


summary(ostr.lm4a)
plot(ostr.lm4a)
avPlots(ostr.lm4a)


summary(ostr.lm5a)
plot(ostr.lm5a)
avPlots(ostr.lm5a)


# 3D preference figures ----

# actual ranges of start values
sal.range<-seq(from=0.55,to=1.05,by=0.1)
clado.range<-seq(from=10,to=190,by=10)
cope.range<-seq(from=40,to=320,by=20)
ostr.range<-seq(from=10,to=230,by=10)

# generate predictions:

# cladocerans ----
# create grid called clado.pred.list
clado.pred.list<-as.data.frame(expand.grid(list(sal.range,clado.range)))
colnames(clado.pred.list)<-c("weight","clado.start")

# populate predictions for preference across this grid
clado.predict<-predict(clado.lm,clado.pred.list,se.fit=TRUE,type="response")

# convert all this business into a single data frame
clado.predict<-data.frame(clado.pred.list,clado.predict)
rm(clado.pred.list)

# copepods ----
cope.pred.list<-as.data.frame(expand.grid(list(sal.range,cope.range)))
colnames(cope.pred.list)<-c("weight","cope.start")
cope.predict<-predict(cope.lm,cope.pred.list,se.fit=TRUE,type="response")

cope.predict<-data.frame(cope.pred.list,cope.predict)
rm(cope.pred.list)

# ostracods ----
ostr.pred.list<-as.data.frame(expand.grid(list(sal.range,ostr.range)))
colnames(ostr.pred.list)<-c("weight","ostr.start")
ostr.predict<-predict(ostr.lm,ostr.pred.list,se.fit=TRUE,type="response")

ostr.predict<-data.frame(ostr.pred.list,ostr.predict)
rm(ostr.pred.list)

# cleanup ----
rm(sal.range,clado.lm,clado.x.pred,clado.y.pred,clado.xy,clado.predict,clado.range)



## sensitivity of R3-T3 and R4-T3 to substitution of 1 leftover cypricercus ----------------------------------
data
r3<-data %>%
  filter(round == 3)
r3 <- data[1:9,]
r4 <- data[10:18,]
plot(r3$clado.e~r3$Substitution,ylim=c(-1,1),col="blue")
points(r3$cope.e~r3$Substitution,col="green")
points(r3$ostr.e~r3$Substitution,col="red")

ggplot(r3, aes(x = Substitution, y = clado.e)) +
  coord_cartesian(ylim = c(-1,1)) +
  geom_point(color="blue")+
  geom_line(color="blue")+
  geom_point(data=r3,aes(x=Substitution,y=cope.e),color="green")+
  geom_line(data=r3,aes(x=Substitution,y=cope.e),color="green")+
  geom_point(data=r3,aes(x=Substitution,y=ostr.e),color="red")+
  geom_line(data=r3,aes(x=Substitution,y=ostr.e),color="red")+
  labs(x="Substitution",y="Preference")+
  ggtitle("R3 - T3")

ggplot(r4, aes(x = Substitution, y = clado.e)) +
  coord_cartesian(ylim = c(-1,1)) +
  geom_point(color="blue")+
  geom_line(color="blue")+
  geom_point(data=r4,aes(x=Substitution,y=cope.e),color="green")+
  geom_line(data=r4,aes(x=Substitution,y=cope.e),color="green")+
  geom_point(data=r4,aes(x=Substitution,y=ostr.e),color="red")+
  geom_line(data=r4,aes(x=Substitution,y=ostr.e),color="red")+
  labs(x="Substitution",y="Preference")+
  ggtitle("R4 - T3")

