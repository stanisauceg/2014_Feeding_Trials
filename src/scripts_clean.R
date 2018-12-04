### Load & prep data ############################# 

library(car)
library(bbmle)
library(lubridate)
library(AICcmodavg)
library(optimx)
library(tidyverse)
library(grid)
library(gridExtra)

path <- file.path("data", "2014_feeding_trials_nearly_raw_data.csv")
data.full<-read.csv(path)

head(data.full)
names(data.full)

# extract the raw data for analysis
data.raw <- data.full %>%
  select(-bucket, -notes, -'X0.5x.lens')
rm(data.full)

names(data.raw)
print(data.raw, n = "all")
min(data.raw$cycl_start - data.raw$cycl_end)
# one observation has 12 cyclopoid copepods remaining when there were 11 at start
# assume that no copepods were eaten, since negative consumption is not possible
# therefore leftover should be equal to start value

# identify this observation
templist <- data.raw$cycl_start - data.raw$cycl_end
index <- which(templist == min(templist))

# make the substitution
data.raw$cycl_end[index] <- data.raw$cycl_start[index]
rm(templist, index)

# retain adjusted scapholeberis leftover counts, 
# rename columns for ease of use
data.raw <- data.raw %>%
  select(-scaph_end) %>%
  rename(simo.start = simo_start,
         scaph.start = scaph_start,
         harp.start = harp_start,
         cycl.start = cycl_start,
         ost.l.start = ost.l_start,
         ost.s.start = ost.s_start,
         simo.end = simo_end,
         scaph.end = scaph_end_adjusted,
         harp.end = harp_end,
         cycl.end = cycl_end,
         ost.l.end = ost.l_end,
         ost.s.end = ost.s_end)

# pool taxa into three main types: clado, cope, ostr
data.pooled <- data.raw %>%
  mutate(clado.start = simo.start + scaph.start,
         cope.start = cycl.start + harp.start,
         ostr.start = ost.l.start + ost.s.start,
         clado.end = simo.end + scaph.end,
         cope.end = cycl.end + harp.end,
         ostr.end = ost.l.end + ost.s.end)
rm(data.raw)

# calculate how many zoop of each type were eaten
data.eaten <- data.pooled %>%
  filter(treatment != "ctrl") %>%
  mutate(simo.eat = simo.start - simo.end,
         scaph.eat = scaph.start - scaph.end,
         harp.eat = harp.start - harp.end,
         cycl.eat = cycl.start - cycl.end,
         ost.l.eat = ost.l.start - ost.l.end,
         ost.s.eat = ost.s.start - ost.s.end,
         clado.eat = clado.start - clado.end,
         cope.eat = cope.start - cope.end,
         ostr.eat = ostr.start - ostr.end,
         total.eat = clado.eat + cope.eat + ostr.eat)
rm(data.pooled)

# note, 0 ostracods remain even after pooling in T3 for rounds 3 & 4

# define function to calculate log of fraction remaining
# if all individuals were eaten, count one leftover to avoid log(0)
log.leftover <- function(start, eat) {
  while(eat == start){
    eat <- eat - 1
  }
  log((start - eat)/start)
}

# This calculation avoids disturbing the actual leftover counts,
# so the functional response analyses later on are unaffected by the substitution.
# We will need to evaluate the effect of the substition on preferences separately.

# calculate log of fraction remaining, to use for Chesson's alpha calculation in next step
data.ln <- data.eaten %>%
  mutate(simo.ln = map2_dbl(simo.start, simo.eat, log.leftover),
         scaph.ln = map2_dbl(scaph.start, scaph.eat, log.leftover),
         harp.ln = map2_dbl(harp.start, harp.eat, log.leftover),
         cycl.ln = map2_dbl(cycl.start, cycl.eat, log.leftover),
         ost.l.ln = map2_dbl(ost.l.start, ost.l.eat, log.leftover),
         ost.s.ln = map2_dbl(ost.s.start, ost.s.eat, log.leftover),
         clado.ln = map2_dbl(clado.start, clado.eat, log.leftover),
         cope.ln = map2_dbl(cope.start, cope.eat, log.leftover),
         ostr.ln = map2_dbl(ostr.start, ostr.eat, log.leftover))
rm(data.eaten)

# calculate Chesson's alpha
data.a <- data.ln %>%
  mutate(sum.log.leftover = simo.ln +scaph.ln + harp.ln + cycl.ln + ost.l.ln + ost.s.ln,
         sum.log.pooled.leftover = clado.ln + cope.ln + ostr.ln,
         simo.a = simo.ln/sum.log.leftover,
         scaph.a = scaph.ln/sum.log.leftover,
         harp.a = harp.ln/sum.log.leftover,
         cycl.a = cycl.ln/sum.log.leftover,
         ost.l.a = ost.l.ln/sum.log.leftover,
         ost.s.a = ost.s.ln/sum.log.leftover,
         clado.a = clado.ln/sum.log.pooled.leftover,
         cope.a = cope.ln/sum.log.pooled.leftover,
         ostr.a = ostr.ln/sum.log.pooled.leftover)
rm(data.ln)

# define function to calculate electivity index epsilon from Chesson's alpha
electivity <- function(taxon.alpha, n.choices) {
  (n.choices * taxon.alpha - 1)/((n.choices - 2) * taxon.alpha + 1)
}

# calculate electivity index epsilon from Chesson's alpha
# for individual taxa, there are 6 choices; for pooled taxa, 3 choices (clado/cope/ostr)
data <- data.a %>%
  mutate(simo.e = electivity(simo.a, 6),
         scaph.e = electivity(scaph.a, 6),
         harp.e = electivity(harp.a, 6),
         cycl.e = electivity(cycl.a, 6),
         ost.l.e = electivity(ost.l.a, 6),
         ost.s.e = electivity(ost.s.a, 6),
         clado.e = electivity(clado.a, 3),
         cope.e = electivity(cope.a, 3),
         ostr.e = electivity(ostr.a, 3))
rm(data.a)

# clean up, prep for further analysis
str(data)

# convert Round to factor
data$round<-as.factor(data$round)

# convert trial.start.date to date format
data$trial.start.date <- mdy(data$trial.start.date)

# update levels of treatment & block factors
data$treatment <- factor(data$treatment)
data$block <- factor(data$block)

# remove intermediate log.leftover columns
data <- data %>%
  select(-ends_with(".ln"), -ends_with("leftover"))

str(data)
# ready to go!


####  Preference Analysis ##########################################################################

## Multi-taxon preferences: multivariate linear modeling ----------------------------------------

attach(data)

# set response vbls to test either Chesson's epsilon or alpha:
response <- cbind(clado.e, cope.e, ostr.e)
# response <- cbind(clado.a, cope.a, ostr.a)

# fit a selection of models w/ frequency dependence
mlm1 <- lm(response ~ clado.start + cope.start + weight + 
             clado.start:cope.start + clado.start:weight + cope.start:weight)
mlm2 <- lm(response ~ clado.start + cope.start + weight + 
             clado.start:weight + cope.start:weight)
mlm3 <- lm(response ~ clado.start + ostr.start + weight + 
             clado.start:ostr.start + ostr.start:weight)
mlm4 <- lm(response ~ clado.start + cope.start + weight)
mlm5 <- lm(response ~ clado.start + cope.start + clado.start:cope.start)
mlm6 <- lm(response ~ clado.start + cope.start)
mlm7 <- lm(response ~ cope.start + weight + cope.start:weight)
mlm8 <- lm(response ~ ostr.start + weight + ostr.start:weight)
mlm9 <- lm(response ~ clado.start + weight + clado.start:weight)
mlm10 <- lm(response ~ cope.start + weight)
mlm11 <- lm(response ~ ostr.start + weight)
mlm12 <- lm(response ~ clado.start + weight)
mlm13 <- lm(response ~ ostr.start)
mlm14 <- lm(response ~ cope.start)
mlm15 <- lm(response ~ clado.start)

# fit models w/o frequency dependence
mlm16 <- lm(response ~ weight)
mlm17 <- lm(response ~ 1)

# n.params: n.terms in each model counted by hand, then +2 for intercept & variance
n.params <- 2 + c(6,5,5,3,3,2,3,3,3,2,2,2,1,1,1,1,0) 

# put all models in a list
models <- list(mlm1, mlm2, mlm3, mlm4, mlm5, mlm6, mlm7, mlm8, mlm9, mlm10, 
               mlm11, mlm12, mlm13, mlm14, mlm15, mlm16, mlm17)

# get AIC (uncorrected), put into a data frame
AIC.vals <- (unlist(lapply(models,extractAIC)))[2*seq_along(models)]
IDs <- seq_along(models)
model_stats <- data.frame(cbind(IDs, AIC.vals, n.params))

# extract adj.r.squared from summary(mlm)

# str(summary(mlm1))
# summary(mlm1)[[1]][[9]]
# map(models, summary)[[1]][[1]][[9]]
# map(models, summary)[[1]][2]

for(i in seq_along(models)) {
  model_stats$clado.adj.r2[i] <- map(models, summary)[[i]][[1]][[9]]
  model_stats$cope.adj.r2[i] <- map(models, summary)[[i]][[2]][[9]]
  model_stats$ostr.adj.r2[i] <- map(models, summary)[[i]][[3]][[9]]
}


# calculate AICc
## AICc = AIC + 2k(k + 1)/(n - k - 1)
## we have n = 20 observations

model_stats <- model_stats %>%
  rename(AIC = AIC.vals) %>%
  mutate(AICc = AIC + 2*n.params*(n.params + 1)/(20 - n.params - 1)) %>%
  arrange(AICc)

model_stats %>% arrange(IDs)

rm(AIC.vals, IDs, n.params, models, response, i,
  mlm1, mlm2, mlm3, mlm4, mlm5, mlm6, mlm7, mlm8, mlm9, mlm10, 
               mlm11, mlm12, mlm13, mlm14, mlm15, mlm16, mlm17)

### Univariate (single taxon) preferences --------------------------------

### Cladocerans

clado.lm<-lm(clado.e~weight*clado.start)
summary(clado.lm)
par(mfrow=c(2,2))
avPlots(clado.lm,main="Partial Regressions")
plot(clado.lm)
# residuals show humped trend so add a quadratic term
rm(clado.lm)


# full set of models
clado.lm1 <- lm(clado.e ~ weight * clado.start + I(clado.start^2))
clado.lm2 <- lm(clado.e ~ weight + clado.start + I(clado.start^2))
clado.lm3 <- lm(clado.e ~ clado.start + I(clado.start^2))
clado.lm4 <- lm(clado.e ~ weight * clado.start)
clado.lm5 <- lm(clado.e ~ weight + clado.start)
clado.lm6 <- lm(clado.e ~ clado.start)
clado.lm7 <- lm(clado.e ~ weight)
clado.lm8 <- lm(clado.e ~ 1)

clado.models <- list(clado.lm1, clado.lm2, clado.lm3, clado.lm4,
                    clado.lm5, clado.lm6, clado.lm7, clado.lm8)

# compare AICc scores

clado.aicc <- clado.models %>% lapply(AICc) %>% unlist()
clado.aicc
which(clado.aicc <= (min(clado.aicc) + 2))
clado.models[which(clado.aicc <= (min(clado.aicc) + 2))]
# note, AICc supports models both w/ and w/o salamander weight; marginally lower AICc for model w/ weight
AICc(clado.lm3)-AICc(clado.lm2)

# clado.lm2 has minimum AICc
summary(clado.lm2)
plot(clado.lm2)
avPlots(clado.lm2)

# clado.lm3 is w/in 2 units of minimum using AICc
summary(clado.lm3)
plot(clado.lm3)
avPlots(clado.lm3)


# cleanup
rm(clado.lm1,clado.lm2,clado.lm3,clado.lm4,clado.lm5,clado.lm6,clado.lm7,clado.lm8, 
   clado.models, clado.aicc)

### Copepods

# full set of models
cope.lm1 <- lm(cope.e ~ weight * cope.start + I(cope.start^2))
cope.lm2 <- lm(cope.e ~ weight + cope.start + I(cope.start^2))
cope.lm3 <- lm(cope.e ~ cope.start + I(cope.start^2))
cope.lm4 <- lm(cope.e ~ weight * cope.start)
cope.lm5 <- lm(cope.e ~ weight + cope.start)
cope.lm6 <- lm(cope.e ~ cope.start)
cope.lm7 <- lm(cope.e ~ weight)
cope.lm8 <- lm(cope.e ~ 1)

cope.models <- list(cope.lm1, cope.lm2, cope.lm3, cope.lm4, cope.lm5, cope.lm6,
                  cope.lm7, cope.lm8)

cope.aicc <- cope.models %>% lapply(AICc) %>% unlist()
cope.aicc
which(cope.aicc <= (min(cope.aicc)+2))
cope.models[which(cope.aicc <= (min(cope.aicc) + 2))]
# once again, weight is in best model but the model w/o weight is also supported

# aic 2.5030, aicc 6.7888
summary(cope.lm2)
plot(cope.lm2)
avPlots(cope.lm2)

# aic 5.7402 (+3.237), aicc 8.4068 (+1.618)
summary(cope.lm3)
plot(cope.lm3)
avPlots(cope.lm3)

rm(cope.lm1,cope.lm2,cope.lm3,cope.lm4,cope.lm5,
   cope.lm6,cope.lm7, cope.lm8, cope.models, cope.aicc)

## Ostracods
ostr.lm1<-lm(ostr.e~weight*ostr.start+I(ostr.start^2))
ostr.lm2<-lm(ostr.e~weight+ostr.start+I(ostr.start^2))
ostr.lm3<-lm(ostr.e~ostr.start+I(ostr.start^2))
ostr.lm4<-lm(ostr.e~weight*ostr.start)
ostr.lm5<-lm(ostr.e~weight+ostr.start)
ostr.lm6<-lm(ostr.e~ostr.start)
ostr.lm7<-lm(ostr.e~weight)
ostr.lm8 <- lm(ostr.e ~ 1)

ostr.models<-list(ostr.lm1,ostr.lm2,ostr.lm3,ostr.lm4,ostr.lm5,ostr.lm6,
                  ostr.lm7, ostr.lm8)

(aicc<-unlist(lapply(ostr.models,AICc)))
which(aicc<=(min(aicc)+2))
ostr.models[which(aicc<=(min(aicc)+2))]
# most complex model is best

summary(ostr.lm1)
plot(ostr.lm1)
avPlots(ostr.lm1)

rm(aicc,ostr.models,ostr.lm1,ostr.lm2,ostr.lm3,ostr.lm4,ostr.lm5,
   ostr.lm6,ostr.lm7)


# best single-taxon models
clado.lm <- lm(clado.e ~ weight + clado.start + I(clado.start^2), data = data)
cope.lm <- lm(cope.e ~ weight + cope.start + I(cope.start^2), data = data)
ostr.lm <- lm(ostr.e ~ weight * ostr.start + I(ostr.start^2), data = data)

# 3D scatterplots of preference w/ prediction surfaces ----
library(plot3D)

# Cladocerans: weight + start + start^2 ----
x<-clado.start
y<-weight
z<-clado.e

fit<-lm(z~x+y+I(x^2))

grid.lines = 26

x.pred<-seq(min(x),max(x),length.out=grid.lines)
y.pred<-seq(min(y),max(y),length.out=grid.lines)
xy<-expand.grid(x=x.pred,y=y.pred)

z.pred<-matrix(predict(fit,newdata=xy),nrow=grid.lines,ncol=grid.lines)

# fitted points for droplines to surface
fitpoints <- predict(fit)

# scatter plot with regression plane
par(mfrow=c(1,1))
plot3D::scatter3D(clado.start, weight, clado.e, pch = 18, cex = 1, 
                  theta = 30, phi = 00, ticktype = "detailed",
                  xlab = "Initial frequency", ylab = "", zlab = "Preference",
                  surf=list(x=x.pred,y=y.pred,z=z.pred,
                            facets = NA, fit = fitpoints), main = "Preference for cladocerans")
dims<-par("usr")
x<-dims[1]+0.85*diff(dims[1:2])
y<-dims[3]+0.1*diff(dims[3:4])
text(x,y,"Salamander weight (g)",srt=52)

# cleanup
rm(x,y,z,xy,x.pred,y.pred,z.pred,fit,fitpoints,grid.lines,dims)


## 2-D preference plots --------------------------

# select quantiles of salamander weights to display fitted curves
weights.lo.hi <- quantile(weight, probs = c(0.1, 0.5, 0.9))

# set size of points across all figures
pt.size <- 3
jitter <- position_jitter(width = 1)

# Cladocerans ----

clado.coeffs <- clado.lm$coefficients 
clado.q1.curve <- function(x) clado.coeffs[[1]] + clado.coeffs[[2]]*weights.lo.hi[[1]] + clado.coeffs[[3]]*x + clado.coeffs[[4]]*x^2
clado.q5.curve <- function(x) clado.coeffs[[1]] + clado.coeffs[[2]]*weights.lo.hi[[2]] + clado.coeffs[[3]]*x + clado.coeffs[[4]]*x^2
clado.q9.curve <- function(x) clado.coeffs[[1]] + clado.coeffs[[2]]*weights.lo.hi[[3]] + clado.coeffs[[3]]*x + clado.coeffs[[4]]*x^2

newx <- seq(min(clado.start), max(clado.start), by=0.1)
clado.q1.pred.interval <- predict(clado.lm, newdata=data.frame(clado.start=newx, weight = weights.lo.hi[[1]]), 
                                 interval="confidence", level = 0.95)
clado.q1.pred.interval <- as.data.frame(clado.q1.pred.interval)
clado.q1.pred.interval$clado.start <- newx
clado.q5.pred.interval <- predict(clado.lm, newdata=data.frame(clado.start=newx, weight = weights.lo.hi[[2]]), 
                                 interval="confidence", level = 0.95)
clado.q5.pred.interval <- as.data.frame(clado.q5.pred.interval)
clado.q5.pred.interval$clado.start <- newx
clado.q9.pred.interval <- predict(clado.lm, newdata=data.frame(clado.start=newx, weight = weights.lo.hi[[3]]), 
                                 interval="confidence", level = 0.95)
clado.q9.pred.interval <- as.data.frame(clado.q9.pred.interval)
clado.q9.pred.interval$clado.start <- newx
rm(newx)

clado.q1.pred.interval$lwr[which(clado.q1.pred.interval$lwr < -1)] <- -1
clado.q5.pred.interval$lwr[which(clado.q5.pred.interval$lwr < -1)] <- -1
clado.q9.pred.interval$lwr[which(clado.q9.pred.interval$lwr < -1)] <- -1


clado.1 <- ggplot(data, aes(x = clado.start, y = clado.e, fill = weight)) +
  # stat_smooth(method="lm",se=FALSE,formula=y~x+I(x^2),color="red") +
  # geom_ribbon(inherit.aes = FALSE, data = clado.q1.pred.interval, aes(x = clado.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "dark goldenrod 2") +
  # geom_ribbon(inherit.aes = FALSE, data = clado.q5.pred.interval, aes(x = clado.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "gray 40") +
  # geom_ribbon(inherit.aes = FALSE, data = clado.q9.pred.interval, aes(x = clado.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "dark blue") +
  stat_function(fun = clado.q1.curve, color = "dark goldenrod 2", size = 1, lty = 2, alpha = .6) +
  stat_function(fun = clado.q5.curve, color = "gray 40", size = 1, lty = 2, alpha = .6) +
  stat_function(fun = clado.q9.curve, color = "dark blue", size = 1, lty = 2, alpha = .6) +
  geom_point(pch = 21, size = pt.size, position = jitter, alpha = 0.9) +
  labs(x="Cladoceran initial frequency",y="Preference") +
  scale_fill_viridis_c("Salamander \nweight (g)", option = "E", direction = -1) +
  theme_classic()

clado.2 <- ggplot(data, aes(x = weight, y = clado.e, color = factor(clado.start))) +
  geom_point(size = pt.size)+
  stat_smooth(aes(group = "identity"), color = "red", method="lm",se=FALSE)+
  labs(x="Salamander weight (g)",y="Preference")+
  scale_color_viridis_d("Cladoceran \ninitial \nfrequency", option = "D", direction = 1) +
  theme_classic()

clado.3 <- ggplot(data, aes(x = clado.start, color = clado.e, y = weight)) +
  geom_point(size = pt.size)+
  labs(x="Cladoceran initial frequency",y="Salamander weight (g)")+
  scale_color_viridis_c("Preference", option = "B", direction = 1) +
  theme_classic()


# Copepods ----

cope.coeffs <- cope.lm$coefficients 
cope.q1.curve <- function(x) cope.coeffs[[1]] + cope.coeffs[[2]]*weights.lo.hi[[1]] + cope.coeffs[[3]]*x + cope.coeffs[[4]]*x^2
cope.q5.curve <- function(x) cope.coeffs[[1]] + cope.coeffs[[2]]*weights.lo.hi[[2]] + cope.coeffs[[3]]*x + cope.coeffs[[4]]*x^2
cope.q9.curve <- function(x) cope.coeffs[[1]] + cope.coeffs[[2]]*weights.lo.hi[[3]] + cope.coeffs[[3]]*x + cope.coeffs[[4]]*x^2


newx <- seq(min(cope.start), max(cope.start), by=0.1)
cope.q1.pred.interval <- predict(cope.lm, newdata=data.frame(cope.start=newx, weight = weights.lo.hi[[1]]), 
                                 interval="confidence", level = 0.95)
cope.q1.pred.interval <- as.data.frame(cope.q1.pred.interval)
cope.q1.pred.interval$cope.start <- newx
cope.q5.pred.interval <- predict(cope.lm, newdata=data.frame(cope.start=newx, weight = weights.lo.hi[[2]]), 
                                 interval="confidence", level = 0.95)
cope.q5.pred.interval <- as.data.frame(cope.q5.pred.interval)
cope.q5.pred.interval$cope.start <- newx
cope.q9.pred.interval <- predict(cope.lm, newdata=data.frame(cope.start=newx, weight = weights.lo.hi[[3]]), 
                                 interval="confidence", level = 0.95)
cope.q9.pred.interval <- as.data.frame(cope.q9.pred.interval)
cope.q9.pred.interval$cope.start <- newx
rm(newx)

cope.q1.pred.interval$lwr[which(cope.q1.pred.interval$lwr < -1)] <- -1
cope.q5.pred.interval$lwr[which(cope.q5.pred.interval$lwr < -1)] <- -1
cope.q9.pred.interval$lwr[which(cope.q9.pred.interval$lwr < -1)] <- -1

cope.1 <- ggplot(data, aes(x = cope.start, y = cope.e, fill = weight)) +
  stat_function(fun = cope.q1.curve, color = "dark goldenrod 2", size = 1, lty = 2, alpha = .6) +
  stat_function(fun = cope.q5.curve, color = "gray 40", size = 1, lty = 2, alpha = .6) +
  stat_function(fun = cope.q9.curve, color = "dark blue", size = 1, lty = 2, alpha = .6) +
  # stat_smooth(method="lm", se=FALSE, formula = y ~ x + I(x^2), color="red")+
  geom_point(pch = 21, size = pt.size, position = jitter, alpha = 0.9)+
  # geom_ribbon(inherit.aes = FALSE, data = cope.q1.pred.interval, aes(x = cope.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "dark goldenrod 2") +
  # geom_ribbon(inherit.aes = FALSE, data = cope.q5.pred.interval, aes(x = cope.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "gray 40") +
  # geom_ribbon(inherit.aes = FALSE, data = cope.q9.pred.interval, aes(x = cope.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "dark blue") +
  labs(x="Copepod initial frequency",y="Preference")+
  scale_fill_viridis_c("Salamander \nweight (g)", option = "E", direction = -1) +
  theme_classic()
# export 700*550

cope.2 <- ggplot(data, aes(x = weight, y = cope.e, color = factor(cope.start))) +
  stat_smooth(method="lm", se=FALSE, color="red")+
  geom_point(size = pt.size)+
  labs(x="Salamander weight (g)",y="Preference")+
  scale_color_viridis_d("Copepod \ninitial \nfrequency", option = "D", direction = 1) +
  theme_classic()

cope.3 <- ggplot(data, aes(x = cope.start, y = weight, color = cope.e)) +
  geom_point(size = pt.size)+
  labs(x="Copepod initial frequency",y="Salamander weight (g)")+
  scale_color_viridis_c("Preference", option = "B", direction = 1) +
  theme_classic()

# Ostracods ----

newx <- seq(min(ostr.start), max(ostr.start), by=0.1)
ostr.q1.pred.interval <- predict(ostr.lm, newdata=data.frame(ostr.start=newx, weight = weights.lo.hi[[1]]), 
                                 interval="confidence", level = 0.95)
ostr.q1.pred.interval <- as.data.frame(ostr.q1.pred.interval)
ostr.q1.pred.interval$ostr.start <- newx
ostr.q5.pred.interval <- predict(ostr.lm, newdata=data.frame(ostr.start=newx, weight = weights.lo.hi[[2]]), 
                                 interval="confidence", level = 0.95)
ostr.q5.pred.interval <- as.data.frame(ostr.q5.pred.interval)
ostr.q5.pred.interval$ostr.start <- newx
ostr.q9.pred.interval <- predict(ostr.lm, newdata=data.frame(ostr.start=newx, weight = weights.lo.hi[[3]]), 
                                 interval="confidence", level = 0.95)
ostr.q9.pred.interval <- as.data.frame(ostr.q9.pred.interval)
ostr.q9.pred.interval$ostr.start <- newx
rm(newx)

ostr.q1.pred.interval$lwr[which(ostr.q1.pred.interval$lwr < -1)] <- -1
ostr.q5.pred.interval$lwr[which(ostr.q5.pred.interval$lwr < -1)] <- -1
ostr.q9.pred.interval$lwr[which(ostr.q9.pred.interval$lwr < -1)] <- -1

ostr.coeffs <- ostr.lm$coefficients 
ostr.q1.curve <- function(x) ostr.coeffs[[1]] + ostr.coeffs[[2]]*weights.lo.hi[[1]] + ostr.coeffs[[3]]*x + 
  ostr.coeffs[[4]]*x^2 + ostr.coeffs[[5]]*weights.lo.hi[[1]]*x
ostr.q5.curve <- function(x) ostr.coeffs[[1]] + ostr.coeffs[[2]]*weights.lo.hi[[2]] + ostr.coeffs[[3]]*x + 
  ostr.coeffs[[4]]*x^2 + ostr.coeffs[[5]]*weights.lo.hi[[2]]*x
ostr.q9.curve <- function(x) ostr.coeffs[[1]] + ostr.coeffs[[2]]*weights.lo.hi[[3]] + ostr.coeffs[[3]]*x + 
  ostr.coeffs[[4]]*x^2 + ostr.coeffs[[5]]*weights.lo.hi[[3]]*x

ostr.1 <- ggplot(data, aes(x = ostr.start, y = ostr.e, fill = weight)) +
  # geom_ribbon(inherit.aes = FALSE, data = ostr.q1.pred.interval, aes(x = ostr.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "dark goldenrod 2") +
  # geom_ribbon(inherit.aes = FALSE, data = ostr.q5.pred.interval, aes(x = ostr.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "gray 40") +
  # geom_ribbon(inherit.aes = FALSE, data = ostr.q9.pred.interval, aes(x = ostr.start, ymin = lwr, ymax = upr), alpha = 0.1, fill = "dark blue") +
  stat_function(fun = ostr.q1.curve, color = "dark goldenrod 2", size = 1, lty = 2, alpha = .6) +
  stat_function(fun = ostr.q5.curve, color = "gray 40", size = 1, lty = 2, alpha = .6) +
  stat_function(fun = ostr.q9.curve, color = "dark blue", size = 1, lty = 2, alpha = .6) +
  # stat_smooth(method="lm", se=TRUE, formula = y ~ x + I(x^2), color="red")+
  # stat_smooth(data = subset(data, weight > mean(weight)), lty = 2,
  #             method="lm", se=FALSE, formula = y ~ x + I(x^2), color = "dark blue")+
  # stat_smooth(data = subset(data, weight < mean(weight)), lty = 2,
  #             method="lm", se=FALSE, formula = y ~ x + I(x^2), color = "dark goldenrod 2")+
  geom_point(pch = 21, size = pt.size, position = jitter, alpha = 0.9)+
  labs(x="Ostracod initial frequency",y="Preference")+
  scale_fill_viridis_c("Salamander \nweight (g)", option = "E", direction = -1) +
  theme_classic()

ostr.2 <- ggplot(data, aes(x = weight, y = ostr.e, color = factor(ostr.start), group = factor(ostr.start))) +
  stat_smooth(method="lm", se=FALSE, alpha = 0.2)+
  geom_point(size = pt.size)+
  labs(x="Salamander weight (g)",y="Preference")+
  scale_color_viridis_d("Ostracod \ninitial \nfrequency", option = "D", direction = 1) +
  theme_classic()

ostr.3 <- ggplot(data, aes(x = ostr.start, y = weight, color = ostr.e)) +
  geom_point(size = pt.size)+
  labs(x="Ostracod initial frequency",y="Salamander weight (g)")+
  scale_color_viridis_c("Preference", option = "B", direction = 1) +
  theme_classic()

# prepare each plot as a panel

clado_1<-arrangeGrob(clado.1, bottom = textGrob("(a)", x = unit(0, "npc")
                                              , y   = unit(2, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=11)))
cope_1<-arrangeGrob(cope.1, bottom = textGrob("(b)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
ostr_1<-arrangeGrob(ostr.1, bottom = textGrob("(c)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
clado_2<-arrangeGrob(clado.2, bottom = textGrob("(d)", x = unit(0, "npc")
                                              , y   = unit(2, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=11)))
cope_2<-arrangeGrob(cope.2, bottom = textGrob("(e)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
ostr_2<-arrangeGrob(ostr.2, bottom = textGrob("(f)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
clado_3<-arrangeGrob(clado.3, bottom = textGrob("(g)", x = unit(0, "npc")
                                               , y   = unit(2, "npc"), just=c("left","top"),
                                               gp=gpar(col="black", fontsize=11)))
cope_3<-arrangeGrob(cope.3, bottom = textGrob("(h)", x = unit(0, "npc")
                                             , y   = unit(2, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=11)))
ostr_3<-arrangeGrob(ostr.3, bottom = textGrob("(i)", x = unit(0, "npc")
                                             , y   = unit(2, "npc"), just=c("left","top"),
                                             gp=gpar(col="black", fontsize=11)))

# combine panels into single grid figure
grid.arrange(clado_1, cope_1, ostr_1)
grid.arrange(clado_1, cope_1, ostr_1, clado_2, cope_2, ostr_2, clado_3, cope_3, ostr_3)
# export in 1250 * 700

rm(pt.size, clado_1, cope_1, ostr_1, clado_2, cope_2, ostr_2, clado_3, cope_3, ostr_3)



### Functional Responses ################################################################

# check out predation on each pooled taxon ----

cladoplot <- ggplot(data, aes(x = clado.start, y = clado.eat)) + 
  geom_point(size = 2, alpha = 0.7) +
  ylab("Cladocerans eaten") +
  xlab("Initial frequency of cladocerans") +
  theme_classic()
cladoplot

copeplot <- ggplot(data, aes(x = cope.start, y = cope.eat)) + 
  geom_point() +
  theme_classic()
copeplot

ostrplot <- ggplot(data, aes(x = ostr.start, y = ostr.eat)) + 
  geom_point() +
  theme_classic()
ostrplot


# fit Type I predation response ----

# nllH1Binom <- function(pars, m, N){
#   # p <- a # get the predicted probability 
#   lL <- dbinom(x = m, prob = a, size = N, log=T) # for given attack rate, how many eaten?
#   nlL <- -sum(lL)
#   return(nlL)}
# 
# # estimate start parameter, a (attack rate)
# (cladostartpar1 <- ("a"=sum(data$clado.eat)/sum(data$clado.start)))
# (copestartpar1 <- ("a"=sum(data$cope.eat)/sum(data$cope.start)))
# (ostrstartpar1 <- ("a"=sum(data$ostr.eat)/sum(data$ostr.start)))
# 
# 
# # fit Type I curves
# (fit1.clado <- optim(fn = nllH1Binom, par = cladostartpar1, m=data$clado.eat, 
#                      N = data$clado.start,method="BFGS"))
# 
# (fit1.cope<- optim(fn = nllH1Binom, par = copestartpar1, m=data$cope.eat, 
#                    N = data$cope.start,method="BFGS"))
# 
# (fit1.ostr<- optim(fn = nllH1Binom, par = ostrstartpar1, m=data$ostr.eat, 
#                    N = data$ostr.start,method="BFGS"))
# 
# # negative log-likelihoods:
# 
# fit1.clado$value
# 
# fit1.cope$value
# 
# fit1.ostr$value
# 
# 
# # visualize Type I curves against data
# cladoplot + 
#   geom_smooth(method='lm',formula=y~x-1,se=FALSE) +
#   labs(x="Cladoceran initial frequency",y="Cladocerans eaten")+
#   ggtitle("Cladocerans - Type I response")
# 
# copeplot + 
#   geom_smooth(method='lm',formula=y~x-1,se=FALSE) +
#   labs(x="Copepod initial frequency",y="Copepods eaten")+
#   ggtitle("Copepods - Type I response")
# 
# ostrplot + 
#   geom_smooth(method='lm',formula=y~x-1,se=FALSE) +
#   labs(x="Ostracod initial frequency",y="Ostracods eaten")+
#   ggtitle("Ostracods - Type I response")

# using mle2

nllH1Binom2 <- function(a, m, N){
  # p <- a # get the predicted probability 
  lL <- dbinom(x = m, prob = a, size = N, log=T) # for given attack rate, how many eaten?
  nlL <- -sum(lL)
  return(nlL)}

# estimate start parameters - mle2 requires a list
cladostartpar1.mle2<-list("a" = sum(data$clado.eat)/sum(data$clado.start))
copestartpar1.mle2<-list("a" = sum(data$cope.eat)/sum(data$cope.start))
ostrstartpar1.mle2<-list("a" = sum(data$ostr.eat)/sum(data$ostr.start))


# fit Type I curves
(cladofit.mle2<-mle2(nllH1Binom2,start=cladostartpar1.mle2,
                     data=data.frame(N=data$clado.start, m=data$clado.eat),
                     optimizer="optimx"))

(copefit.mle2<-mle2(nllH1Binom2,start=copestartpar1.mle2,
                    data=data.frame(N=data$cope.start, m=data$cope.eat),
                    optimizer="optimx"))

(ostrfit.mle2<-mle2(nllH1Binom2,start=ostrstartpar1.mle2,
                    data=data.frame(N=data$ostr.start, m=data$ostr.eat),
                    optimizer="optimx"))

# fit Type II predation response ----

# holling2 <- function(N, a, h){
#   m <- a/(1 + a*N*h) # this is the type II fxn response
#   return(m)}
# 
# 
# nllH2Binom <- function(pars, m, N){
#   p <- holling2(N=N, a=pars['a'], h=pars['h']) # get the predicted probability 
#   lL <- dbinom(x = m, prob = p, size = N, log=T) # for given attack rate, how many eaten?
#   nlL <- -sum(lL)
#   return(nlL)}
# 
# 
# cladostartpars2 <- c("a" = sum(data$clado.eat)/sum(data$clado.start), 
#                      "h" = 1/100)
# cladostartpars2
# 
# copestartpars2 <- c("a" = sum(data$cope.eat)/sum(data$cope.start), 
#                     "h" = 1/270)
# copestartpars2
# 
# ostrstartpars2 <- c("a" = sum(data$ostr.eat)/sum(data$ostr.start), 
#                     "h" = 1/100)
# ostrstartpars2
# 
# # use optim()
# fit2.clado <- optim(fn = nllH2Binom, par = cladostartpars2, m=data$clado.eat, 
#                     N = data$clado.start, method="L-BFGS-B",lower=list(a=1e-10, h=1e-10), 
#                     upper=list(a=1, h=1e10))
# fit2.cope<- optim(fn = nllH2Binom, par = copestartpars2, m=data$cope.eat, 
#                   N = data$cope.start, method="L-BFGS-B",lower=list(a=1e-5, h=1e-15), 
#                   upper=list(a=1, h=1e10))
# fit2.ostr<-  optim(fn = nllH2Binom, par = ostrstartpars2, m=data$ostr.eat, 
#                    N = data$ostr.start, method="L-BFGS-B",lower=list(a=1e-5, h=1e-15), 
#                    upper=list(a=1, h=1e15))
# fit2.clado
# # nll = 151.8
# fit2.cope
# # nll = 942.7
# fit2.ostr
# # nll = 219.8
# # convergence error message


# use mle2() ----

# mle2 takes list of start parameters, not vector
cladostartpars.mle2 <- list("a" = sum(data$clado.eat)/sum(data$clado.start), 
                            "h" = 1/100)
copestartpars.mle2 <- list("a" = sum(data$cope.eat)/sum(data$cope.start), 
                           "h" = 1/270)
ostrstartpars.mle2 <- list("a" = sum(data$ostr.eat)/sum(data$ostr.start), 
                           "h" = 1/100)

holling2 <- function(N, a, h){
  m <- a/(1 + a*N*h) # this is the type II fxn response
  return(m)}


nllH2Binom2 <- function(a, h, m, N){
  p <- holling2(N=N, a=a, h=h) # get the predicted probability 
  lL <- dbinom(x = m, prob = p, size = N, log=T) # for given attack rate, how many eaten?
  nlL <- -sum(lL)
  return(nlL)}

cladofit1.mle2<-mle2(nllH2Binom2,start=cladostartpars.mle2,
                     data=data.frame(N=data$clado.start, m=data$clado.eat),
                     method="L-BFGS-B",lower=list(a=1e-10,h=1e-15),
                     upper=list(a=1,h=1e10))
cladofit1.mle2 # pretty close to optim() result

copefit1.mle2<-mle2(nllH2Binom2,start=copestartpars.mle2,
                    data=data.frame(N=data$cope.start, m=data$cope.eat),
                    method="L-BFGS-B",lower=list(a=1e-5,h=1e-15),
                    upper=list(a=1,h=1e10))
copefit1.mle2 # pretty close to optim() result

ostrfit1.mle2<-mle2(nllH2Binom2,start=ostrstartpars.mle2,
                    data=data.frame(N=data$ostr.start, m=data$ostr.eat),
                    method="L-BFGS-B",lower=list(a=1e-5,h=1e-15),
                    upper=list(a=1,h=1e10))
ostrfit1.mle2 # pretty close to optim() result
# still the warning about lack of convergence. 

# redefine the function b/c taking log functions on a & h,
# so transform them to keep within a sensible range
holling2bound <- function(N, a, h){
  a <- plogis(a)
  h <- exp(h)
  p <- a/(1 + a*h*N)
  return(p)}

# check math
curve(holling2(x, a=0.5, h=0.1), from=0, to=100)
curve(holling2bound(x, a=qlogis(0.5), h=log(0.1)), from=0, to=100,
      add=T, col=2, lty=2)

# now redefine nllH2Binom to use the new function.
nllH2boundBinom <- function(a, h, m, N){
  p <- (holling2bound(N=N, a=a, h=h)) ## get the predicted probability 
  lL <- dbinom(x = m, prob = p, size = N, log=T)
  nlL <- -sum(lL)
  return(nlL)}

# transform starting parameters
cladostartpars2.mle2 <- list(a=qlogis(cladostartpars.mle2$a), h=log(cladostartpars.mle2$h))
copestartpars2.mle2 <- list(a=qlogis(copestartpars.mle2$a), h=log(copestartpars.mle2$h))
ostrstartpars2.mle2 <- list(a=qlogis(ostrstartpars.mle2$a), h=log(ostrstartpars.mle2$h))

# calculate fits using mle2
(fit2.clado2 <- mle2(nllH2boundBinom, start=cladostartpars2.mle2,optimizer="optimx",
                     method="Nelder-Mead",
                     #method="L-BFGS-B",
                     #lower=list(a=1e-5, h=1e-5),upper=list(a=1, h=1e10),
                     data=data.frame(N=data$clado.start, m=data$clado.eat))) 
(fit2.cope2 <- mle2(nllH2boundBinom, start=copestartpars2.mle2, optimizer="optimx",
                    data=data.frame(N=data$cope.start, m=data$cope.eat),method='Nelder-Mead'
)) 
(fit2.ostr2 <- mle2(nllH2boundBinom, start=ostrstartpars2.mle2, optimizer="optimx",
                    data=data.frame(N=data$ostr.start, m=data$ostr.eat),method='Nelder'
)) 

# no more convergence problems

coef(fit2.clado2)
plogis(coef(fit2.clado2)['a'])
exp(coef(fit2.clado2)['h'])
fit2.clado$par # for comparison w/ optim()

coef(fit2.cope2)
plogis(coef(fit2.cope2)['a'])
exp(coef(fit2.cope2)['h'])
fit2.cope$par # for comparison w/ optim()

coef(fit2.ostr2)
plogis(coef(fit2.ostr2)['a'])
exp(coef(fit2.ostr2)['h'])
fit2.ostr$par # for comparison w/ optim()


fit2.clado2
fit2.cope2
fit2.ostr2

# obtain likelihood profiles, for CIs
(prof.clado2<-profile(fit2.clado2))
(prof.cope2<-profile(fit2.cope2))
(prof.ostr2<-profile(fit2.ostr2))

plot(prof.clado2)
plogis(confint(prof.clado2))
exp(confint(prof.clado2))

plogis(confint(prof.cope2))

plot(prof.ostr2)
confint(prof.ostr2)
plogis(confint(prof.ostr2))
exp(confint(prof.ostr2))

holling2N <- function(N, a, h){ 
  predN <- holling2bound(N=N, a=a, h=h)*N
  return(predN)}

cladoplot + stat_function(fun=holling2N, 
                          args=list(a=coef(fit2.clado2)["a"], h=coef(fit2.clado2)["h"])) +
  labs(x="Cladoceran initial frequency",y="Cladocerans eaten")#+
# ggtitle("Cladocerans - Type II response")
# Hessian is ill-behaved. This feels pretty judge-y. Thanks, RStudio.

copeplot + stat_function(fun=holling2N, 
                         args=list(a=coef(fit2.cope2)["a"], h=coef(fit2.cope2)["h"])) +
  labs(x="Copepod initial frequency",y="Copepods eaten")#+
# ggtitle("Copepods - Type II response")

ostrplot + stat_function(fun=holling2N, 
                         args=list(a=coef(fit2.ostr2)["a"], h=coef(fit2.ostr2)["h"])) +
  labs(x="Ostracod initial frequency",y="Ostracods eaten")#+
# ggtitle("Ostracods - Type II response")


# fit Type III predation response ----

# starting parameter estimates: a = fraction eaten, h = 1/max eaten
cladostartpars.mle2 <- list("a" = sum(data$clado.eat)/sum(data$clado.start), 
                            "h" = 1/100)

copestartpars.mle2 <- list("a" = sum(data$cope.eat)/sum(data$cope.start), 
                           "h" = 1/270)

ostrstartpars.mle2 <- list("a" = sum(data$cope.eat)/sum(data$cope.start), 
                           "h" = 1/100)

# starting parameter estimates transformed for sensible parameter ranges
cladostartparsB.mle2 <- list(a=qlogis(cladostartpars.mle2$a), h=log(cladostartpars.mle2$h))
copestartparsB.mle2 <- list(a=qlogis(copestartpars.mle2$a), h=log(copestartpars.mle2$h))
ostrstartparsB.mle2 <- list(a=qlogis(ostrstartpars.mle2$a), h=log(ostrstartpars.mle2$h))

holling3 <- function(N, a, h){
  c<- 1/h # c = 1/h = max#eat
  d<- c/(a) # d=c/a = max#eat/a
  mort <- c*N/(d^2 + N^2)
  return(mort)} # gives death rate @ given pop size & param values
# c = 1/h = max#eat, d=c/a = max#eat/a

# nlL eqn for holling3
nllH3Binom2 <- function(a, h, m, N){ # mle2() method
  p <- holling3(N=N, a=a, h=h) # get the predicted probability 
  lL <- dbinom(x = m, prob = p, size = N, log=T) # for given attack rate, how many eaten?
  nlL <- -sum(lL)
  return(nlL)}


# Holling Type III eqn w/ parameters bound to sensible ranges via transformation
holling3bound <- function(N, a, h){
  a <- plogis(a) # transforms to range (0 - 1)
  h <- exp(h) # transforms to range (0-infinity)
  c<- 1/h # c = 1/h = max#eat
  d<- c/(a) # d=c/a = max#eat/a
  mort <- c*N/(d^2 + N^2)
  return(mort)} # gives death rate @ given pop size & param values
# c = 1/h = max#eat, d=c/a = max#eat/a

# nlL eqn for holling3bound
nllH3boundBinom2 <- function(a, h, m, N){ # mle2() method
  p <- holling3bound(N=N, a=a, h=h) # get the predicted probability 
  lL <- dbinom(x = m, prob = p, size = N, log=T) # for given attack rate, how many eaten?
  nlL <- -sum(lL)
  return(nlL)}

#check math
curve(holling3(x, a=0.5, h=0.1), from=0, to=100)
curve(holling3bound(x, a=qlogis(0.5), h=log(0.1)), from=0, to=100,
      add=T, col=2, lty=2)

# cladocerans: fit holling3bound
fit3.clado <- mle2(nllH3Binom2, optimizer="optimx",
                   start = cladostartpars.mle2, data=data.frame(
                     N = data$clado.start, m=data$clado.eat))
fit3.clado
confint(fit3.clado)

fit3_mle.clado <- mle2(nllH3boundBinom2, optimizer="optimx",
                       start = cladostartparsB.mle2, data=data.frame(
                         N = data$clado.start, m=data$clado.eat))

fit3_mle.clado
coef(fit3_mle.clado)
plogis(coef(fit3_mle.clado)['a']) # a
exp(coef(fit3_mle.clado)['h']) # h
1/exp(coef(fit3_mle.clado)['h']) # c = 1/h
1/(exp(coef(fit3_mle.clado)['h'])*plogis(coef(fit3_mle.clado)['a'])) # d = c/a = 1/(h*a)
fit3.clado$par # for comparison w/ optim()
clado.confint<- confint(fit3_mle.clado)
plogis(clado.confint) #a
exp(clado.confint) #h
plot(clado.confint)

fit3.cope <- mle2(nllH3Binom2, optimizer="optimx",method="Nelder",
                  start = copestartpars.mle2, data=data.frame(
                    N = data$cope.start, m=data$cope.eat))
fit3.cope
confint(fit3.cope)
# no convergence

fit3_mle.cope <- mle2(nllH3boundBinom2, optimizer="optimx", method="Nelder",
                      start = copestartparsB.mle2, data=data.frame(
                        N = data$cope.start, m=data$cope.eat))

fit3_mle.cope

coef(fit3_mle.cope)
plogis(coef(fit3_mle.cope)['a']) # a
exp(coef(fit3_mle.cope)['h']) # h
1/exp(coef(fit3_mle.cope)['h']) # c
1/exp(coef(fit3_mle.cope)['h'])/plogis(coef(fit3_mle.cope)['a']) # d
fit3.cope$par # for comparison w/ optim()
cope.confint <- confint(fit3_mle.cope)
plogis(cope.confint) #a
exp(cope.confint) #h
profile(fit3_mle.cope)
plot(exp(profile(fit3_mle.cope)))

fit3.ostr <- mle2(nllH3Binom2, optimizer="optimx",method="Nelder",
                  start = ostrstartpars.mle2, data=data.frame(
                    N = data$ostr.start, m=data$ostr.eat))
fit3.ostr
confint(fit3.ostr)

fit3_mle.ostr <- mle2(nllH3boundBinom2, optimizer="optimx", method="Nelder-Mead",
                      start = ostrstartparsB.mle2, data=data.frame(
                        N = data$ostr.start, m=data$ostr.eat))
fit3_mle.ostr
coef(fit3_mle.ostr)
plogis(coef(fit3_mle.ostr)['a'])
exp(coef(fit3_mle.ostr)['h'])
1/exp(coef(fit3_mle.ostr)['h']) # c
1/exp(coef(fit3_mle.ostr)['h'])/plogis(coef(fit3_mle.ostr)['a']) # d
fit3.ostr$par # for comparison w/ optim()
ostr.confint <- confint(fit3_mle.ostr)
plogis(ostr.confint) #a
exp(ostr.confint) #h
plot(profile(fit3_mle.ostr))


holling3N <- function(N, a, h){ 
  predN <- holling3bound(N=N, a=a, h=h)*N
  return(predN)}

cladoplot + stat_function(fun=holling3N, 
                          args=list(a=coef(fit3_mle.clado)["a"], h=coef(fit3_mle.clado)["h"]), color = "blue") +
  labs(x="Cladoceran initial frequency",y="Cladocerans eaten")#+
# ggtitle("Cladocerans - Type III response")

copeplot + stat_function(fun=holling3N, 
                         args=list(a=coef(fit3_mle.cope)["a"], h=coef(fit3_mle.cope)["h"]), color = "blue") +
  labs(x="Copepod initial frequency",y="Copepods eaten")#+
# ggtitle("Copepods - Type III response")

ostrplot + stat_function(fun=holling3N, 
                         args=list(a=coef(fit3_mle.ostr)["a"], h=coef(fit3_mle.ostr)["h"]), color = "blue") +
  labs(x="Ostracod initial frequency",y="Ostracods eaten")#+
# ggtitle("Ostracods - Type III response")


nllH3boundBinom <- function(a, h, m, N){
  p <- holling3bound(N=N, a=a, h=h) # get the predicted probability 
  lL <- dbinom(x = m, size = N, prob = p, log=T) # for given attack rate, how many eaten?
  nlL <- -sum(lL)
  return(nlL)}

### functional response AICc values ----------

## log likelihoods
# mle2 does not produce subsettable output, so save log likelihoods manually

# Type 1

cladofit.mle2
clado1.nll <- 151.79

copefit.mle2
cope1.nll <- 940.91

ostrfit.mle2
ostr1.nll <- 360.58

# Type 2

fit2.clado2
clado2.nll <- 151.79

fit2.cope2
cope2.nll <- 940.91

fit2.ostr2
ostr2.nll <- 219.74

# Type 3

fit3_mle.clado
clado3.nll <- 148.95

fit3_mle.cope
cope3.nll <- 860.08

fit3_mle.ostr
ostr3.nll <- 406.2

## k = n terms for each model +1 for variance (intercept is always zero)
k.1<-1 +1
k.2<-2 +1
k.3<-2 +1


# AIC = 2* (negative log likelihood) + 2K
clado1.aic<-2*clado1.nll +2*k.1
cope1.aic<-2*cope1.nll +2*k.1
ostr1.aic<-2*ostr1.nll +2*k.1

clado2.aic<-2*clado2.nll +2*k.2
cope2.aic<-2*cope2.nll +2*k.2
ostr2.aic<-2*ostr2.nll +2*k.2

clado3.aic<-2*clado3.nll +2*k.3
cope3.aic<-2*cope3.nll +2*k.3
ostr3.aic<-2*ostr3.nll +2*k.3

print(c(clado1.aic,clado2.aic,clado3.aic))
print(c(cope1.aic,cope2.aic,cope3.aic))
print(c(ostr1.aic,ostr2.aic,ostr3.aic))

# AICc: AICc = AIC + 2K(K+1)/(n-K-1)
clado1.aicc<-clado1.aic+2*k.1*(k.1+1)/(20-k.1-1)
cope1.aicc<-cope1.aic+2*k.1*(k.1+1)/(20-k.1-1)
ostr1.aicc<-ostr1.aic+2*k.1*(k.1+1)/(20-k.1-1)

clado2.aicc<-clado2.aic+2*k.2*(k.2+1)/(20-k.2-1)
cope2.aicc<-cope2.aic+2*k.2*(k.2+1)/(20-k.2-1)
ostr2.aicc<-ostr2.aic+2*k.2*(k.2+1)/(20-k.2-1)

clado3.aicc<-clado3.aic+2*k.3*(k.3+1)/(20-k.3-1)
cope3.aicc<-cope3.aic+2*k.3*(k.3+1)/(20-k.3-1)
ostr3.aicc<-ostr3.aic+2*k.3*(k.3+1)/(20-k.3-1)

print(c(clado1.aicc,clado2.aicc,clado3.aicc))
print(c(cope1.aicc,cope2.aicc,cope3.aicc))
print(c(ostr1.aicc,ostr2.aicc,ostr3.aicc))


# Calculating R2 for best models ----

# R2 = 1 - SSres/SStot
# SSres = sum((obs-pred)^2)
# SStot = sum((obs-mean)^2)

# Cladocerans
clado.pred<-clado.start*holling3bound(a=coef(fit3_mle.clado)['a'],h=coef(fit3_mle.clado)['h'],N=clado.start)
clado.mean<-mean(clado.eat)

clado.SSres<-sum((clado.eat-clado.pred)^2)
clado.SStot<-sum((clado.eat-clado.mean)^2)
(clado.R2<-1-clado.SSres/clado.SStot)

clado.pred1<-clado.start*fit1.clado$par
clado.mean<-mean(clado.eat)

clado.SSres1<-sum((clado.eat-clado.pred1)^2)
(clado1.R2<-1-clado.SSres1/clado.SStot)


# Copepods
cope.pred<-cope.start*holling3bound(a=coef(fit3_mle.cope)['a'],h=coef(fit3_mle.cope)['h'],N=cope.start)
cope.mean<-mean(cope.eat)

cope.SSres<-sum((cope.eat-cope.pred)^2)
cope.SStot<-sum((cope.eat-cope.mean)^2)
(cope.R2<-1-cope.SSres/cope.SStot)

# Ostracods
ostr.pred<-ostr.start*holling2bound(a=coef(fit2.ostr2)['a'],h=coef(fit2.ostr2)['h'],N=ostr.start)
ostr.mean<-mean(ostr.eat)

ostr.SSres<-sum((ostr.eat-ostr.pred)^2)
ostr.SStot<-sum((ostr.eat-ostr.mean)^2)
(ostr.R2<-1-ostr.SSres/ostr.SStot)


ostr.pred3<-ostr.start*holling3bound(a=coef(fit3_mle.ostr)['a'],h=coef(fit3_mle.ostr)['h'],N=ostr.start)

ostr.SSres3<-sum((ostr.eat-ostr.pred3)^2)
(ostr3.R2<-1-ostr.SSres3/ostr.SStot)




### manuscript_fxn_response_figures.R ########################################

names(data)

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
plot(y~x,ylim=c(0,1),main="",
     xlab="Cypricercus initial frequency",ylab="Cypricercus per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}

avg.mean<-mean(total.eat[c(1,5,9,13,17)]-ost.l.eat[c(1,5,9,13,17)])
avg.sem<-sd(total.eat[c(1,5,9,13,17)]-ost.l.eat[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(total.eat[c(2,6,10,14,18)]-ost.l.eat[c(2,6,10,14,18)])
t1.sem<-sd(total.eat[c(2,6,10,14,18)]-ost.l.eat[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(total.eat[c(3,7,11,15,19)]-ost.l.eat[c(3,7,11,15,19)])
t2.sem<-sd(total.eat[c(3,7,11,15,19)]-ost.l.eat[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(total.eat[c(4,8,12,16,20)]-ost.l.eat[c(4,8,12,16,20)])
t3.sem<-sd(total.eat[c(4,8,12,16,20)]-ost.l.eat[c(4,8,12,16,20)])/sqrt(5)

y=c(t3.mean,avg.mean,t1.mean,t2.mean)
x=ost.l.start[c(4,1:3)]
sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
plot(y~x,ylim=c(0,450),xlim=c(0,200),main="",
     xlab="Cypricercus Initial Frequency",ylab="Other Zooplankton Consumed", 
     pch=16,bty="L")
for(i in 1:4){segments(x[i],y[i]-sem[i],x[i],y[i]+sem[i])}

rm(i,sem,x,y)
# export 700*350

library(bbmle)
library(optimx)

# Holling Type II eqn w/ parameters bound to sensible ranges via transformation
holling2bound <- function(N, a, h){
  a <- plogis(a)
  h <- exp(h)
  p <- a/(1 + a*h*N)
  return(p)}

# redefine nllH2Binom to use the new function.
nllH2boundBinom <- function(a, h, m, N){
  p <- (holling2bound(N=N, a=a, h=h)) ## get the predicted probability 
  lL <- dbinom(x = m, prob = p, size = N, log=T)
  nlL <- -sum(lL)
  return(nlL)}

ostrstartpars.mle2 <- list("a" = sum(data$ostr.eat)/sum(data$ostr.start), 
                           "h" = 1/100)
# transform start pars for bounded fxn
ostrstartpars2.mle2 <- list(a=qlogis(ostrstartpars.mle2$a), h=log(ostrstartpars.mle2$h))

(fit2.ostr2 <- mle2(nllH2boundBinom, start=ostrstartpars2.mle2, optimizer="optimx",
                    data=data.frame(N=data$ostr.start, m=data$ostr.eat),method='Nelder'
)) 

# fit Type III predation response ----

# starting parameter estimates: a = fraction eaten, h = 1/max eaten
cladostartpars.mle2 <- list("a" = sum(data$clado.eat)/sum(data$clado.start), 
                            "h" = 1/100)

copestartpars.mle2 <- list("a" = sum(data$cope.eat)/sum(data$cope.start), 
                           "h" = 1/270)

# starting parameter estimates transformed for sensible parameter ranges
cladostartparsB.mle2 <- list(a=qlogis(cladostartpars.mle2$a), h=log(cladostartpars.mle2$h))
copestartparsB.mle2 <- list(a=qlogis(copestartpars.mle2$a), h=log(copestartpars.mle2$h))


# Holling Type III eqn w/ parameters bound to sensible ranges via transformation
holling3bound <- function(N, a, h){
  a <- plogis(a) # transforms to range (0 - 1)
  h <- exp(h) # transforms to range (0-infinity)
  c<- 1/h # c = 1/h = max#eat
  d<- c/(a) # d=c/a = max#eat/a
  mort <- c*N/(d^2 + N^2)
  return(mort)} # gives death rate @ given pop size & param values
# c = 1/h = max#eat, d=c/a = max#eat/a

# nlL eqn for holling3bound
nllH3boundBinom2 <- function(a, h, m, N){ # mle2() method
  p <- holling3bound(N=N, a=a, h=h) # get the predicted probability 
  lL <- dbinom(x = m, prob = p, size = N, log=T) # for given attack rate, how many eaten?
  nlL <- -sum(lL)
  return(nlL)}

# fit holling3bound
fit3_mle.clado <- mle2(nllH3boundBinom2, optimizer="optimx",
                       start = cladostartparsB.mle2, data=data.frame(
                         N = data$clado.start, m=data$clado.eat))
fit3_mle.clado

fit3_mle.cope <- mle2(nllH3boundBinom2, optimizer="optimx", method="Nelder",
                      start = copestartparsB.mle2, data=data.frame(
                        N = data$cope.start, m=data$cope.eat))
fit3_mle.cope


# fxn for graphing Type II response curves
holling2N <- function(N, a, h){ 
  predN <- holling2bound(N=N, a=a, h=h)*N
  return(predN)}

# fxn for graphing Type III response curves
holling3N <- function(N, a, h){ 
  predN <- holling3bound(N=N, a=a, h=h)*N
  return(predN)}



cladoplot <- qplot(data=data, x=clado.start, y=clado.eat, geom=c('point'))
copeplot <- qplot(data=data, x=cope.start, y=cope.eat, geom=c('point'))
ostrplot <- qplot(data=data, x=ostr.start, y=ostr.eat, geom=c('point'))


# library(gridBase)


# Cladoceran figs
avg.mean<-mean(clado.eat[c(1,5,9,13,17)]/clado.start[c(1,5,9,13,17)])
avg.sem<-sd(clado.eat[c(1,5,9,13,17)]/clado.start[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(clado.eat[c(2,6,10,14,18)]/clado.start[c(2,6,10,14,18)])
t1.sem<-sd(clado.eat[c(2,6,10,14,18)]/clado.start[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(clado.eat[c(3,7,11,15,19)]/clado.start[c(3,7,11,15,19)])
t2.sem<-sd(clado.eat[c(3,7,11,15,19)]/clado.start[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(clado.eat[c(4,8,12,16,20)]/clado.start[c(4,8,12,16,20)])
t3.sem<-sd(clado.eat[c(4,8,12,16,20)]/clado.start[c(4,8,12,16,20)])/sqrt(5)

(clado.y<-c(t3.mean,avg.mean,t1.mean,t2.mean))
(clado.x<-clado.start[c(4,1:3)])
(clado.sem <- c(t3.sem,avg.sem,t1.sem,t2.sem))
(plot.dat<-data.frame(cbind(clado.x,clado.y,clado.sem)))
plot(clado.y~clado.x,ylim=c(min(clado.y-clado.sem),max(clado.y+clado.sem)),main="",
     xlab="Cladoceran initial frequency",ylab="Per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(clado.x[i],clado.y[i]-clado.sem[i],clado.x[i],clado.y[i]+clado.sem[i])}

clado.percap<-ggplot(plot.dat, aes(x=clado.x, y=clado.y))

clado1<-clado.percap+geom_point(size=2)+
  coord_cartesian(ylim = c(max(clado.y+clado.sem),0)) +
  labs(x="Cladoceran initial frequency", y="Per capita mortality")+
  geom_errorbar(ymin=clado.y-clado.sem, ymax=clado.y+clado.sem,width=2)+
  stat_function(fun=holling3bound,
                args=list(a=coef(fit3_mle.clado)["a"], h=coef(fit3_mle.clado)["h"]))+
  theme_classic()

clado2<-cladoplot + stat_function(fun=holling3N, 
                                  args=list(a=coef(fit3_mle.clado)["a"], h=coef(fit3_mle.clado)["h"])) +
  theme_classic()+
  labs(x="Cladoceran initial frequency",y="Cladocerans eaten")#+
# ggtitle("Cladocerans - Type III response")


# Copepod figs
avg.mean<-mean(cope.eat[c(1,5,9,13,17)]/cope.start[c(1,5,9,13,17)])
avg.sem<-sd(cope.eat[c(1,5,9,13,17)]/cope.start[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(cope.eat[c(2,6,10,14,18)]/cope.start[c(2,6,10,14,18)])
t1.sem<-sd(cope.eat[c(2,6,10,14,18)]/cope.start[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(cope.eat[c(3,7,11,15,19)]/cope.start[c(3,7,11,15,19)])
t2.sem<-sd(cope.eat[c(3,7,11,15,19)]/cope.start[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(cope.eat[c(4,8,12,16,20)]/cope.start[c(4,8,12,16,20)])
t3.sem<-sd(cope.eat[c(4,8,12,16,20)]/cope.start[c(4,8,12,16,20)])/sqrt(5)

cope.y<- c(t3.mean,avg.mean,t1.mean,t2.mean)
cope.x<-cope.start[c(4,1:3)]
cope.sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
plot.dat$cope.y<-cope.y
plot.dat$cope.x<-cope.x
plot.dat$cope.sem<-cope.sem
plot(cope.y~cope.x,ylim=c(min(cope.y-cope.sem),max(cope.y+cope.sem)),main="",
     xlab="Copepod initial frequency",ylab="Per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(cope.x[i],cope.y[i]-cope.sem[i],cope.x[i],cope.y[i]+cope.sem[i])}
plot.dat

cope.percap<-ggplot(plot.dat, aes(x=cope.x, y=cope.y))

cope1<-cope.percap+geom_point(size=2)+
  coord_cartesian(ylim = c(max(cope.y+cope.sem),0)) +
  labs(x="Copepod initial frequency", y="Per capita mortality")+
  geom_errorbar(ymin=cope.y-cope.sem, ymax=cope.y+cope.sem,width=2)+
  stat_function(fun=holling3bound,
                args=list(a=coef(fit3_mle.cope)["a"], h=coef(fit3_mle.cope)["h"]))+
  theme_classic()

cope2<-copeplot + stat_function(fun=holling3N, 
                                args=list(a=coef(fit3_mle.cope)["a"], h=coef(fit3_mle.cope)["h"])) +
  theme_classic()+
  labs(x="Copepod initial frequency",y="Copepods eaten")#+
# ggtitle("Copepods - Type III response")


# Ostracod figs
avg.mean<-mean(ostr.eat[c(1,5,9,13,17)]/ostr.start[c(1,5,9,13,17)])
avg.sem<-sd(ostr.eat[c(1,5,9,13,17)]/ostr.start[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(ostr.eat[c(2,6,10,14,18)]/ostr.start[c(2,6,10,14,18)])
t1.sem<-sd(ostr.eat[c(2,6,10,14,18)]/ostr.start[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(ostr.eat[c(3,7,11,15,19)]/ostr.start[c(3,7,11,15,19)])
t2.sem<-sd(ostr.eat[c(3,7,11,15,19)]/ostr.start[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(ostr.eat[c(4,8,12,16,20)]/ostr.start[c(4,8,12,16,20)])
t3.sem<-sd(ostr.eat[c(4,8,12,16,20)]/ostr.start[c(4,8,12,16,20)])/sqrt(5)

ostr.y<-c(t3.mean,avg.mean,t1.mean,t2.mean)
ostr.x<-ostr.start[c(4,1:3)]
ostr.sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
plot.dat$ostr.x<-ostr.x
plot.dat$ostr.y<-ostr.y
plot.dat$ostr.sem<-ostr.sem




plot(ostr.y~ostr.x,ylim=c(min(ostr.y-ostr.sem),max(ostr.y+ostr.sem)),main="",
     xlab="Ostracod initial frequency",ylab="Per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(ostr.x[i],ostr.y[i]-ostr.sem[i],ostr.x[i],ostr.y[i]+ostr.sem[i])}

ostr.percap<-ggplot(plot.dat, aes(x=ostr.x, y=ostr.y))

ostr1<-ostr.percap+geom_point(size=2)+
  coord_cartesian(ylim = c(max(ostr.y+ostr.sem),0)) +
  labs(x="Ostracod initial frequency", y="Per capita mortality")+
  geom_errorbar(ymin=ostr.y-ostr.sem, ymax=ostr.y+ostr.sem,width=2,
                position=position_dodge())+
  stat_function(fun=holling2bound,
                args=list(a=coef(fit2.ostr2)["a"], h=coef(fit2.ostr2)["h"]))+
  theme_classic()


ostr2<-ostrplot + stat_function(fun=holling2N, 
                                args=list(a=coef(fit2.ostr2)["a"], h=coef(fit2.ostr2)["h"])) +
  theme_classic()+
  labs(x="Ostracod initial frequency",y="Ostracods eaten")#+
# ggtitle("Ostracods - Type II response")




# Cladoceran 2 figs
avg.mean<-mean(clado.eat[c(1,5,9,13,17)])
avg.sem<-sd(clado.eat[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(clado.eat[c(2,6,10,14,18)])
t1.sem<-sd(clado.eat[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(clado.eat[c(3,7,11,15,19)])
t2.sem<-sd(clado.eat[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(clado.eat[c(4,8,12,16,20)])
t3.sem<-sd(clado.eat[c(4,8,12,16,20)])/sqrt(5)

(clado.y<-c(t3.mean,avg.mean,t1.mean,t2.mean))
(clado.x<-clado.start[c(4,1:3)])
(clado.sem <- c(t3.sem,avg.sem,t1.sem,t2.sem))
(plot.dat<-data.frame(cbind(clado.x,clado.y,clado.sem)))
plot(clado.y~clado.x,ylim=c(min(clado.y-clado.sem),max(clado.y+clado.sem)),main="",
     xlab="Cladocera initial frequency",ylab="Per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(clado.x[i],clado.y[i]-clado.sem[i],clado.x[i],clado.y[i]+clado.sem[i])}

clado.percap<-ggplot(plot.dat, aes(x=clado.x, y=clado.y))

# Copepod 2 figs
avg.mean<-mean(cope.eat[c(1,5,9,13,17)])
avg.sem<-sd(cope.eat[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(cope.eat[c(2,6,10,14,18)])
t1.sem<-sd(cope.eat[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(cope.eat[c(3,7,11,15,19)])
t2.sem<-sd(cope.eat[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(cope.eat[c(4,8,12,16,20)])
t3.sem<-sd(cope.eat[c(4,8,12,16,20)])/sqrt(5)

cope.y<- c(t3.mean,avg.mean,t1.mean,t2.mean)
cope.x<-cope.start[c(4,1:3)]
cope.sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
plot.dat$cope.y<-cope.y
plot.dat$cope.x<-cope.x
plot.dat$cope.sem<-cope.sem
plot(cope.y~cope.x,ylim=c(min(cope.y-cope.sem),max(cope.y+cope.sem)),main="",
     xlab="Copepoda initial frequency",ylab="Per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(cope.x[i],cope.y[i]-cope.sem[i],cope.x[i],cope.y[i]+cope.sem[i])}
plot.dat

cope.percap<-ggplot(plot.dat, aes(x=cope.x, y=cope.y))


# Ostracod 2 figs
avg.mean<-mean(ostr.eat[c(1,5,9,13,17)])
avg.sem<-sd(ostr.eat[c(1,5,9,13,17)])/sqrt(5)
t1.mean<-mean(ostr.eat[c(2,6,10,14,18)])
t1.sem<-sd(ostr.eat[c(2,6,10,14,18)])/sqrt(5)
t2.mean<-mean(ostr.eat[c(3,7,11,15,19)])
t2.sem<-sd(ostr.eat[c(3,7,11,15,19)])/sqrt(5)
t3.mean<-mean(ostr.eat[c(4,8,12,16,20)])
t3.sem<-sd(ostr.eat[c(4,8,12,16,20)])/sqrt(5)

ostr.y<-c(t3.mean,avg.mean,t1.mean,t2.mean)
ostr.x<-ostr.start[c(4,1:3)]
ostr.sem <- c(t3.sem,avg.sem,t1.sem,t2.sem)
plot.dat$ostr.x<-ostr.x
plot.dat$ostr.y<-ostr.y
plot.dat$ostr.sem<-ostr.sem




plot(ostr.y~ostr.x,ylim=c(min(ostr.y-ostr.sem),max(ostr.y+ostr.sem)),main="",
     xlab="Ostracod initial frequency",ylab="Per capita mortality", 
     pch=16,bty="L")
for(i in 1:4){segments(ostr.x[i],ostr.y[i]-ostr.sem[i],ostr.x[i],ostr.y[i]+ostr.sem[i])}

ostr.percap<-ggplot(plot.dat, aes(x=ostr.x, y=ostr.y))


# edit original fxn response plots to show means & sem's
cladoplot <- qplot(data=data, x=clado.start, y=clado.eat, geom=c('point'))
copeplot <- qplot(data=data, x=cope.start, y=cope.eat, geom=c('point'))
ostrplot <- qplot(data=data, x=ostr.start, y=ostr.eat, geom=c('point'))

library(ggplot2)
library(grid)
library(gridExtra)

clado2<-clado.percap+geom_point(size=2)+
  coord_cartesian(ylim = c(max(clado.y+clado.sem),0)) +
  labs(x="Cladoceran initial frequency", y="Cladocerans eaten")+
  geom_errorbar(ymin=clado.y-clado.sem, ymax=clado.y+clado.sem,width=2)+
  stat_function(fun=holling3N,
                args=list(a=coef(fit3_mle.clado)["a"], h=coef(fit3_mle.clado)["h"]))+
  theme_classic()

cope2<-cope.percap+geom_point(size=2)+
  coord_cartesian(ylim = c(max(cope.y+cope.sem),0)) +
  labs(x="Copepod initial frequency", y="Copepods eaten")+
  geom_errorbar(ymin=cope.y-cope.sem, ymax=cope.y+cope.sem,width=2)+
  stat_function(fun=holling3N,
                args=list(a=coef(fit3_mle.cope)["a"], h=coef(fit3_mle.cope)["h"]))+
  theme_classic()

ostr2<- ostr.percap+geom_point(size=2)+
  coord_cartesian(ylim = c(max(ostr.y+ostr.sem),0)) +
  geom_errorbar(ymin=ostr.y-ostr.sem, ymax=ostr.y+ostr.sem,width=2,
                position=position_dodge())+
  stat_function(fun=holling2N,
                args=list(a=coef(fit2.ostr2)["a"], h=coef(fit2.ostr2)["h"])) +
  theme_classic()+
  labs(x="Ostracod initial frequency",y="Ostracods eaten")

# add panel labels
clado1<-arrangeGrob(clado1, bottom = textGrob("(d)", x = unit(0, "npc")
                                              , y   = unit(2, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=11)))
cope1<-arrangeGrob(cope1, bottom = textGrob("(e)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
ostr1<-arrangeGrob(ostr1, bottom = textGrob("(f)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
clado2<-arrangeGrob(clado2, bottom = textGrob("(a)", x = unit(0, "npc")
                                              , y   = unit(2, "npc"), just=c("left","top"),
                                              gp=gpar(col="black", fontsize=11)))
cope2<-arrangeGrob(cope2, bottom = textGrob("(b)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))
ostr2<-arrangeGrob(ostr2, bottom = textGrob("(c)", x = unit(0, "npc")
                                            , y   = unit(2, "npc"), just=c("left","top"),
                                            gp=gpar(col="black", fontsize=11)))


grid.arrange(clado2,clado1,cope2,cope1,ostr2,ostr1)
# export in 800*XXX (after locking aspect ratio)