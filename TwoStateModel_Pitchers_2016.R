
library(rjags)
library(lme4)
load("/Users/gregorymatthews/Dropbox/HMMbaseball/data2016pitchersPitchInfo.RData")
#163638 obervations


#Remove missing data 
data<-data[!is.na(data$elevation),] 
data<-data[!is.na(data$gtemp),]
data<-data[!is.na(data$start_speed_corr),]
data<-data[!is.na(data$batter_count),]
data<-data[!is.na(data$order_count),]
data<-data[!is.na(data$pitch_count),]
data<-data[!is.na(data$strike),]
data<-data[!is.na(data$park_id),]
data<-data[!is.na(data$on_1b),]
data<-data[!is.na(data$on_2b),]
data<-data[!is.na(data$on_3b),]
data<-data[!is.na(data$outs),]
data<-data[!is.na(data$inning),]
data<-data[!is.na(data$home_runs),]
data<-data[!is.na(data$away_runs),]
data<-data[!is.na(data$outs),]

#Remove pitcher effect
modPicher<-lmer(start_speed_corr ~ pitch_count+(1|mlbid),data=data)
data$resid1<-resid(modPicher)
#Remove Context effect
modContext<-lm(resid1 ~ elevation*gtemp + (!on_1b==0 | !on_2b==0 | !on_3b==0) ,data=data)
#modContext<-lm(resid1 ~ pi_pitch_type + elevation*gtemp + (!on_1b==0 | !on_2b==0 | !on_3b==0) ,data=data)
data$resid2<-resid(modContext)

#Order the pitches chronologically
data<-data[order(data$mlbid,data$UTC),]

#Add an indicator for new games.  
data$newGame<-c(1,diff(data$pitch_count)<0+0)

#Remove pitchers with fewer than 800 pitches.
keep<-rownames(table(data$mlbid))[(table(data$mlbid)>=800)]
data<-data[as.character(data$mlbid)%in%keep,] #98239 observations in 2016

#Pull out all the ids
ids<-sort(unique(data$mlbid)) #81 unique ids 
small<-data[data$mlbid%in%ids,c("mlbid","resid2","pitch_id","UTC","newGame","first","last")]

un<-sort(unique(small$mlbid))
small$ind <- NA
for (i in 1:length(un)){print(i)
  small$ind[small$mlbid==un[i]] <- i
}

small$UTC<-as.POSIXct(small$UTC)


#Get JAGS 4.x.y!  That's what makes it work!
small<-small[order(small$ind,small$UTC),]
n<-length(ids)
small<-small[small$ind<=n,]
mlbkey<-small[!duplicated(small$mlbid),c("mlbid","ind")]
ind<-c(0,cumsum(as.vector(table(small$ind))))
dat <- small$resid2

Sinit <- 0 

#Two state HMM model
model.str<-"model { 
for (j in 1:n){
dat[ind[j]+1] ~ dnorm(mu[ind[j]+1], sigma)
mu[ind[j]+1] <- (beta[2]+b[j,2])*(S[ind[j]+1]==1) + (beta[1]+b[j,1])*(S[ind[j]+1]==0)
S[ind[j]+1] ~ dbern(p[ind[j]+1])
logit(p[ind[j]+1]) <- (gamma[1]+g[j,1])*(Sinit==1) + (gamma[2]+g[j,2])*(Sinit==0)


for (t in (ind[j]+2):ind[j+1])
{
  dat[t] ~ dnorm(mu[t], sigma)
  mu[t] <- (beta[2]+b[j,2])*(S[t]==1) + (beta[1]+b[j,1])*(S[t]==0)
  S[t] ~ dbern(p[t])
  logit(p[t]) <- (gamma[1]+g[j,1])*(S[t-1]==1) + (gamma[2]+g[j,2])*(S[t-1]==0)
}

}

beta[1] ~ dunif(-100, 100)
beta[2] ~ dunif(beta[1], 100)

for (j in 1:n){
b[j,1] ~ dunif(-100, 100)
b[j,2] ~ dunif(b[j,1]+beta[1]-beta[2], 100)
g[j,1] ~ dunif(-10, 10)
g[j,2] ~ dunif(-10, 10)
}

for (i in 1:2){
gamma[i] ~ dunif(-10,10)
}

sigma ~ dunif(0,10000)


} "

setwd("/Users/gregorymatthews/Dropbox/HMMbaseball/")
write(model.str,"jags_model_twoStateModel_pitchers.bug")



greg<-function(seed){
  jags<-jags.model('/Users/gregorymatthews/Dropbox/HMMbaseball/jags_model_twoStateModel_pitchers.bug',data=list('dat'=dat, 'ind' = ind, 'Sinit' = Sinit, 'n' = n),n.chains=1,n.adapt=100, inits=list(.RNG.name="base::Wichmann-Hill", .RNG.seed=seed))
  print("adapted")
  print(seed)
  update(jags,500)
  out<-jags.samples(jags,c('beta','gamma','b','g','sigma',"S"),5000,thin=5)
  out
}
library(parallel)
start<-proc.time()
test<-mclapply(as.list(c(1:4)),greg,mc.cores = 4)
end<-proc.time()
end-start
save.image("/Users/gregorymatthews/Dropbox/HMMbaseball/results/twoStateModel_Pitchers_2016.RData")
