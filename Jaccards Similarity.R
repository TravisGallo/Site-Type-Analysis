###Estimating Jaccard's Similartity indices between sites#####
### BY Travis Gallo & Mason Fidino ###
## Code associated with a project we affectionatly call "Site Types" ##

# BETA DIVERSITY
# Create an array of jaccards indices by site, sample, season, and site comparison
# Pull out samples
zobs=model1$sims.list$z[,1:8,,]
# Set up things for loop
nspec=dim(zobs)[2]
nsamp=dim(zobs)[1]
forest.vec=which(st==1) # We are going to compare everything to forest preserves
nforest=length(forest.vec)
nseason=dim(zobs)[4]
# Rearrange array
zobs.new <- aperm(zobs, c(3,2,1,4))
# Loop using Jaccards simmilarity formula
Jsite=array(NA,dim=c(nsite,nsamp,nseason,nforest))
dim(Jsite)
for(k in 1:nsamp){
  for(i in 1:nsite){
    for(t in 1:nseason){
      for(f in 1:nforest){
        Jsite[i,k,t,f]=sum(zobs.new[forest.vec[f],,k,t]*zobs.new[i,,k,t])/(sum(zobs.new[forest.vec[f],,k,t])+sum(zobs.new[i,,k,t])-sum(zobs.new[i,,k,t]*zobs.new[forest.vec[f],,k,t]))
      }
    }
  }
}
# Process Data
cem.vec=which(st==3)
i=12
y=apply(Jsite[i,,1,],1,mean)
hist(y)
# Difference between Cemeteries and Forest Preserves
cem.j=array(NA,dim=c(nsamp,length(cem.vec),nseason))
for(i in 1:length(cem.vec)){
  for(t in 1:nseason){
    cem.j[,i,t]=apply(Jsite[cem.vec[i],,t,],1,mean)
  }
}
cem.mean=matrix(NA,nrow=nsamp,ncol=nseason)
for(t in 1:nseason){
  cem.mean[,t]=apply(cem.j[,,t],1,mean)
}
cem.quants=apply(cem.mean,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)
# Difference between City Parks and Forest Preserves
park.vec=which(st==2)
park.j=array(NA,dim=c(nsamp,length(park.vec),nseason))
for(i in 1:length(park.vec)){
  for(t in 1:nseason){
    park.j[,i,t]=apply(Jsite[park.vec[i],,t,],1,mean)
  }
}
park.mean=matrix(NA,nrow=nsamp,ncol=nseason)
for(t in 1:nseason){
  park.mean[,t]=apply(park.j[,,t],1,mean)
}
park.quants=apply(park.mean,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)
# Difference between Golf Courses and Forest Preserves
golf.vec=which(st==4)
golf.j=array(NA,dim=c(nsamp,length(golf.vec),nseason))
for(i in 1:length(golf.vec)){
  for(t in 1:nseason){
    golf.j[,i,t]=apply(Jsite[golf.vec[i],,t,],1,mean)
  }
}
golf.mean=matrix(NA,nrow=nsamp,ncol=nseason)
for(t in 1:nseason){
  golf.mean[,t]=apply(golf.j[,,t],1,mean)
}
golf.quants=apply(golf.mean,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)
