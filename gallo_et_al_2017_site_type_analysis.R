###############################################################
###### COMMUNITY MODEL USING DYNAMIC OCCUPANCY APPROACH #######
############ BY TRAVIS GALLO and MASON FIDINO #################
############### Urban Wildlife Institute ######################

# load required packages
package_load <- function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  pkgsToDownload <- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

package_load("jagsUI")

# load and format all data for JAGS model
# read in data
y_array <- readRDS("gallo_et_al_2017_y_array.RDS")
z_array <- readRDS("gallo_et_al_2017_z_array.RDS")
j_matrix <- readRDS("gallo_et_al_2017_j_matrix.RDS")


# read in the species names. Note: The 1st dimension of the y and z array
### are ordered this way. As is covariate data.
species <- read.table("gallo_et_al_2017_species_names.txt", header=TRUE)

# data Setup
# number of observed species
n.spec=nrow(species)

# set up data augmentation
# we add 8 unseen species as an informative prior on our community
nz <- 8
# m = metacommunity (8 observed species + 8 theoretical unobserved species)
m <- n.spec+nz
# create new y-array with augemented data
yaug <- array(0,dim=c(m, nrow(j_matrix),ncol(j_matrix)))
yaug[1:n.spec,,] <- y_array

# make the augmented data have the same pattern of missing data
missings <- is.na(yaug[1,,])
for(k in (n.spec+1):(n.spec+nz)){
  yaug[k,,][missings] <- NA
}

# number of sites and seasons
nsite <- dim(yaug)[2]
nseason <- dim(yaug)[3]

# set up site type covariates
covs <- read.table("gallo_et_al_2017_site_type.txt", header=TRUE)
park <- covs$park
forest <- covs$forest
cem <- covs$cemet
golf <- covs$golf

# create site type vector for indexing in model
st <- rep(0,dim(yaug)[2])
st[which(forest>0)] <- 1
st[which(park>0)] <- 2
st[which(cem>0)] <- 3
st[which(golf>0)] <- 4

# urban Covariate
pc <- read.table("gallo_et_al_2017_urban_cov.txt", header=TRUE)

# set up data for JAGS model
model.data <- list(y=yaug,
                nseason=dim(yaug)[3],
                nsite=dim(yaug)[2],
                J=as.matrix(j_matrix),
                nspec=n.spec,
                nz=nz,
                M=n.spec+nz,
                park=park,
                forest=forest,
                cem=cem,
                golf=golf,
                season_vec=c(3,4,1,2,3,4,1,2,3,4,1),
                st=st,
                pc=pc$pc1
                )

# set up initial values for JAGS model
wst <- rep(1,n.spec+nz)
zst <- array(1, dim=c(n.spec+nz, dim(yaug)[2],dim(yaug)[3]))
inits <- function() {list(z=zst,
                      w=wst,
                      lpsi=rnorm(n=n.spec+nz,sd=0.25),
                      spe.gamma=matrix(rnorm(n=nseason*m,sd=0.25),ncol=nseason,nrow=m),
                      spe.phi=matrix(rnorm(n=nseason*m,sd=0.25),ncol=nseason,nrow=m),
                      sp.st.gamma=matrix(rnorm(n=nsite*m,sd=0.25),ncol=m,nrow=nsite),
                      sp.st.phi=matrix(rnorm(n=nsite*m,sd=0.25),ncol=m,nrow=nsite),
                      theta=rnorm(n=n.spec+nz,sd=0.25),
                      alpha=rnorm(n=n.spec+nz,sd=0.25)
                      )}

# model parameters to save
params <- c("theta.phi","theta.gamma","spe.gam.sd","spe.phi.sd","spe.gamma","spe.phi",
         "sd.mu.st.gamma","sd.mu.st.phi","omega","mu.lpsi","sd.lpsi","sp.st.gamma",
         "sp.st.phi","mu.sp.st.gamma","mu.sp.st.phi","mu.st.gamma","mu.st.phi",
         "sd.sp.st.phi","sd.sp.st.gamma","mu.alpha","sd.alpha","Ntotal","Nsite",
         "mu.season","sd.alpha","z","theta","mu.global.phi","mu.global.gamma")

#MCMC settings
#this is running all but 2 cores on a computer. Must adjust "nc" according to number of cores.
ni <- 150000; nt <- 20; nb <- 50000; nc <- 6

#run JAGS
model1 <- jags(model.data, inits, params, "gallo_et_al_2017_site_type_model.R",n.chains=nc,n.thin=nt,n.iter=ni,n.burnin=nb, parallel=TRUE)

###############################
##### DERIVED  PARAMETERS #####
###############################

# BETA DIVERSITY
# Create an array of jaccards indices by site, sample, season, and site comparison
# Pull out samples
zobs <- model1$sims.list$z[,1:8,,]
# Set up things for loop
nspec <- dim(zobs)[2]
nsamp <- dim(zobs)[1]
forest.vec <- which(st==1) # We are going to compare everything to forest preserves
nforest <- length(forest.vec)
nseason <- dim(zobs)[4]
# Rearrange array
zobs.new <- aperm(zobs, c(3,2,1,4))
# Loop using Jaccards simmilarity formula
Jsite <- array(NA,dim=c(nsite,nsamp,nseason,nforest))
dim(Jsite)
for(k in 1:nsamp){
  for(i in 1:nsite){
    for(t in 1:nseason){
      for(f in 1:nforest){
      Jsite[i,k,t,f] <- sum(zobs.new[forest.vec[f],,k,t]*zobs.new[i,,k,t])/(sum(zobs.new[forest.vec[f],,k,t])+sum(zobs.new[i,,k,t])-sum(zobs.new[i,,k,t]*zobs.new[forest.vec[f],,k,t]))
      }
    }
  }
}
# Process Data
cem.vec <- which(st==3)
i <- 12
y <- apply(Jsite[i,,1,],1,mean)

# Difference between Cemeteries and Forest Preserves
cem.j <- array(NA,dim=c(nsamp,length(cem.vec),nseason))
for(i in 1:length(cem.vec)){
  for(t in 1:nseason){
    cem.j[,i,t] <- apply(Jsite[cem.vec[i],,t,],1,mean)
  }
}
cem.mean=matrix(NA,nrow=nsamp,ncol=nseason)
for(t in 1:nseason){
  cem.mean[,t] <- apply(cem.j[,,t],1,mean)
}
cem.quants <- apply(cem.mean,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)

# Difference between City Parks and Forest Preserves
park.vec <- which(st==2)
park.j <- array(NA,dim=c(nsamp,length(park.vec),nseason))
for(i in 1:length(park.vec)){
  for(t in 1:nseason){
    park.j[,i,t] <- apply(Jsite[park.vec[i],,t,],1,mean)
  }
}
park.mean <- matrix(NA,nrow=nsamp,ncol=nseason)
for(t in 1:nseason){
  park.mean[,t] <- apply(park.j[,,t],1,mean)
}
park.quants <- apply(park.mean,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)

# Difference between Golf Courses and Forest Preserves
golf.vec <- which(st==4)
golf.j <- array(NA,dim=c(nsamp,length(golf.vec),nseason))
for(i in 1:length(golf.vec)){
  for(t in 1:nseason){
    golf.j[,i,t] <- apply(Jsite[golf.vec[i],,t,],1,mean)
  }
}
golf.mean <- matrix(NA,nrow=nsamp,ncol=nseason)
for(t in 1:nseason){
  golf.mean[,t] <- apply(golf.j[,,t],1,mean)
}
golf.quants <- apply(golf.mean,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)

# SPECIES RICHNESS VALUES
Ntotal.mean <- model.vague$mean$Ntotal
Ntotal.sd <- model.vague$sd$Ntotal
Ntotal.median <- model.vague$q50$Ntotal
Ntotal.low <- model.vague$q2.5$Ntotal
Ntotal.hi <- model.vague$q97.5$Ntotal

# Variance to Mean Rations
# Global
med.VMR.global.gamma <- NULL
VMR.global.gamma.less <- NULL
VMR.global.gamma.great <- NULL
for(i in 1:4){
mean.global.gamma <- model1$sims.list$mu.global.gamma[,i]
var.global.gamma <- model1$sims.list$sd.mu.st.gamma[,i]^2
VMR.global.gamma <- var.global.gamma/abs(mean.global.gamma)
med.VMR.global.gamma[i] <- median(VMR.global.gamma)
VMR.global.gamma.great[i] <- length(which(VMR.global.gamma>1))/length(VMR.global.gamma)
VMR.global.gamma.less[i] <- length(which(VMR.global.gamma<1))/length(VMR.global.gamma)
}
med.VMR.global.phi <- NULL
VMR.global.phi.less <- NULL
VMR.global.phi.great <- NULL
for(i in 1:4){
  mean.global.phi <- model1$sims.list$mu.global.phi[,i]
  var.global.phi <- model1$sims.list$sd.mu.st.phi[,i]^2
  VMR.global.phi <- var.global.phi/abs(mean.global.phi)
  med.VMR.global.phi[i] <- median(VMR.global.phi)
  VMR.global.phi.great[i] <- length(which(VMR.global.phi>1))/length(VMR.global.phi)
  VMR.global.phi.less[i] <- length(which(VMR.global.phi<1))/length(VMR.global.phi)
}

table.VMR.global <- t(rbind(med.VMR.global.gamma,VMR.global.gamma.less,VMR.global.gamma.great,med.VMR.global.phi,VMR.global.phi.less,VMR.global.phi.great))

# Species specific
med.VMR.gamma <- matrix(NA,nrow=4,ncol=8)
VMR.gamma.less <- matrix(NA,nrow=4,ncol=8)
VMR.gamma.great <- matrix(NA,nrow=4,ncol=8)
for(j in 1:4){
  for(i in 1:8){
    mean.gamma <- model1$sims.list$mu.sp.st.gamma[,j,i]
    var.gamma <- model1$sims.list$sd.sp.st.gamma[,j,i]^2
    VMR.gamma <- var.gamma/abs(mean.gamma)
    med.VMR.gamma[j,i] <- median(VMR.gamma)
    VMR.gamma.less[j,i] <- length(which(VMR.gamma<1))/length(VMR.gamma)
    VMR.gamma.great[j,i] <- length(which(VMR.gamma>1))/length(VMR.gamma)
  }
}
med.VMR.phi <- matrix(NA,nrow=4,ncol=8)
VMR.phi.less <- matrix(NA,nrow=4,ncol=8)
VMR.phi.great <- matrix(NA,nrow=4,ncol=8)
for(j in 1:4){
  for(i in 1:8){
    mean.phi <- model1$sims.list$mu.sp.st.phi[,j,i]
    var.phi <- model1$sims.list$sd.sp.st.phi[,j,i]^2
    VMR.phi <- var.phi/abs(mean.phi)
    med.VMR.phi[j,i] <- median(VMR.phi)
    VMR.phi.less[j,i] <- length(which(VMR.phi<1))/length(VMR.phi)
    VMR.phi.great[j,i] <- length(which(VMR.phi>1))/length(VMR.phi)
  }
}  

# Deriving the 10 most common community make ups at each site type
get_com <- function(zobs = NULL, sv = NULL){
  coms <- apply(zobs[,,sv,], c(1,3,4), paste, collapse = "") 
  most_common <- sort(table(coms), decreasing = TRUE)[1:10]/prod(dim(coms))
  return(cbind(names(most_common), as.numeric(most_common)))
}

fcom <- get_com(zobs, forest.vec) #forest preserves
pcom <- get_com(zobs, park.vec) # city parks
ccom <- get_com(zobs, cem.vec) # cemeteries
gcom <- get_com(zobs, golf.vec) # golf courses