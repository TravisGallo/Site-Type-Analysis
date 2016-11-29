model{
  
  # Priors
  omega~dunif(0,1)
  
  #Priors for species by site effects in occupancy and detection
  # Intercepts can vary by season
  for(k in 1:M){
    for(t in 1:nseason){
      spe.gamma[k,t]~dnorm(0,spe.gam.tau[t]) 
      spe.phi[k,t]~dnorm(0,spe.phi.tau[t])
    }
  }
  # Hyper params that describe community
  
  for(k in 1:M){
    lpsi[k]~dnorm(0,tau.lpsi)
    theta.phi[k]~dnorm(0,tau.theta.phi)
    theta.gamma[k]~dnorm(0,tau.theta.gamma)
    for(i in 1:nsite){
      sp.st.gamma[i,k]~dnorm(mu.sp.st.gamma[st[i],k],tau.sp.st.gamma[st[i],k])
      sp.st.phi[i,k]~dnorm(mu.sp.st.phi[st[i],k],tau.sp.st.phi[st[i],k])
    }
  # For Detection
    alpha[k]~dnorm(0,tau.alpha)
    season[k,1]<-0
    for(t in 2:4){
      season[k,t]~dnorm(0,tau.season[t-1])
    }
  }
  
  # HyperPriors for global response
  # For Occupancy
  for(j in 1:4){
    mu.tau.gamma[j]~dgamma(3,2)
    sd.mu.st.gamma[j]<-1/(sqrt(mu.tau.gamma[j]))
    mu.tau.phi[j]~dgamma(3,2)
    sd.mu.st.phi[j]<-1/(sqrt(mu.tau.phi[j]))
    mu.global.gamma[j]~dnorm(0,0.33)
    mu.global.phi[j]~dnorm(0,0.33)
    for(k in 1:M){
      mu.sp.st.gamma[j,k]~dnorm(mu.global.gamma[j],mu.tau.gamma[j])
      mu.sp.st.phi[j,k]~dnorm(mu.global.phi[j],mu.tau.phi[j])
      tau.sp.st.gamma[j,k]~dgamma(3,2)
      sd.sp.st.gamma[j,k]<-1/(sqrt(tau.sp.st.gamma[j,k]))
      tau.sp.st.phi[j,k]~dgamma(3,2)
      sd.sp.st.phi[j,k]<-1/(sqrt(tau.sp.st.phi[j,k]))
    }
  }
  tau.lpsi~dgamma(3,2)
  sd.lpsi<-1/(sqrt(tau.lpsi))
  tau.theta.gamma~dgamma(3,2)
  sd.theta.gamma<-1/(sqrt(tau.theta.gamma))
  tau.theta.phi~dgamma(3,2)
  sd.theta.phi<-1/(sqrt(tau.theta.phi))
  
  # For intercept variance
  for(t in 1:nseason){
    spe.gam.sd[t]<-1/(sqrt(spe.gam.tau[t]))
    spe.gam.tau[t]~dgamma(3,2)
    spe.phi.sd[t]<-1/(sqrt(spe.phi.tau[t]))
    spe.phi.tau[t]~dgamma(3,2)
  }
  
  # For Detection
  tau.alpha~dgamma(3,2)
  sd.alpha<-1/(sqrt(tau.alpha))
  for(i in 1:3){
    tau.season[i]~dgamma(3,2)
    sd.season[i]<-1/(sqrt(tau.season[i]))
  }
  
  #Superpopulation process
  for(k in 1:M){
    w[k]~dbern(omega)
  }
  
  #Ecological Model
  for(k in 1:M){
    for(i in 1:nsite){
      z[k,i,1]~dbern(psi[k,i]*w[k])
      logit(psi[k,i])<-lpsi[k]
      for(t in 2:nseason){
        logit(pi[k,i,t])<-z[k,i,t-1]*(spe.phi[k,t] + sp.st.phi[i,k] + theta.phi[k]*pc[i]) +
          (1-z[k,i,t-1])*(spe.gamma[k,t] + sp.st.gamma[i,k] + theta.gamma[k]*pc[i])
        z[k,i,t]~dbern(pi[k,i,t]*w[k])
      }
    }
  }
  
  #Observational Model
  for(k in 1:M){
    for(t in 1:nseason){
      logit(p[k,t])<-alpha[k] + season[k,season_vec[t]]
      for(i in 1:nsite){
        mu.p[k,i,t]<-z[k,i,t]*p[k,t]
        y[k,i,t]~dbin(mu.p[k,i,t],J[i,t])
      }
    }
  }
  
  
  #Derived Params
  for(k in 1:M){
    Nocc.fs[k]<-sum(z[k,,]) # number of sites each species occupies
  }
  for(i in 1:nsite){
    for(t in 1:nseason){
      Nsite[i,t]<-sum(z[,i,t]) # number of species at each site each year
    }
  }
  n0<-sum(w[(nspec+1):(nspec+nz)]) # number of species not seen
  Ntotal<-sum(w[]) # metacommunity
  for(i in 1:nsite){
    for(t in 1:nseason){
      Nst[i,t]<-sum(z[,st[i],t]) # number of species at each site type
    }
  }
}