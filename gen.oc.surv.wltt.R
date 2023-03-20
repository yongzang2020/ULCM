#' Obtain the operating characteristics of the umbrella design using the latent class model (ULCM  design) by simulating trials.
#'
#'
#' @param ncohort1 the total number of cohorts in stage I
#' @param ncohort2 the total number of cohorts in stage II
#' @param cohortsize the cohort size
#' @param maxt the maximum follow-up time 
#' @param accrual the accrual rate, i.e., the number of patients accrued in 1 unit of time
#' @int the time interval for the underlying piecewise exponential distribution of PFS
#' @p_B1 the prevalence of biomarker mutation on B1 
#' @p_B2 the prevalence of biomarker mutation on B2
#' @p the matrix of short-term response rate, row for subgroup and column for treatment
#' @lambda the matrix of the base rates for the underlying piecewise exponential distribution of PFS, row for 
#'         interval and column for treatments
#' @lambda_z the matrix of the coefficient of short-term responder rates for the underlying piecewise exponential  
#'           distribution of PFS, row for subgroup and column for treatment 
#' @tau the minimum difference achieved between treatment and control to calculate posterior probabilities Q_gt
#' @delta_s the cutoff to graduate a treatment for superiority.
#'          We recommend the default value of (\code{delta_s=0.99}) for general use.          
#' @delta_f the cutoff to graduate a treatment for futility.
#'          We recommend the default value of (\code{delta_s=0.05}) for general use.          
#' @delta_fn the cutoff to make final decision of efficacy. 
#'          We recommend the default value of (\code{delta_s=0.85}) for general use.          
#' @param ntrial the total number of trials to be simulated.
#' 
#' @return \code{get.oc()} returns the operating characteristics of the ULCM design as a data frame,
#'         including: (1) the short-term response rates in each subgroup by treatment (\code{StResp}),
#'         (2) the probabilities of survival in 12 months in each subgroup by treatment (\code{SurvProb}),
#'         (3) the percentages of the treatment selection in each subgroup (\code{SelPercent}),
#'         (4) the average number of patients (\code{TotalN}),
#'         (5) the average number of patients in each subgroup by treatment (\code{PtNum}),
#'         (6) the percentages of number of patients in each subgroup by treatment (\code{PtPercent}).   
               

Gen.OC.surv.wltt=function(ncohort1,ncohort2,cohortsize=10,accrual=3,maxt=12,int=c(6,6),p_B1=0.4,p_B2=0.5,
                          p,lambda,lambda.z,tau=0,delta_s,delta_f,delta_fn,ntrial=5000){
  
  #model with latent class M
  model_trt_ltt<-"model{
    for (i in 1:n0){
    z0[i]~dbern(p0)
    is.censored0[i]~dinterval(t0[i],t0.cen[i])
    t0[i]~dweib(alpha0,lambda0[i])
    lambda0[i]<-exp(beta.i0+beta.z0*z0[i])
    } 
    
    for (i in 1:n10){
    z10[i]~dbern(p10)
    is.censored10[i]~dinterval(t10[i],t10.cen[i])
    t10[i]~dweib(alpha1,lambda10[i])
    lambda10[i]<-exp(beta.i1+beta.z10*z10[i])
    } 
    for (i in 1:n11){
    z11[i]~dbern(p11)
    is.censored11[i]~dinterval(t11[i],t11.cen[i])
    t11[i]~dweib(alpha1,lambda11[i])
    lambda11[i]<-exp(beta.i1+beta.z11*z11[i])
    } 
    for (i in 1:n12){
    z12[i]~dbern(p12)
    is.censored12[i]~dinterval(t12[i],t12.cen[i])
    t12[i]~dweib(alpha1,lambda12[i])
    lambda12[i]<-exp(beta.i1+beta.z12*z12[i])
    } 
    for (i in 1:n13){
    z13[i]~dbern(p13)
    is.censored13[i]~dinterval(t13[i],t13.cen[i])
    t13[i]~dweib(alpha1,lambda13[i])
    lambda13[i]<-exp(beta.i1+beta.z13*z13[i])
    } 
    
    for (i in 1:n20){
    z20[i]~dbern(p20)
    is.censored20[i]~dinterval(t20[i],t20.cen[i])
    t20[i]~dweib(alpha2,lambda20[i])
    lambda20[i]<-exp(beta.i2+beta.z20*z20[i])
    } 
    for (i in 1:n21){
    z21[i]~dbern(p21)
    is.censored21[i]~dinterval(t21[i],t21.cen[i])
    t21[i]~dweib(alpha2,lambda21[i])
    lambda21[i]<-exp(beta.i2+beta.z21*z21[i])
    } 
    for (i in 1:n22){
    z22[i]~dbern(p22)
    is.censored22[i]~dinterval(t22[i],t22.cen[i])
    t22[i]~dweib(alpha2,lambda22[i])
    lambda22[i]<-exp(beta.i2+beta.z22*z22[i])
    } 
    for (i in 1:n23){
    z23[i]~dbern(p23)
    is.censored23[i]~dinterval(t23[i],t23.cen[i])
    t23[i]~dweib(alpha2,lambda23[i])
    lambda23[i]<-exp(beta.i2+beta.z23*z23[i])
    }
  
    p10=p_10
    p11=p_10*(M1==1)+p_11*(M1==2)+p_11*(M1==3)
    p12=p_10*(M1==1)+p_10*(M1==2)+p_12*(M1==3)
    p13=p_10*(M1==1)+p_11*(M1==2)+p_11*(M1==3)
    
    beta.z10=beta.z_10
    beta.z11=beta.z_10*(M1==1)+beta.z_11*(M1==2)+beta.z_11*(M1==3)
    beta.z12=beta.z_10*(M1==1)+beta.z_10*(M1==2)+beta.z_12*(M1==3)
    beta.z13=beta.z_10*(M1==1)+beta.z_11*(M1==2)+beta.z_11*(M1==3)
    
    p20=p_20
    p21=p_20*(M2==1)+p_20*(M2==2)+p_21*(M2==3)
    p22=p_20*(M2==1)+p_22*(M2==2)+p_22*(M2==3)
    p23=p_20*(M2==1)+p_22*(M2==2)+p_22*(M2==3)
    
    beta.z20=beta.z_20
    beta.z21=beta.z_20*(M2==1)+beta.z_20*(M2==2)+beta.z_21*(M2==3)
    beta.z22=beta.z_20*(M2==1)+beta.z_22*(M2==2)+beta.z_22*(M2==3)
    beta.z23=beta.z_20*(M2==1)+beta.z_22*(M2==2)+beta.z_22*(M2==3)
    
    M1~dcat(c(1/3, 1/3, 1/3))
    M2~dcat(c(1/3, 1/3, 1/3))
    
    alpha0~dgamma(0.001,0.001)T(0.01,1000)
    alpha1~dgamma(0.001,0.001)T(0.01,1000)
    alpha2~dgamma(0.001,0.001)T(0.01,1000)
    
    beta.i0~dnorm(0.0,0.01)
    beta.i1~dnorm(0.0,0.01)
    beta.i2~dnorm(0.0,0.01)
    
    p0~dbeta(0.5, 0.5)T(0.01,0.99)
    p_10~dbeta(0.5, 0.5)T(0.01,0.99)
    p_11~dbeta(0.5, 0.5)T(0.01,0.99)
    p_12~dbeta(0.5, 0.5)T(0.01,0.99)
    p_20~dbeta(0.5, 0.5)T(0.01,0.99)
    p_21~dbeta(0.5, 0.5)T(0.01,0.99)
    p_22~dbeta(0.5, 0.5)T(0.01,0.99)
    
    beta.z0~dnorm(0.0,0.01)
    beta.z_10~dnorm(0.0,0.01)
    beta.z_11~dnorm(0.0,0.01)
    beta.z_12~dnorm(0.0,0.01)
    beta.z_20~dnorm(0.0,0.01)
    beta.z_21~dnorm(0.0,0.01)
    beta.z_22~dnorm(0.0,0.01)
    
  }"
  
  
  #generate piece-wise exponential distribution;
  rpwexp <- function(n, rate, intervals=NULL, cumulative=FALSE){
    if(is.null(intervals)){
      if (cumulative){return(cumsum(rexp(n,rate[1])))}else
        return(rexp(n,rate[1]))}
    k <- length(rate)
    if (k==1){
      if(cumulative){return(cumsum(rexp(n,rate)))}else
        return(rexp(n,rate))
    }
    if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
    tx <- 0
    j <- 1
    times <- array(0,n)
    timex <- cumsum(intervals)
    indx <- array(TRUE,n)
    for(i in 1:k){
      nindx <- sum(indx)
      if (nindx==0) break
      increment <- rexp(nindx,rate[i])
      if (cumulative) times[indx] <- tx + cumsum(increment)
      else times[indx] <- tx + increment
      if (i<k){
        tx <- timex[i]
        indx <- (times > timex[i])
      }
    }
    return(times)
  }
  
  mcmc<-function(maxt,dur,enttime,trt,g,z,pfs){
    lof=pmin(dur-enttime,maxt)
    status=rep(1,length(pfs))
    status[which(pfs>lof)]=0
    pfs[which(pfs>lof)]=lof[which(pfs>lof)]
    t=pfs
    is.na(t)<-status==0
    is.censored<-1-status
    t.cen<-pfs+status
    #Create initial values
    tinits<-pfs+5
    is.na(tinits)<-status==1
    
    s0=trt==0
    s10=trt==1&g==0
    s11=trt==1&g==1
    s12=trt==1&g==2
    s13=trt==1&g==3
    s20=trt==2&g==0
    s21=trt==2&g==1
    s22=trt==2&g==2
    s23=trt==2&g==3
    
    #siminits<-list(t0=tinits[s0],t10=tinits[s10],t11=tinits[s11],t12=tinits[s12],t13=tinits[s13],
    #t20=tinits[s20],t21=tinits[s21],t22=tinits[s22],t23=tinits[s23])
    siminits<-list(p0=0.5,alpha0=1,beta.i0=0,beta.z0=0,
                   M1=1,p_10=0.5,p_11=0.5,p_12=0.5,alpha1=1,beta.i1=0,beta.z_10=0,beta.z_11=0,beta.z_12=0,
                   M2=1,p_20=0.5,p_21=0.5,p_22=0.5,alpha2=1,beta.i2=0,beta.z_20=0,beta.z_21=0,beta.z_22=0)
    
    simdata.temp<-list(n0=sum(s0),  t0=t[s0],  is.censored0=is.censored[s0],  t0.cen=t.cen[s0],  z0=z[s0],
                       n10=sum(s10),t10=t[s10],is.censored10=is.censored[s10],t10.cen=t.cen[s10],z10=z[s10],
                       n11=sum(s11),t11=t[s11],is.censored11=is.censored[s11],t11.cen=t.cen[s11],z11=z[s11],
                       n12=sum(s12),t12=t[s12],is.censored12=is.censored[s12],t12.cen=t.cen[s12],z12=z[s12],
                       n13=sum(s13),t13=t[s13],is.censored13=is.censored[s13],t13.cen=t.cen[s13],z13=z[s13],
                       n20=sum(s20),t20=t[s20],is.censored20=is.censored[s20],t20.cen=t.cen[s20],z20=z[s20],
                       n21=sum(s21),t21=t[s21],is.censored21=is.censored[s21],t21.cen=t.cen[s21],z21=z[s21],
                       n22=sum(s22),t22=t[s22],is.censored22=is.censored[s22],t22.cen=t.cen[s22],z22=z[s22],
                       n23=sum(s23),t23=t[s23],is.censored23=is.censored[s23],t23.cen=t.cen[s23],z23=z[s23]
    )
    
    simjags<-jags.model(textConnection(model_trt_ltt),data=simdata.temp,inits=siminits,quiet=T)
    update(simjags,2000, progress.bar="none")
    post.temp<- jags.samples(simjags, variable.names=c('p0','alpha0','beta.i0','beta.z0',
                                                       'p10','p11','p12','p13','alpha1','beta.i1','beta.z10','beta.z11','beta.z12','beta.z13',
                                                       'p20','p21','p22','p23','alpha2','beta.i2','beta.z20','beta.z21','beta.z22','beta.z23'), 
                             n.iter=2000,thin=2,progress.bar="none")
    #samples1<-coda.samples(simjags,c('p0','alpha0','beta.i0','beta.z0','p10','p11','p12','p13','alpha1','beta.i1','beta.z10','beta.z11','beta.z12','beta.z13','p20','p21','p22','p23','alpha2','beta.i2','beta.z20','beta.z21','beta.z22','beta.z23'),5000,thin=5)
    #summary(samples1)
    return(post.temp)
  }
  
  post.pfs=function(tau,post,trt,g){
    nmcmc=length(post$alpha0)
    alpha0=post[[1]]
    beta.i0=post[[4]]
    beta.z0=post[[7]]
    p0=post[[16]]
    alpha=post[[trt+1]]
    beta.i=post[[trt+4]]
    if(trt==1){beta.z=post[[g+8]];p=post[[g+17]]}
    if(trt==2){beta.z=post[[g+12]];p=post[[g+21]]}
    
    post.mean=p/((exp(beta.i+beta.z))^(1/alpha))*gamma(1+1/alpha)+(1-p)/((exp(beta.i))^(1/alpha))*gamma(1+1/alpha)
    post.mean0=p0/((exp(beta.i0+beta.z0))^(1/alpha0))*gamma(1+1/alpha0)+(1-p0)/((exp(beta.i0))^(1/alpha0))*gamma(1+1/alpha0)
    cond=(post.mean-post.mean0>tau)
    return(sum(cond)/nmcmc)
  }
  
  post.tprob=function(post){
    tprob=vector(mode = "list", length = 4)
    for (g in 0:3) {
      nmcmc=length(post$alpha0)
      alpha0=post[[1]]
      beta.i0=post[[4]]
      beta.z0=post[[7]]
      p0=post[[16]]
      alpha1=post[[2]]
      beta.i1=post[[5]]
      beta.z1=post[[g+8]]
      p1=post[[g+17]]
      alpha2=post[[3]]
      beta.i2=post[[6]]
      beta.z2=post[[g+12]]
      p2=post[[g+21]]
      post.med0=median(p0/((exp(beta.i0+beta.z0))^(1/alpha0))*gamma(1+1/alpha0)+(1-p0)/((exp(beta.i0))^(1/alpha0))*gamma(1+1/alpha0))
      post.med1=median(p1/((exp(beta.i1+beta.z1))^(1/alpha1))*gamma(1+1/alpha1)+(1-p1)/((exp(beta.i1))^(1/alpha1))*gamma(1+1/alpha1))
      post.med2=median(p2/((exp(beta.i2+beta.z2))^(1/alpha2))*gamma(1+1/alpha2)+(1-p2)/((exp(beta.i2))^(1/alpha2))*gamma(1+1/alpha2))
      tprob012=c(post.med0,post.med1,post.med2)/(post.med0+post.med1+post.med2)
      tprob01=c(post.med0,post.med1)/(post.med0+post.med1)
      tprob02=c(post.med0,post.med2)/(post.med0+post.med2)
      tprob[[g+1]]=list(tprob012=tprob012,tprob01=tprob01,tprob02=tprob02)
    }
    return(tprob)
  } 
  
  #simulate stage I patients: equal randomization
  stageI=function(ncohort1){
    
    for (i in 1:ncohort1){
      enttime.temp=rep(dur,cohortsize)
      dur=dur+accrual
      g.temp=sample(0:3, size=cohortsize, replace=TRUE, prob=c((1-p_B1)*(1-p_B2),p_B1*(1-p_B2),(1-p_B1)*p_B2,p_B1*p_B2))#assign subgroup 0:3
      trt.temp=sample(0:2, size=cohortsize, replace=TRUE, c(1/3,1/3,1/3))#assign treatment 0:2
      z.temp=rep(0,cohortsize)
      pfs.temp=rep(0,cohortsize)
      for (j in 1:cohortsize){
        z.temp[j]=rbinom(1,1,p[g.temp[j]+1,trt.temp[j]+1])
        rate1=lambda[1,trt.temp[j]+1]+lambda.z[g.temp[j]+1,trt.temp[j]+1]*z.temp[j]
        rate2=lambda[2,trt.temp[j]+1]+lambda.z[g.temp[j]+1,trt.temp[j]+1]*z.temp[j]
        rate=c(rate1,rate2,1/12)
        pfs.temp[j]=rpwexp(1,rate,int)
      }
      enttime=c(enttime,enttime.temp)
      g=c(g,g.temp)
      trt=c(trt,trt.temp)
      z=c(z,z.temp)
      pfs=c(pfs,pfs.temp)
    }
    
    return(list(dur=dur,enttime=enttime,trt=trt,g=g,z=z,pfs=pfs))
    
  }
  
  #simulate stage II patients: adaptive randomization
  stageII=function(ncohort1,ncohort2,delta_s,delta_f){ 
    dataI=stageI(ncohort1+1)
    dur=dataI$dur
    enttime=dataI$enttime
    g=dataI$g
    trt=dataI$trt
    z=dataI$z
    pfs=dataI$pfs
    ncohort2=ncohort2-1
    
    for (i in 1:ncohort2){
      
      post.temp=mcmc(maxt,dur,enttime,trt,g,z,pfs)
      tprob=post.tprob(post.temp)
      
      if (i==1) {
        for (m in 1:2){
          for (n in 0:3){
            prob=post.pfs(tau,post.temp,m,n)
            if (prob>delta_s) {decision[m,n+1]=1
            }else if (prob<delta_f) {decision[m,n+1]=-1
            }
          }
        }
        tg=(decision==0)
        
      } else {
        for (m in 1:2){
          for (n in 0:3){
            if (tg[m,n+1]) {
              prob=post.pfs(tau,post.temp,m,n)
              if (prob>delta_s) {decision[m,n+1]=1
              }else if (prob<delta_f) {decision[m,n+1]=-1}
              tg[m,n+1]=(decision[m,n+1]==0)
            }
          }
        }
        
      }
      
      if(sum(tg)==0) {stop=1;break}
      
      enttime.temp=rep(dur,cohortsize)
      dur=dur+accrual
      g.temp=sample(0:3, size=cohortsize, replace=TRUE, prob=c((1-p_B1)*(1-p_B2),p_B1*(1-p_B2),(1-p_B1)*p_B2,p_B1*p_B2))#assign subgroup 0:3
      trt.temp=rep(0,cohortsize)
      z.temp=rep(0,cohortsize)
      pfs.temp=rep(0,cohortsize)
      for (j in 1:cohortsize){
        tprob.temp=tprob[[g.temp[j]+1]]
        if (tg[1,g.temp[j]+1]&tg[2,g.temp[j]+1]) {trt.temp[j]=sample(0:2, size=1, replace=TRUE, tprob.temp$tprob012)
        }else if (tg[1,g.temp[j]+1]) {trt.temp[j]=sample(c(0,1), size=1, replace=TRUE, tprob.temp$tprob01)
        }else if (tg[2,g.temp[j]+1]) {trt.temp[j]=sample(c(0,2), size=1, replace=TRUE, tprob.temp$tprob02)
        }else {trt.temp[j]=99;next;}
        z.temp[j]=rbinom(1,1,p[g.temp[j]+1,trt.temp[j]+1])
        rate1=lambda[1,trt.temp[j]+1]+lambda.z[g.temp[j]+1,trt.temp[j]+1]*z.temp[j]
        rate2=lambda[2,trt.temp[j]+1]+lambda.z[g.temp[j]+1,trt.temp[j]+1]*z.temp[j]
        rate=c(rate1,rate2,1/12)
        pfs.temp[j]=rpwexp(1,rate,int)
      }
      
      enttime=c(enttime,enttime.temp[which(trt.temp!=99)])
      g=c(g,g.temp[which(trt.temp!=99)])
      z=c(z,z.temp[which(trt.temp!=99)])
      pfs=c(pfs,pfs.temp[which(trt.temp!=99)])
      trt=c(trt,trt.temp[which(trt.temp!=99)])
    }
    
    return(list(dur=dur,enttime=enttime,trt=trt,g=g,z=z,pfs=pfs,decision=decision,tg=tg))
  }
  
  
  
  set.seed(seed);
  Y = matrix(rep(0, 4 * ntrial), nrow=ntrial);
  N = matrix(rep(0, 12 * ntrial), nrow=ntrial);
  selpercent = matrix(rep(0, 4 * 4), nrow=4);
  SurvProb = matrix(rep(0, 3 * 4), nrow=3);
  npts = (ncohort1+ncohort2)*cohortsize;
  
  for(trial in 1:ntrial)
  {
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    setSessionTimeLimit(cpu = Inf, elapsed = Inf)
    
    dur=0 ## trial duration
    enttime=NULL ##entering time
    g=NULL ##subgroup
    trt=NULL ##treatment
    z=NULL ##z
    pfs=NULL ##pfs
    decision=matrix(rep(0, 8), ncol = 4) #decision matrix: row for treatment;column for subgroup
    stop=0
    
    data=stageII(ncohort1,ncohort2,delta_s,delta_f)
    decision=data$decision
    dur=data$dur
    enttime=data$enttime
    g=data$g
    trt=data$trt
    z=data$z
    pfs=data$pfs
    tg=data$tg
    
    if(sum(tg)>0){
      post=post.temp=mcmc(maxt,dur,enttime,trt,g,z,pfs)
      for (m in 1:2){
        for (n in 0:3){
          if (tg[m,n+1]) {
            prob=post.pfs(tau,post,m,n)
            if (prob>delta_fn) {decision[m,n+1]=1}
          }
        }
      }
    }
    
    for (i in 1:4){
      if (decision[1,i]==1 & decision[2,i]==1) {Y[trial,i]=1
      } else if (decision[1,i]==1 & decision[2,i]!=1) {Y[trial,i]=2
      } else if (decision[1,i]!=1 & decision[2,i]==1) {Y[trial,i]=3
      } else if (decision[1,i]!=1 & decision[2,i]!=1) {Y[trial,i]=4
      }
    }
    
    k=0
    for (i in 0:2){
      for (j in 0:3){
        k=k+1
        N[trial,k]=sum(trt==i&g==j)
      }
    }
    
  }
  
  for (i in 1:4) {
    selpercent[1,i]=length(which(Y[,i] == 1))/ntrial*100
    selpercent[2,i]=length(which(Y[,i] == 2))/ntrial*100
    selpercent[3,i]=length(which(Y[,i] == 3))/ntrial*100
    selpercent[4,i]=length(which(Y[,i] == 4))/ntrial*100
  }
  
  for (i in 1:3) {
    l1=lambda[1,i];l2=lambda[2,i];
    for (j in 1:4) {
      l1_sr=l1+lambda.z[j,i];l2_sr=l2+lambda.z[j,i];
      SurvProb[i,j]=round((1-p[j,i])*(exp(-6*l1)-exp(-6*l2)+exp(-12*l2))+p[j,i]*(exp(-6*l1_sr)-exp(-6*l2_sr)+exp(-12*l2_sr)),digits=3)
    }
  }
  
  
  out=list(StResp=t(p),SurvProb=SurvProb,N=N,Y=Y,SelPercent=selpercent,TotalN=mean(rowSums(N)),PtNum=matrix(colMeans(N),nrow=3,byrow=TRUE),PtPercent=matrix(colMeans(N/rowSums(N)*100),nrow=3,byrow=TRUE))
  return(out)
  
}


#example
ntrial=5000
ncohort1=25
ncohort2=25
delta_s=0.99
delta_f=0.05
delta_fn=0.85
p=matrix(c(0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6),nrow=4) #row for subgroup; column for treatment
lambda=matrix(c(1/13,1/11,1/13.2,1/11.2,1/13.1,1/11.1),nrow=2) #row for interval, column for treatments
lambda.z=matrix(c(-1/50,-1/50,-1/50,-1/50,-1/18,-1/18,-1/18,-1/18,-1/17.7,-1/17.7,-1/17.7,-1/17.7),nrow=4) #row for subgroup; column for treatment
Gen.OC.surv.wltt(ncohort1=ncohort1,ncohort2=ncohort2,p=p,lambda=lambda,lambda.z=lambda.z,delta_s=delta_s,delta_f=delta_f,delta_fn=delta_fn,ntrial=ntrial)
