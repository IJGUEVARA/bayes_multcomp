# Univariate normal Dirichlet process mixture example.
# Run from the directory containing this script and partition_proposals.R.
project_dir <- "."

library(ggplot2)
library(tidyr)
library(dplyr)

source(file.path(project_dir, "partition_proposals.R"))

# 1. Data simulation --------

# MEANS AND STD DEV FOR TREATMENT/MIXTURE
mu11=0
mu21=1
mu31=0

mu12=1
mu22=-1
mu32=1

s11=0.2
s21=0.3
s31=0.2

s12=0.3
s22=0.4
s32=0.3

# SAMPLE SIZE

n1=200
n2=200
n3=200

n=n1+n2+n3

# Treatment design matrix
x=cbind(1,rep(c(1,2,3),c(n1,n2,n3)))
p=3 # Number of treatment groups

# Draw each treatment sample from an equally weighted two-component mixture.
y1=rep(NA,n1)
y2=rep(NA,n2)
y3=rep(NA,n3)

U=runif(n1)
y1[U<=0.5]=rnorm(sum(U<=0.5),mu11,sd=s11)
y1[U>0.5]=rnorm(sum(U>0.5),mu12,sd=s12)

U=runif(n2)
y2[U<=0.5]=rnorm(sum(U<=0.5),mu21,sd=s21)
y2[U>0.5]=rnorm(sum(U>0.5),mu22,sd=s22)

U=runif(n3)
y3[U<=0.5]=rnorm(sum(U<=0.5),mu31,sd=s31)
y3[U>0.5]=rnorm(sum(U>0.5),mu32,sd=s32)

y=c(y1,y2,y3)

# Evaluate the true mixture densities on the plotting grid.

ygrid=seq(-3,3,length.out=100) # grid
tdens1=0.5*dnorm(ygrid,mu11,sd=s11)+0.5*dnorm(ygrid,mu12,sd=s12)
tdens2=0.5*dnorm(ygrid,mu21,sd=s21)+0.5*dnorm(ygrid,mu22,sd=s22)
tdens3=0.5*dnorm(ygrid,mu31,sd=s31)+0.5*dnorm(ygrid,mu32,sd=s32)

# 2. MCMC sampler --------

# Hyperparameters

pro=1 # prior z and model size
M=3 # Dirichlet concentration

# HYPERPARAMETERS NORMAL-GAMMA
a0=2
b0=0.5
mu0=0
k0=0.1

# Gamma prior for the Dirichlet process concentration parameter
am=1 # SHAPE
bm=0.1 # RATE

# Initialize mixture allocations and the treatment partition.

ci=sample(1:M,size=n,replace = TRUE) # CLUSTER ALLOCATION MEMBERSHIPS
etzi=c(1,1,1) # Initial treatment partition
zi <- rep(etzi, c(n1, n2, n3))

# Allocate storage for the posterior summaries.

numMCMC <- 10000L
nburn <- round(numMCMC * 0.5) + 1L
# DENSITIES FOR TREATMENT 1, 2 AND 3
fyest1=matrix(0,nrow=numMCMC,ncol=100)
fyest2=matrix(0,nrow=numMCMC,ncol=100)
fyest3=matrix(0,nrow=numMCMC,ncol=100)
# SAVE zi's
sampleszi=matrix(NA,ncol=p,nrow=numMCMC)

# Run the MCMC chain.

for(it in 1:numMCMC)
{
  cat(it," of ",numMCMC,"\r")

  labs=unique(ci) # Occupied mixture labels
  K=length(labs) # Number of occupied mixture clusters
  L=length(unique(etzi)) # Number of treatment blocks

  # Record the same state used for the predictive densities below.
  sampleszi[it, ] <- etzi

  # Allocate posterior parameters for each mixture cluster and treatment block.

  muni=matrix(NA,nrow=K,ncol = L)
  nci=rep(NA,K)
  kni=matrix(NA,nrow=K,ncol = L)
  ani=matrix(NA,nrow=K,ncol = L)
  bni=matrix(NA,nrow=K,ncol = L)
  sni=matrix(NA,nrow=K,ncol = L)

  # Store unnormalized allocation weights, including a new mixture cluster.
  psi=matrix(0,ncol=K+1,nrow=n) # probs of cj for yi

  for(i in 1:K)
  {
    nci[i]=sum(ci==labs[i])

    # Calculate predictive allocation weights under the current partition.
    for(z in unique(etzi))
    {
      nji=sum(ci==labs[i] & zi==z)
      yi=y[ci==labs[i] & zi==z]
      myi=ifelse(nji>0,mean(yi),0)

      # Posterior parameters
      kni[i,z]=k0+nji
      ani[i,z]=a0+nji/2
      bni[i,z]=b0+0.5*sum((yi-myi)^2)+(k0*nji*(myi-mu0)^2)/(2*(k0+nji))
      muni[i,z]=(k0*mu0+nji*myi)/(k0+nji)

      sni[i,z]=bni[i,z]*(kni[i,z]+1)/(ani[i,z]*kni[i,z])

      # Student t posterior predictive density
      psi[zi==z,i]=nci[i]*dt((y[zi==z] - muni[i,z])/sqrt(sni[i,z]), df = 2*ani[i,z]) / sqrt(sni[i,z])
    }

    # Add the cluster contribution for each treatment.
    #xi=1
    z=unique(zi[x[,2]==1])
    fyest1[it,]=fyest1[it,]+nci[i]*dt((ygrid - muni[i,z])/sqrt(sni[i,z]), df = 2*ani[i,z]) / sqrt(sni[i,z])
    #xi=2
    z=unique(zi[x[,2]==2])
    fyest2[it,]=fyest2[it,]+nci[i]*dt((ygrid - muni[i,z])/sqrt(sni[i,z]), df = 2*ani[i,z]) / sqrt(sni[i,z])
    #xi=3
    z=unique(zi[x[,2]==3])
    fyest3[it,]=fyest3[it,]+nci[i]*dt((ygrid - muni[i,z])/sqrt(sni[i,z]), df = 2*ani[i,z]) / sqrt(sni[i,z])

  }

  nlab=max(ci)+1

  # PRED DENSITY
  sn=b0*(k0+1)/(a0*k0) # scale t pred distribution
  psi[,K+1]=M*dt((y - mu0)/sqrt(sn), df = 2*a0) / sqrt(sn)
  fyest1[it,]=fyest1[it,]+M*dt((ygrid - mu0)/sqrt(sn), df = 2*a0) / sqrt(sn)
  fyest2[it,]=fyest2[it,]+M*dt((ygrid - mu0)/sqrt(sn), df = 2*a0) / sqrt(sn)
  fyest3[it,]=fyest3[it,]+M*dt((ygrid - mu0)/sqrt(sn), df = 2*a0) / sqrt(sn)

  fyest1[it,]=fyest1[it,]/(n+M)
  fyest2[it,]=fyest2[it,]/(n+M)
  fyest3[it,]=fyest3[it,]/(n+M)

  # Update the mixture allocation of each observation.
  for(j in 1:n)
  {
    ci[j]=sample(c(labs,nlab),size = 1,prob = psi[j,])
  }

  # Draw a size change, then a local move conditional on that change.
  lold <- L
  moves <- partition_size_moves(lold, p)
  k <- moves[sample.int(length(moves), 1L)]
  lnew <- lold + k
  etzin <- propose_local_partition(etzi, k)

  # Expand the proposed treatment labels to observation-level labels.
  zin=rep(etzin,c(n1,n2,n3))

  # Condition both partition likelihoods on the newly updated allocations.
  if (!identical(canonicalize_partition(etzi), etzin)) {
    liko <- univariate_partition_log_likelihood(y, ci, zi, k0, a0, b0, mu0)
    likn <- univariate_partition_log_likelihood(y, ci, zin, k0, a0, b0, mu0)
    qprop <- proposal_log_probs(etzi, etzin, k)
    log_size_ratio <- log(size_move_prob(lold, lnew, p)) -
      log(size_move_prob(lnew, lold, p))
    log_alpha <- log_partition_prior(lnew, pro, p) -
      log_partition_prior(lold, pro, p) + likn - liko +
      qprop$log_hastings + log_size_ratio

    if (log(runif(1L)) < min(0, log_alpha)) {
      zi <- zin
      etzi <- etzin
    }
  }

  # Use the occupied cluster count from the updated allocations.
  K <- length(unique(ci))
  # Update the Dirichlet process concentration with the auxiliary beta draw.
  phi=rbeta(1,M+1,n)
  val=(am+K-1)/(n*(bm-log(phi)))
  probpi=val/(1+val)
  M=ifelse(runif(1)<=probpi,
           rgamma(1,shape = am+K,rate = bm-log(phi)),
           rgamma(1,shape = am+K-1,rate = bm-log(phi)))
}

# 3. Model table with posterior probabilities --------

# Estimate partition probabilities from the retained iterations.
retained <- nburn:numMCMC
model <- apply(sampleszi[retained, , drop = FALSE], 1L, paste, collapse = ",")
postprob <- sort(prop.table(table(model)), decreasing = TRUE)
print(postprob)

# 4. Plots --------

fnew1 <- fyest1[retained, , drop = FALSE]
fnew2 <- fyest2[retained, , drop = FALSE]
fnew3 <- fyest3[retained, , drop = FALSE]

cint1 <- density_quartiles(fnew1)
cint2 <- density_quartiles(fnew2)
cint3 <- density_quartiles(fnew3)

# Each row of a retained density matrix is one draw, evaluated at all grid points.
# The interval endpoints are pointwise posterior quartiles.

df <- data.frame(
  y = ygrid,
  true1 = tdens1,
  true2 = tdens2,
  true3 = tdens3,
  est1 = colMeans(fnew1),
  est2 = colMeans(fnew2),
  est3 = colMeans(fnew3),
  lower1 = cint1[,1],
  upper1 = cint1[,2],
  lower2 = cint2[,1],
  upper2 = cint2[,2],
  lower3 = cint3[,1],
  upper3 = cint3[,2]
)

# Reshape data to long format
df_long <- df %>%
  pivot_longer(
    cols = starts_with(c("true","est","lower","upper")),
    names_to = c(".value","level"),
    names_pattern = "([a-z]+)([1-3])"
  )

df_long$density_group <- posterior_density_groups(postprob)[as.integer(df_long$level)]
pdens <- ggplot(df_long, aes(x = y, group = level)) +
  # Pointwise 50% credible intervals
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = density_group), alpha = 0.5) +
  # True densities
  geom_line(aes(y = true, color = density_group, linetype = "True"), linewidth = 0.5) +
  # Estimated densities
  geom_line(aes(y = est, color = density_group, linetype = "Estimate"), linewidth = 0.5) +
  # Scale settings
  scale_color_discrete(name = "Density group") +
  scale_fill_discrete(name = "Density group") +
  scale_linetype_manual(name = "Type", values = c("True"="solid","Estimate"="dashed")) +
  theme_minimal() +
  theme(panel.grid = element_blank())+
  labs(x = "y", y = "Posterior Predictive Density", title = " ")

if (interactive()) print(pdens)
