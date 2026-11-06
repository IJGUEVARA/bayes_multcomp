# Bivariate normal Dirichlet process mixture example.
# Run from the directory containing this script and partition_proposals.R.
project_dir <- "."

library(ggplot2)
library(dplyr)
library(mvtnorm)
library(MASS)

source(file.path(project_dir, "partition_proposals.R"))

# 1. Data simulation --------
# Component means by treatment (first index) and mixture component (second).
mu11 <- c(0, 0)
mu21 <- c(0, 0)
mu31 <- c(0, 0)

mu12 <- c(-2, -2)
mu22 <- c(1, 1)
mu32 <- c(-2, -2)

# Both components in all three treatments share this covariance matrix.
simulation_covariance <- matrix(c(1, 0.6, 0.6, 1), nrow = 2)
m <- length(mu11)

# SAMPLE SIZE

n1=48
n2=48
n3=48

n=n1+n2+n3

# Treatment design matrix
x=cbind(1,rep(c(1,2,3),c(n1,n2,n3)))
p=3 # Number of treatment groups

# Draw each treatment sample from an equally weighted two-component mixture.
y1=matrix(NA,n1,m)
y2=matrix(NA,n2,m)
y3=matrix(NA,n3,m)

U=runif(n1)
y1[U<=0.5,]=mvrnorm(sum(U<=0.5), mu11, simulation_covariance)
y1[U>0.5,]=mvrnorm(sum(U>0.5), mu12, simulation_covariance)

U=runif(n2)
y2[U<=0.5,]=mvrnorm(sum(U<=0.5),mu21,simulation_covariance)
y2[U>0.5,]=mvrnorm(sum(U>0.5),mu22,simulation_covariance)

U=runif(n3)
y3[U<=0.5,]=mvrnorm(sum(U<=0.5),mu31,simulation_covariance)
y3[U>0.5,]=mvrnorm(sum(U>0.5),mu32,simulation_covariance)

y0=rbind(y1,y2,y3)
y=scale(y0)
att=attributes(y)

# Evaluate the true mixture densities on the plotting grid.

xgrid = seq(-4, 4, length = 200)
ygrid = seq(-4, 4, length = 200)

xngrid=xgrid*att$`scaled:scale`[1]+att$`scaled:center`[1]
yngrid=ygrid*att$`scaled:scale`[2]+att$`scaled:center`[2]

gridxy = expand.grid(x = xgrid, y = ygrid)
gridxyn = expand.grid(x = xngrid, y = yngrid)

tdens1=0.5*dmvnorm(gridxyn,mu11,simulation_covariance)+0.5*dmvnorm(gridxyn,mu12,simulation_covariance)
tdens2=0.5*dmvnorm(gridxyn,mu21,simulation_covariance)+0.5*dmvnorm(gridxyn,mu22,simulation_covariance)
tdens3=0.5*dmvnorm(gridxyn,mu31,simulation_covariance)+0.5*dmvnorm(gridxyn,mu32,simulation_covariance)

# 2. MCMC sampler --------

# Hyperparameters

pro=1 # prior z and model size
M=3 # Dirichlet concentration

# Normal–inverse-Wishart hyperparameters
a0=m+2
Psi0=diag(m)
mu0=rep(0,m)
k0=1

# Gamma prior for the Dirichlet process concentration parameter
am=1 # SHAPE
bm=0.1 # RATE

# Initialize mixture allocations and the treatment partition.

ci=sample(1:M,size=n,replace = TRUE) # CLUSTER ALLOCATION MEMBERSHIPS
etzi=c(1,2,1) # Initial treatment partition
zi <- rep(etzi, c(n1, n2, n3))

# Allocate storage for the posterior summaries.

numMCMC <- 10000L
nburn <- round(numMCMC * 0.5) + 1L
# Accumulate retained predictive densities for the posterior mean contours.
density_sum1 <- numeric(nrow(gridxy))
density_sum2 <- numeric(nrow(gridxy))
density_sum3 <- numeric(nrow(gridxy))
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

  muni=array(NA,dim=c(m,1,K,L))
  nci=rep(NA,K)
  # Posterior mean precision and inverse-Wishart degrees of freedom
  kni=matrix(NA,nrow=K,ncol = L)
  ani=matrix(NA,nrow=K,ncol = L)
  # covariance matrices
  Psini=array(NA,dim=c(m,m,K,L))
  # Posterior predictive Student t scale matrices
  Sni=array(NA,dim=c(m,m,K,L))

  # Store unnormalized allocation weights, including a new mixture cluster.
  psi=matrix(0,ncol=K+1,nrow=n) # probs of cj for yi

  if (it >= nburn) {
    density1 <- density2 <- density3 <- numeric(nrow(gridxy))
  }

  for(i in 1:K)
  {
    nci[i]=sum(ci==labs[i])

    # Calculate predictive allocation weights under the current partition.
    for(z in unique(etzi))
    {
      nji=sum(ci==labs[i] & zi==z)
      yi=y[ci==labs[i] & zi==z,]
      myi=rep(0,m)
      if(nji>0){if(nji>1){myi=colMeans(yi)}else{myi=yi}}

      # Posterior parameters
      kni[i,z]=k0+nji
      ani[i,z]=a0+nji
      # COMPUTE SUM OF SQUARE MATRICES
      if(nji>1)
      {
        Caux=crossprod(scale(yi, scale = FALSE))
      }else{
        Caux=(yi-myi)%*%t(yi-myi)
      }

      if(nji==0){Caux=matrix(0,m,m)}
      Psini[,,i,z]=Psi0+Caux+(k0*nji*(myi-mu0)%*%t(myi-mu0))/(k0+nji)
      muni[,,i,z]=(k0*mu0+nji*myi)/(k0+nji)

      Sni[,,i,z]=(kni[i,z]+1)/((ani[i,z]-m+1)*kni[i,z])*Psini[,,i,z]

      # Student t posterior predictive density
      psi[zi==z,i]=nci[i]*dmvt(y[zi==z,], delta = muni[,,i,z], sigma = Sni[,,i,z], df = ani[i,z]-m+1,log = FALSE)
    }

    # Add the cluster contribution for each treatment.
    #xi=1
    if (it >= nburn) {
      z=unique(zi[x[,2]==1])
      density1=density1+nci[i]*dmvt(gridxy, delta = muni[,,i,z], sigma = Sni[,,i,z], df = ani[i,z]-m+1,log = F)
      z=unique(zi[x[,2]==2])
      density2=density2+nci[i]*dmvt(gridxy, delta = muni[,,i,z], sigma = Sni[,,i,z], df = ani[i,z]-m+1,log = F)
      z=unique(zi[x[,2]==3])
      density3=density3+nci[i]*dmvt(gridxy, delta = muni[,,i,z], sigma = Sni[,,i,z], df = ani[i,z]-m+1,log = F)
    }

  }

  nlab=max(ci)+1

  # PRED DENSITY
  Psin=Psi0*(k0+1)/(k0*(a0-m+1))
  psi[,K+1]=M*dmvt(y, delta = mu0, sigma = Psin, df = a0-m+1, log = FALSE)

  if (it >= nburn) {
    base_density <- M * dmvt(gridxy, delta = mu0, sigma = Psin, df = a0 - m + 1, log = FALSE)
    density_sum1 <- density_sum1 + (density1 + base_density) / (n + M)
    density_sum2 <- density_sum2 + (density2 + base_density) / (n + M)
    density_sum3 <- density_sum3 + (density3 + base_density) / (n + M)
  }

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
    liko <- multivariate_partition_log_likelihood(y, ci, zi, k0, a0, Psi0, mu0, m)
    likn <- multivariate_partition_log_likelihood(y, ci, zin, k0, a0, Psi0, mu0, m)
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

z_mat1 <- matrix(tdens1, nrow = 200, ncol = 200)
z_mat2 <- matrix(tdens2, nrow = 200, ncol = 200)
z_mat3 <- matrix(tdens3, nrow = 200, ncol = 200)

# Convert the average density back to the scale of the simulated responses.
density_divisor <- length(retained) * prod(att$`scaled:scale`)

z_est1 <- matrix(density_sum1 / density_divisor, nrow = 200, ncol = 200)
z_est2 <- matrix(density_sum2 / density_divisor, nrow = 200, ncol = 200)
z_est3 <- matrix(density_sum3 / density_divisor, nrow = 200, ncol = 200)

# Prepare the data for the faceted contour plot.
# combine all components
df_all <- bind_rows(
  contour_grid_data(xngrid, yngrid, z_mat1, "1", "True"),
  contour_grid_data(xngrid, yngrid, z_est1, "1", "Estimated"),
  contour_grid_data(xngrid, yngrid, z_mat2, "2", "True"),
  contour_grid_data(xngrid, yngrid, z_est2, "2", "Estimated"),
  contour_grid_data(xngrid, yngrid, z_mat3, "3", "True"),
  contour_grid_data(xngrid, yngrid, z_est3, "3", "Estimated")
)

df_all$density_group <- posterior_density_groups(postprob)[as.integer(df_all$comp)]

# Combine points
df_points <- bind_rows(
  data.frame(x = y1[,1], y = y1[,2], comp = "1"),
  data.frame(x = y2[,1], y = y2[,2], comp = "2"),
  data.frame(x = y3[,1], y = y3[,2], comp = "3")
)

# custom labels for facets
facet_labels <- c("1" = "j = 1", "2" = "j = 2", "3" = "j = 3")

pdens=ggplot() +
  geom_contour(
    data = df_all,
    aes(x, y, z = z, colour = density_group, linetype = type, group = type),
    bins = 8
  ) +
  geom_point(data = df_points, aes(x, y), size = 0.7) +
  facet_wrap(~comp, labeller = labeller(comp = facet_labels)) +

  scale_colour_discrete(name = "Density group") +
  scale_linetype_manual(
    name = "Density",
    values = c("True" = "solid", "Estimated" = "dashed")
  ) +

  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(size = 12, face = "bold"),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  ) +
  labs(x = NULL, y = NULL)
if (interactive()) print(pdens)
