#----------------------------------#
# Simulation Study: RRmix Model    #
# Stephen Salerno Jr. (ss2658)     #
# Last Update: September 23, 2016  #
#----------------------------------#

#----------------------------#
# Inverse-Gamma Optimization # 
#----------------------------#


load("AX-RRmix.RData")  # Change Path Based on Directory

sample_sig2g <- result.AX.RRmix$sig2_g

hist(sqrt(sample_sig2g))

IGFunction <- function(x,A,B){(B^A)/gamma(A) * x^(-A-1) * exp(-B/x)}

library(MASS)

IGParams <- fitdistr(sample_sig2g, IGFunction, list(A=1, B=1))$estimate


#------------------#
# is.wholenumber() #
#------------------#


is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){   # From integer{base} help file
  
  abs(x - round(x)) < tol

}


#---------------------#
# Simulation Function #
#---------------------#


simRRmix <- function(nsims=1, n=1, G=1, A=3, B=1, p=0.05, psi=1.5, QC=NULL,
                     
                     sig20=1, sig21=0.1, trmt=0.5, mu=0, q=0, Lam=NA){
  
  
  #-----------------------------------------------------------------------------#
  #                                                                             #
  #  Function to simulate 'nsims' data sets of size 'n' x 'G data from the      #
  #  RRmix Model Distributional Assumptions.                                    #
  #                                                                             #
  #  nsims - Number of simulated data sets generated        - Positive Integer  #
  #  n     - Number of observations per simulated data set  - Positive Integer  #
  #  G     - Number of genes per simulated data set         - Positive Integer  #
  #  A     - Shape parameter for Inverse-Gamma Distribution - [1, inf)          #
  #  B     - Scale parameter for Inverse-Gamma Distribution - [1, inf)          #
  #  p     - Proportion of genes differentially expressed   - [0, 1]            #
  #  psi   - Second mean value for Bg ~ Multivariate Normal - [0, inf)          #
  #  QC    - Proportion of genes set as quality controls    - [0, 1]            #
  #  sig20 - First variance component for Bg ~ MV Normal    - [0, inf)          #
  #  sig21 - Second variance component for Bg ~ MV Normal   - [0, inf)          #
  #  trmt  - Proportion of observations in treatment group  - [0, 1]            #
  #  mu    - Overall mean abundance                         - (-inf, inf)       #
  #  q     - Number of latent factors                       - Positive Integer  #
  #  Lam   - Loading matrix of size n x q                   - (-inf, inf)       #
  #                                                                             #
  #-----------------------------------------------------------------------------#
  
  
  ## ASSERTIONS ##
  
  
  if (mode(c(nsims, n, G, A, B, p, psi, QC, sig20,            # Assert parameters are numeric
             sig21, trmt, q, mu)) != "numeric"){            
    
    stop("Argument is not numeric")
    
  }
  
  if ((nsims<=0) | (n<=0) | (G<=0) |                          # Assert parameters are positive
      (A<=0) | (B<=0) | (p<0 | p>1) |
      (psi<0) | (sig20<0) | (sig21<0) |
      (trmt<0 | trmt>1) | (q<0)){                          
    
    stop("Argument out of range")
  
  }
  
  integer.check <- is.wholenumber(c(nsims, n, G, q))          # Check if nsims, n, G, & q are integers
  
  if (all.equal(integer.check, c(T,T,T,T)) != TRUE){          # Assert all are integers
    
    stop("nsims, n, or G are not all integers")
    
  }
  
  if ((q==0) & any(!is.na(Lam))){                                  # Assertions on Lam
    
    stop("Set number of factors to ncol(Lam)")
    
  }
  
  if ((q>0) & (mode(Lam) != "numeric")){            
    
    stop("Lam is not a numeric matrix")
    
  }
  
  if ((q > 0) & (class(Lam) != "matrix")){
    
    stop("Lam is not a numeric matrix")
    
  }
  
  if (!is.null(QC)){
    
    if ((QC < 0) | (QC > 1-p)){                                  # Assertion on QC
      
      stop("Argument out of range:  QC < 0 or QC > 1 - p")
      
    }
    
  }

  
  ## SIMULATE DATA ##
  
  
  library(MCMCpack)                                           # Package for rinvgamma() function
  
  library(MASS)                                               # Package for mvrnorm() function
  
  trmt.vect <- c(rep(1, n*trmt),                              # Treatment Status
                 rep(0, n*(1-trmt)))
  
  diff.vects <- as.list(rep(NA, nsims))                       # Container for latent indicators
  
  QC.vects <- as.list(rep(NA, nsims))                         # Container for quality controls
  
  sets <- as.list(rep(NA, nsims))                             # Container for simulated data sets
  
  signal.mats <- as.list(rep(NA, nsims))                      # Container for true biological signal
  
  factor.mats <- as.list(rep(NA, nsims))                      # Container for latent variation
    
  noise.mats  <- as.list(rep(NA, nsims))                      # Container for random noise
    
  for (i in 1:nsims) {                                        # For-loop to generate each data set
    
    diff.genes <- rep(0, G)                                   # Gene-specific latent indicators 

    set <- matrix(nrow=n, ncol=G)                             # Initialize each data set as matrix
    
    mu.set <- matrix(nrow=n, ncol=G)                          # Initialize overall mean component
    
    XB.set <- matrix(nrow=n, ncol=G)                          # Initialize biological signal component
    
    LamF.set <- matrix(nrow=n, ncol=G)                        # Initialize latent factor component
    
    W.set <- matrix(nrow=n, ncol=G)                           # Initialize random noise component
    
    for (j in 1:G) {                                          # For-loop for gene-specific observations
      
      sig2g <- rinvgamma(1, A, B)                             # Simulate gene-specific error variance
      
      bg <- rbinom(1, 1, p)                                   # Simulate gene-specific latent indicator
      
      diff.genes[j] <- bg                                     # Store gene-specific latent indicators 
      
      mu.Bg <- c(0, bg*psi)                                   # Mean vector for simulated Bg
      
      Sigma.Bg <- ((1-bg)*matrix(c(sig20,0,0,0),2,2)) +       # Covariance matrix for simulated Bg
                  ((bg)*matrix(c(sig20,0,0,sig21),2,2))
      
      Bg <- mvrnorm(1, mu.Bg, Sigma.Bg)                       # Simulated Bg vector
      
      X <- as.matrix(cbind(rep(1,n),                          # Simulated X matrix
                           c(rep(1, n*trmt),
                             rep(0, n*(1-trmt)))))              
      
      XBg <- X%*%Bg                                           # Gene indicator term in model
      
      XB.set[,j] <- XBg                                       # Save gene-specific signal to matrix
      
      if (q > 0){                                             # Simulate latent factors/loadings
        
        Fg <- rnorm(q)                                        # Simulated latent factors 
          
        LamFg <- Lam%*%Fg                                     # Latent term in model

        LamF.set[,j] <- LamFg                                 # Save gene-specific latent var. to matrix
        
      }
                         
      Wg <- rnorm(n, mean=0, sd=sig2g)                        # Simulate residual terms

      W.set[,j] <- Wg                                         # Save random noise to matrix
      
      mu.vect <- rep(mu, n)                                   # Overall mean vector

      mu.set[,j] <- mu.vect                                   # Save mean vector to matrix

      if(!is.null(QC)){                                       # "Spike-In" Quality Controls
        
        G.QC <- round(QC*G) 
        
        QC.index <- which(diff.genes == 0)[1:G.QC]
        
        XB.set[,QC.index] <- 0
        
      }

      if (q > 0){
        
        set <- mu.set + XB.set + LamF.set + W.set               # Create Simulated Data With Loadings
        
      }
      
      else{
        
        set <- mu.set + XB.set + W.set                          # Create Simulated Data Without Loadings
        
      }
      
      
      
    }
    
    diff.vects[[i]] <- diff.genes                             # Store differential genes for set
    
    if (!is.null(QC)){
      
      QC.vects[[i]]   <- QC.index                             # Store quality controls for set 
      
    }
    
    set <- data.frame(set)                                    # Convert data set to data frame
        
    sets[[i]] <- set                                          # Append data frame to list
    
    signal.mats[[i]] <- XB.set                                # Append biological signal to list
    
    factor.mats[[i]] <- LamF.set                              # Append latent variation to list
    
    noise.mats[[i]]  <- W.set                                 # Append random noise to list

    
 }
 
 result <- list(Treatment.Groups = trmt.vect,                 # Store Results
                Differential.Compounds = diff.vects,
                Quality.Controls = QC.vects,
                Simulated.Data  = sets,
                Biological.Signal = signal.mats,
                Latent.Variation = factor.mats,
                Random.Noise = noise.mats)
 
 return(result)                                               # Return Results
  
}




#------------------------------------------#
# Small Data Simulation - 0 Latent Factors #
#------------------------------------------#


set.seed (1212)


simulations <- simRRmix(nsims=50, n=6, G=265, A=3, B=1, p=0.05, psi=0.5, QC=0.05,
                        sig20=.55, sig21=.23, trmt=0.5, mu=13, q=0, Lam=NA)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)



#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)


## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)


## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)


## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=2)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)

## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)

## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 0)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)



## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)


## FDR Averages


FAMT.NBF.FDR.avmat <- FAMT.NBF.FDR[[1]]

for (i in (2:length(FAMT.NBF.FDR))){
  
  FAMT.NBF.FDR.avmat <- cbind(FAMT.NBF.FDR.avmat, FAMT.NBF.FDR[[i]])
  
}

FAMT.NBF.FDR.avg <- rowMeans(FAMT.NBF.FDR.avmat, na.rm=T)


## null Averages (By Rank)
 

FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, 
                               FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]], decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors=4)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)

FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)


## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)


## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, 
                               FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]], decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)


## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)



## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]], 
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, 
                                PCA_limma.null[[i]][order(PCA_limma.null[[i]], 
                                                          decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, 
                                RUV_limma.null[[i]][order(RUV_limma.null[[i]], 
                                                          decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)


#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


library(ggplot2)

dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p6x265.0F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("6 x 265 - 0 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))


p6x265.0F


png(filename="6x265_0F_ROC.png")

p6x265.0F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d6x265.0F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("6 x 265 - 0 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d6x265.0F

png(filename="6x265_0F_DET.png", width=1000, height=500)

d6x265.0F

dev.off()


#------------------------------------------#
# Small Data Simulation - 4 Latent Factors #
#------------------------------------------#

set.seed (1212)

Lam <- matrix(c(sample(rep(c(1, -1),each=6), 12),
                sample(rep(c(1, -1),each=6), 12),
                sample(rep(c(1, -1),each=6), 12),
                sample(rep(c(1, -1),each=6), 12)),
              ncol=4)

Lam <- Lam + rnorm(dim(Lam)[1]*dim(Lam)[2], 0, 0.1)


simulations <- simRRmix(nsims=50, n=12, G=265, A=3, B=1, p=0.05, psi=0.5, QC=0.05,
                        sig20=.55, sig21=.23, trmt=0.5, mu=13, q=4, Lam=Lam)



#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)

## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)

## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)


## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=4)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 

 
#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 4)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors=8)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)

## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)


## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)


## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)


## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]], 
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]], 
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]],
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)



#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p12x265.4F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("12 x 265 - 4 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p12x265.4F

png(filename="12x265_4F_ROC.png")

p12x265.4F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d12x265.4F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("12 x 265 - 4 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d12x265.4F

png(filename="12x265_4F_DET.png", width=1000, height=500)

d12x265.4F

dev.off()


#-------------------------------------------------#
# Medium-Small Data Simulation - 0 Latent Factors #
#-------------------------------------------------#

set.seed (1212)

simulations <- simRRmix(nsims=50, n=50, G=265, A=3, B=1, QC=0.05,     # Simulate Data
                        p=0.05, psi=0.5, sig20=.55, sig21=.23, 
                        trmt=0.5, mu=13, q=0, Lam=NULL)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)



#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)

## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)

## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]],
                                                                    decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=2)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 0)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors = 4)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- 1 - fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)



## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)


## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]],
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]],
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]],
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)



#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p50x265.0F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("50 x 265 - 0 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p50x265.0F

png(filename="50x265_0F_ROC.png")

p50x265.0F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d50x265.0F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("50 x 265 - 0 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d50x265.0F

png(filename="50x265_0F_DET.png", width=1000, height=500)

d50x265.0F

dev.off()
 


#-------------------------------------------------#
# Medium-Small Data Simulation - 4 Latent Factors #
#-------------------------------------------------#

set.seed (1212)

Lam <- matrix(c(sample(rep(c(1, -1),each=50), 100),
                sample(rep(c(1, -1),each=50), 100),
                sample(rep(c(1, -1),each=50), 100),
                sample(rep(c(1, -1),each=50), 100)),
              ncol=4)

Lam <- Lam + rnorm(dim(Lam)[1]*dim(Lam)[2], 0, 0.1)

simulations <- simRRmix(nsims=50, n=100, G=265, A=3, B=1, QC=0.05,     # Simulate Data
                        p=0.05, psi=0.5, sig20=.55, sig21=.23, 
                        trmt=0.5, mu=13, q=4, Lam=Lam)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)


#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){

    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)


## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)

## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=4)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 4)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]], 
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors = 8)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- 1 - fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)

## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)


## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)


## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]], 
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray", 
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]],
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]],
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]],
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)



#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p100x265.4F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("100 x 265 - 4 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p100x265.4F

png(filename="100x265_4F_ROC.png")

p100x265.4F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d100x265.4F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("100 x 265 - 4 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d100x265.4F

png(filename="100x265_4F_DET.png", width=1000, height=500)

d100x265.4F

dev.off()



#-------------------------------------------#
# Medium Data Simulation - 0 Latent Factors #
#-------------------------------------------#

set.seed (1212)

simulations <- simRRmix(nsims=50, n=100, G=265, A=3, B=1, QC=0.05,     # Simulate Data
                        p=0.05, psi=0.5, sig20=.55, sig21=.23, 
                        trmt=0.5, mu=13, q=0, Lam=NULL)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)


#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)

## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)

## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)


## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=2)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 0)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)

  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2,
                        test=2,
                        maxnbfactors = 4)$optimalnbfactors
  
    fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf = nbf.FAMT)
    
    FAMT.DEF.adjtest <- fit.FAMT.DEF$adjtest
    
    FAMT.DEF.adjpval <- fit.FAMT.DEF$adjpval
  
  FAMT.DEF.Fstats[[i]] <- FAMT.DEF.adjtest
  
  FAMT.DEF.null[[i]]   <- 1 - FAMT.DEF.adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)


## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)



## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]], 
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]],
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]],
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)



#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p100x265.0F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("100 x 265 - 0 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p100x265.0F

png(filename="100x265_0F_ROC.png")

p100x265.0F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d100x265.0F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("100 x 265 - 0 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d100x265.0F

png(filename="100x265_0F_DET.png", width=1000, height=500)

d100x265.0F

dev.off()


 
#-------------------------------------------#
# Medium Data Simulation - 2 Latent Factors #
#-------------------------------------------#

set.seed (1212)

Lam <- matrix(c(sample(rep(c(1, -1),each=100), 200),
                sample(rep(c(1, -1),each=100), 200),
                sample(rep(c(1, -1),each=100), 200),
                sample(rep(c(1, -1),each=100), 200)),
              ncol=4)

Lam <- Lam + rnorm(dim(Lam)[1]*dim(Lam)[2], 0, 0.1)

simulations <- simRRmix(nsims=50, n=200, G=265, A=3, B=1, QC=0.05,     # Simulate Data
                        p=0.05, psi=0.5, sig20=.55, sig21=.23, 
                        trmt=0.5, mu=13, q=4, Lam=Lam)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)



#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)

## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)

## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)


## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))

RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=4)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 4)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)



#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors = 8)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- 1 - fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)

## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)


## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)


## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]],
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]], 
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]], 
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)



#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)

dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p200x265.4F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("200 x 265 - 4 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p200x265.4F

png(filename="200x265_4F_ROC.png")

p200x265.4F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d200x265.4F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("200 x 265 - 4 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d200x265.4F

png(filename="200x265_4F_DET.png", width=1000, height=500)

d200x265.4F

dev.off()



#------------------------------------------#
# Large Data Simulation - 0 Latent Factors #
#------------------------------------------#

set.seed (1212)

simulations <- simRRmix(nsims=50, n=100, G=500, A=3, B=1, QC=0.05,     # Simulate Data
                        p=0.05, psi=0.5, sig20=.55, sig21=.23, 
                        trmt=0.5, mu=13, q=0, Lam=NULL)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)



#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect
  
  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)


## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect
  
  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)


## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=2)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 
  
  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))
  
  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect
  
  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 0)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect
  
  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors = 4)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- 1 - fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 
  
  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))
  
  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect
  
  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)

## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)


## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)


## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]],
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]],
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}


RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
  
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]],
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)


#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p100x500.0F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("100 x 500 - 0 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p100x500.0F

png(filename="100x500_0F_ROC.png")

p100x500.0F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d100x500.0F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("100 x 500 - 0 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d100x500.0F

png(filename="100x500_0F_DET.png", width=1000, height=500)

d100x500.0F

dev.off()


#------------------------------------------#
# Large Data Simulation - 2 Latent Factors #
#------------------------------------------#

set.seed (1212)

Lam <- matrix(c(sample(rep(c(1, -1),each=100), 200),
                sample(rep(c(1, -1),each=100), 200),
                sample(rep(c(1, -1),each=100), 200),
                sample(rep(c(1, -1),each=100), 200)),
              ncol=4)

Lam <- Lam + rnorm(dim(Lam)[1]*dim(Lam)[2], 0, 0.1)

simulations <- simRRmix(nsims=50, n=200, G=500, A=3, B=1, QC=0.05,     # Simulate Data
                        p=0.05, psi=0.5, sig20=.55, sig21=.23, 
                        trmt=0.5, mu=13, q=4, Lam=Lam)


trmt.id     <- which(simulations$Treatment.Groups == 1)
cont.id     <- which(simulations$Treatment.Groups == 0)
nsims       <- length(simulations$Simulated.Data)


#--------------------#
# Individual t-Tests #
#--------------------#


ttest.tstats <- as.list(rep(NA, nsims))

ttest.null   <- as.list(rep(NA, nsims))

ttest.TPR    <- as.list(rep(NA, nsims))
ttest.FPR    <- as.list(rep(NA, nsims))

ttest.PPV    <- as.list(rep(NA, nsims))
ttest.FNR    <- as.list(rep(NA, nsims))

ttest.FDR    <- as.list(rep(NA, nsims))
ttest.PWR    <- as.list(rep(NA, nsims))


i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)
  tstat <- rep(NA, G)
  nullp <- rep(NA, G)
  
  
  for (j in seq(G)){
    
    tstat[j] <- t.test(set[trmt.id, j], set[cont.id, j])$statistic
    
    nullp[j] <- 1 - t.test(set[trmt.id, j], set[cont.id, j])$p.value
    
  }
  
  ttest.tstats[[i]] <- tstat
  
  ttest.null[[i]]   <- nullp
  
  i <- i + 1  
  
}

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds                   

i <- 1

for (tstats in ttest.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits)) 
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 

  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))

  
  j <- 1                              
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)                                         # Truth = +
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))                       # Test  = +
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR

    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR

    
    j <- j + 1
    
  }
  
  ttest.TPR[[i]] <- TPR.vect
  ttest.FPR[[i]] <- FPR.vect
  
  ttest.PPV[[i]] <- PPV.vect
  ttest.FNR[[i]] <- FNR.vect

  ttest.FDR[[i]] <- FDR.vect
  ttest.PWR[[i]] <- PWR.vect

  
  i <- i + 1
  
}


## ttest Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="Independent t-test ROC Curves",
     asp=1)

for (i in seq(length(ttest.TPR))){
  
  lines(ttest.FPR[[i]], ttest.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "t-test"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


ttest.TPR.avmat <- ttest.TPR[[1]]

for (i in (2:length(ttest.TPR))){
  
  ttest.TPR.avmat <- cbind(ttest.TPR.avmat, ttest.TPR[[i]])
  
}

ttest.TPR.avg <- rowMeans(ttest.TPR.avmat, na.rm=T)


## FPR Averages


ttest.FPR.avmat <- ttest.FPR[[1]]

for (i in (2:length(ttest.FPR))){
  
  ttest.FPR.avmat <- cbind(ttest.FPR.avmat, ttest.FPR[[i]])
  
}

ttest.FPR.avg <- rowMeans(ttest.FPR.avmat, na.rm=T)


## PPV Averages


ttest.PPV.avmat <- ttest.PPV[[1]]

for (i in (2:length(ttest.PPV))){
  
  ttest.PPV.avmat <- cbind(ttest.PPV.avmat, ttest.PPV[[i]])
  
}

ttest.PPV.avg <- rowMeans(ttest.PPV.avmat, na.rm=T)


## FNR Averages


ttest.FNR.avmat <- ttest.FNR[[1]]

for (i in (2:length(ttest.FNR))){
  
  ttest.FNR.avmat <- cbind(ttest.FNR.avmat, ttest.FNR[[i]])
  
}

ttest.FNR.avg <- rowMeans(ttest.FNR.avmat, na.rm=T)

## null Averages (By Rank)


ttest.null.avmat <- ttest.null[[1]][order(ttest.null[[1]], decreasing=T)]

for (i in (2:length(ttest.null))){
  
  ttest.null.avmat <- cbind(ttest.null.avmat, ttest.null[[i]][order(ttest.null[[i]], decreasing=T)])
  
}

ttest.null.avg <- rowMeans(ttest.null.avmat, na.rm=T)


#-------#
# LIMMA #
#-------#


library(limma)

limma.tstats <- as.list(rep(NA, nsims))

limma.null   <- as.list(rep(NA, nsims))

limma.TPR    <- as.list(rep(NA, nsims))
limma.FPR    <- as.list(rep(NA, nsims))

limma.PPV    <- as.list(rep(NA, nsims))
limma.FNR    <- as.list(rep(NA, nsims))

limma.FDR    <- as.list(rep(NA, nsims))
limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  limma.tstats[[i]] <- ebayes$t[,2]
  
  limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 

  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))

  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR

    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR

    
    j <- j + 1
    
  }
  
  limma.TPR[[i]] <- TPR.vect
  limma.FPR[[i]] <- FPR.vect
  
  limma.PPV[[i]] <- PPV.vect
  limma.FNR[[i]] <- FNR.vect

  limma.FDR[[i]] <- FDR.vect
  limma.PWR[[i]] <- PWR.vect

  
  i <- i + 1
  
}


## LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="LIMMA ROC Curves", asp=1)

for (i in seq(length(limma.TPR))){
  
  lines(limma.FPR[[i]], limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


limma.TPR.avmat <- limma.TPR[[1]]

for (i in (2:length(limma.TPR))){
  
  limma.TPR.avmat <- cbind(limma.TPR.avmat, limma.TPR[[i]])
  
}

limma.TPR.avg <- rowMeans(limma.TPR.avmat, na.rm=T)


## FPR Averages


limma.FPR.avmat <- limma.FPR[[1]]

for (i in (2:length(limma.FPR))){
  
  limma.FPR.avmat <- cbind(limma.FPR.avmat, limma.FPR[[i]])
  
}

limma.FPR.avg <- rowMeans(limma.FPR.avmat, na.rm=T)

## PPV Averages


limma.PPV.avmat <- limma.PPV[[1]]

for (i in (2:length(limma.PPV))){
  
  limma.PPV.avmat <- cbind(limma.PPV.avmat, limma.PPV[[i]])
  
}

limma.PPV.avg <- rowMeans(limma.PPV.avmat, na.rm=T)


## FNR Averages


limma.FNR.avmat <- limma.FNR[[1]]

for (i in (2:length(limma.FNR))){
  
  limma.FNR.avmat <- cbind(limma.FNR.avmat, limma.FNR[[i]])
  
}

limma.FNR.avg <- rowMeans(limma.FNR.avmat, na.rm=T)


## null Averages (By Rank)


limma.null.avmat <- limma.null[[1]][order(limma.null[[1]], decreasing=T)]

for (i in (2:length(limma.null))){
  
  limma.null.avmat <- cbind(limma.null.avmat, limma.null[[i]][order(limma.null[[i]], decreasing=T)])
  
}

limma.null.avg <- rowMeans(limma.null.avmat, na.rm=T)


#-------#
# RRmix #
#-------#


RRmix.post <- as.list(rep(NA, nsims))
RRmix.TPR  <- as.list(rep(NA, nsims))
RRmix.FPR  <- as.list(rep(NA, nsims))

RRmix.PPV    <- as.list(rep(NA, nsims))
RRmix.FNR    <- as.list(rep(NA, nsims))

RRmix.FDR    <- as.list(rep(NA, nsims))
RRmix.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                        # Set Treatment Groups

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

i <- 1

for (set in simulations$Simulated.Data){
  
  G     <- ncol(set)                                             # Number of Metabolites
  n     <- nrow(set)                                             # Number of Observations
  Xc    <- matrix(nrow=0, ncol=0)                                # Covariate Matrix
  mu.0  <- 1/G * as.matrix(set) %*% rep(1,G)                     # Initialize mu
  eta.0 <- matrix(0, 2+ncol(Xc), G)                              # Initialize eta
  
  betac.0 <- matrix(nrow=0, ncol=0)                              # Initialize beta_c
  sig20.0 <- 1                                                   # Initialize sig^2_0
  sig21.0 <- 0.1                                                 # Initialize sig^2_1
  
  result <- runHEFTmix(G.in=G,                                   # Run RRmix Model
                       n.in=n, 
                       Xc.in=Xc, 
                       Y.in=as.matrix(set), 
                       SNP.in=trmt.ind,
                       mu.0=mu.0, 
                       betac.0=betac.0, 
                       sig20.0=sig20.0, 
                       sig21.0=sig21.0, 
                       p.0=0.05, 
                       er_tol.in=10^(-3),   
                       q.in=4)  
  
  RRmix.post[[i]] <- result[['b_g']]
  
  i <- i + 1  
  
}  

post.probs <- seq(0.0, 1.0, by=0.0001)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (posts in RRmix.post){
  
  TPR.vect <- rep(NA, length(post.probs))
  FPR.vect <- rep(NA, length(post.probs))
  
  PPV.vect <- rep(NA, length(post.probs))
  FNR.vect <- rep(NA, length(post.probs)) 

  FDR.vect <- rep(NA, length(post.probs)) 
  PWR.vect <- rep(NA, length(post.probs))

  
  j <- 1
  
  for (post.prob in post.probs){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(posts > post.prob)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR

    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR

    
    j <- j + 1
    
  }
  
  RRmix.TPR[[i]] <- TPR.vect
  RRmix.FPR[[i]] <- FPR.vect
  
  RRmix.PPV[[i]] <- PPV.vect
  RRmix.FNR[[i]] <- FNR.vect

  RRmix.FDR[[i]] <- FDR.vect
  RRmix.PWR[[i]] <- PWR.vect

  
  i <- i + 1
  
}


## RRmix Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RRmix ROC Curves", asp=1)

for (i in seq(length(RRmix.TPR))){
  
  lines(RRmix.FPR[[i]], RRmix.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RRmix"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RRmix.TPR.avmat <- RRmix.TPR[[1]]

for (i in (2:length(RRmix.TPR))){
  
  RRmix.TPR.avmat <- cbind(RRmix.TPR.avmat, RRmix.TPR[[i]])
  
}

RRmix.TPR.avg <- rowMeans(RRmix.TPR.avmat, na.rm=T)


## FPR Averages


RRmix.FPR.avmat <- RRmix.FPR[[1]]

for (i in (2:length(RRmix.FPR))){
  
  RRmix.FPR.avmat <- cbind(RRmix.FPR.avmat, RRmix.FPR[[i]])
  
}

RRmix.FPR.avg <- rowMeans(RRmix.FPR.avmat, na.rm=T)


## PPV Averages


RRmix.PPV.avmat <- RRmix.PPV[[1]]

for (i in (2:length(RRmix.PPV))){
  
  RRmix.PPV.avmat <- cbind(RRmix.PPV.avmat, RRmix.PPV[[i]])
  
}

RRmix.PPV.avg <- rowMeans(RRmix.PPV.avmat, na.rm=T)


## FNR Averages


RRmix.FNR.avmat <- RRmix.FNR[[1]]

for (i in (2:length(RRmix.FNR))){
  
  RRmix.FNR.avmat <- cbind(RRmix.FNR.avmat, RRmix.FNR[[i]])
  
}

RRmix.FNR.avg <- rowMeans(RRmix.FNR.avmat, na.rm=T)


## FDR Averages

RRmix.FDR.avmat <- RRmix.FDR[[1]]

for (i in (2:length(RRmix.FDR))){
  
  RRmix.FDR.avmat <- cbind(RRmix.FDR.avmat, RRmix.FDR[[i]])
  
}

RRmix.FDR.avg <- rowMeans(RRmix.FDR.avmat, na.rm=T)


## PWR Averages

RRmix.PWR.avmat <- RRmix.PWR[[1]]

for (i in (2:length(RRmix.PWR))){
  
  RRmix.PWR.avmat <- cbind(RRmix.PWR.avmat, RRmix.PWR[[i]])
  
}

RRmix.PWR.avg <- rowMeans(RRmix.PWR.avmat, na.rm=T)


## post Averages (By Rank, Not Index)

RRmix.post.avmat <- RRmix.post[[1]][order(RRmix.post[[1]])]

for (i in (2:length(RRmix.post))){
  
  RRmix.post.avmat <- cbind(RRmix.post.avmat, RRmix.post[[i]][order(RRmix.post[[i]])])
  
}

RRmix.post.avg <- rowMeans(RRmix.post.avmat, na.rm=T)

RRmix.null.avg <- 1 - RRmix.post.avg 


#----------------#
# FAMT - NBF Set #
#----------------#

library(FAMT)

FAMT.NBF.Fstats <- as.list(rep(NA, nsims))

FAMT.NBF.null   <- as.list(rep(NA, nsims))

FAMT.NBF.TPR    <- as.list(rep(NA, nsims))
FAMT.NBF.FPR    <- as.list(rep(NA, nsims))

FAMT.NBF.PPV    <- as.list(rep(NA, nsims))
FAMT.NBF.FNR    <- as.list(rep(NA, nsims))

FAMT.NBF.FDR    <- as.list(rep(NA, nsims))
FAMT.NBF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  expr.FAMT.NBF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.NBF) <- (1:ncol(expr.FAMT.NBF))
  
  cov.FAMT.NBF  <- data.frame(id    = colnames(expr.FAMT.NBF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.NBF <- as.FAMTdata(expression = expr.FAMT.NBF,     # Make Data Structure
                               covariates = cov.FAMT.NBF, 
                               idcovar    = 1)
  
  fit.FAMT.NBF  <- modelFAMT(data.FAMT.NBF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2, 
                             nbf  = 4)
  
  FAMT.NBF.Fstats[[i]] <- fit.FAMT.NBF$adjtest
  
  FAMT.NBF.null[[i]]   <- 1 - fit.FAMT.NBF$adjpval
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.NBF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 

  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))

  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR

    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR

    
    j <- j + 1
    
  }
  
  FAMT.NBF.TPR[[i]] <- TPR.vect
  FAMT.NBF.FPR[[i]] <- FPR.vect
  
  FAMT.NBF.PPV[[i]] <- PPV.vect
  FAMT.NBF.FNR[[i]] <- FNR.vect

  FAMT.NBF.FDR[[i]] <- FDR.vect
  FAMT.NBF.PWR[[i]] <- PWR.vect

  
  i <- i + 1
  
}


## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT NBF ROC Curves", asp=1)

for (i in seq(length(FAMT.NBF.TPR))){
  
  lines(FAMT.NBF.FPR[[i]], FAMT.NBF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.NBF.TPR.avmat <- FAMT.NBF.TPR[[1]]

for (i in (2:length(FAMT.NBF.TPR))){
  
  FAMT.NBF.TPR.avmat <- cbind(FAMT.NBF.TPR.avmat, FAMT.NBF.TPR[[i]])
  
}

FAMT.NBF.TPR.avg <- rowMeans(FAMT.NBF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.NBF.FPR.avmat <- FAMT.NBF.FPR[[1]]

for (i in (2:length(FAMT.NBF.FPR))){
  
  FAMT.NBF.FPR.avmat <- cbind(FAMT.NBF.FPR.avmat, FAMT.NBF.FPR[[i]])
  
}

FAMT.NBF.FPR.avg <- rowMeans(FAMT.NBF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.NBF.PPV.avmat <- FAMT.NBF.PPV[[1]]

for (i in (2:length(FAMT.NBF.PPV))){
  
  FAMT.NBF.PPV.avmat <- cbind(FAMT.NBF.PPV.avmat, FAMT.NBF.PPV[[i]])
  
}

FAMT.NBF.PPV.avg <- rowMeans(FAMT.NBF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.NBF.FNR.avmat <- FAMT.NBF.FNR[[1]]

for (i in (2:length(FAMT.NBF.FNR))){
  
  FAMT.NBF.FNR.avmat <- cbind(FAMT.NBF.FNR.avmat, FAMT.NBF.FNR[[i]])
  
}

FAMT.NBF.FNR.avg <- rowMeans(FAMT.NBF.FNR.avmat, na.rm=T)

## null Averages (By Rank)


FAMT.NBF.null.avmat <- FAMT.NBF.null[[1]][order(FAMT.NBF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.NBF.null))){
  
  FAMT.NBF.null.avmat <- cbind(FAMT.NBF.null.avmat, FAMT.NBF.null[[i]][order(FAMT.NBF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.NBF.null.avg <- rowMeans(FAMT.NBF.null.avmat, na.rm=T)


#----------------#
# FAMT - Default #
#----------------#


library(FAMT)

FAMT.DEF.Fstats <- as.list(rep(NA, nsims))

FAMT.DEF.null   <- as.list(rep(NA, nsims))

FAMT.DEF.TPR    <- as.list(rep(NA, nsims))
FAMT.DEF.FPR    <- as.list(rep(NA, nsims))

FAMT.DEF.PPV    <- as.list(rep(NA, nsims))
FAMT.DEF.FNR    <- as.list(rep(NA, nsims))

FAMT.DEF.FDR    <- as.list(rep(NA, nsims))
FAMT.DEF.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                      # Set Treatment Groups

i <- 1

for (set in simulations$Simulated.Data){
  
  tryCatch({
  
  expr.FAMT.DEF <- t(set)                                      # Set Data Matrix
  colnames(expr.FAMT.DEF) <- (1:ncol(expr.FAMT.DEF))
  
  cov.FAMT.DEF  <- data.frame(id    = colnames(expr.FAMT.DEF), # Set Covariates Matrix
                              trmt  = as.factor(trmt.ind)) 
  
  data.FAMT.DEF <- as.FAMTdata(expression = expr.FAMT.DEF,     # Make Data Structure
                               covariates = cov.FAMT.DEF, 
                               idcovar    = 1)
  
  nbf.FAMT <- nbfactors(data.FAMT.DEF,                         # Determine Number of Factors
                        x=2, 
                        test=2,
                        maxnbfactors = 8)$optimalnbfactors
  
  fit.FAMT.DEF  <- modelFAMT(data.FAMT.DEF,                    # Fit FAMT Model
                             x    = 2, 
                             test = 2,
                             nbf  = nbf.FAMT)
  
  FAMT.DEF.Fstats[[i]] <- fit.FAMT.DEF$adjtest
  
  FAMT.DEF.null[[i]]   <- 1 - fit.FAMT.DEF$adjpval
  
  }, error=function(e){"FAMT FAILED TO CONVERGE ON A SOLUTION"})
  
  i <- i + 1  
  
}  

Fcrits <- seq(0, 1000.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (Fstats in FAMT.DEF.Fstats){
  
  TPR.vect <- rep(NA, length(Fcrits))
  FPR.vect <- rep(NA, length(Fcrits))
  
  PPV.vect <- rep(NA, length(Fcrits))
  FNR.vect <- rep(NA, length(Fcrits)) 

  FDR.vect <- rep(NA, length(Fcrits)) 
  PWR.vect <- rep(NA, length(Fcrits))

  
  j <- 1
  
  for (Fcrit in Fcrits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which(Fstats > Fcrit)
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR

    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR

    
    j <- j + 1
    
  }
  
  FAMT.DEF.TPR[[i]] <- TPR.vect
  FAMT.DEF.FPR[[i]] <- FPR.vect
  
  FAMT.DEF.PPV[[i]] <- PPV.vect
  FAMT.DEF.FNR[[i]] <- FNR.vect

  FAMT.DEF.FDR[[i]] <- FDR.vect
  FAMT.DEF.PWR[[i]] <- PWR.vect

  
  i <- i + 1
  
}

FAMT.DEF.TPR <- lapply(FAMT.DEF.TPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.TPR[[1]])))])
FAMT.DEF.TPR <- Filter(length, FAMT.DEF.TPR)


FAMT.DEF.FPR <- lapply(FAMT.DEF.FPR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FPR[[1]])))])
FAMT.DEF.FPR <- Filter(length, FAMT.DEF.FPR)

FAMT.DEF.PPV <- lapply(FAMT.DEF.PPV, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PPV[[1]])))])
FAMT.DEF.PPV <- Filter(length, FAMT.DEF.PPV)

FAMT.DEF.FNR <- lapply(FAMT.DEF.FNR, 
                       function(x) x[!identical(x, rep(1, length(FAMT.DEF.FNR[[1]])))])
FAMT.DEF.FNR <- Filter(length, FAMT.DEF.FNR)

FAMT.DEF.FDR <- lapply(FAMT.DEF.FDR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.FDR[[1]])))])
FAMT.DEF.FDR <- Filter(length, FAMT.DEF.FDR)

FAMT.DEF.PWR <- lapply(FAMT.DEF.PWR, 
                       function(x) x[!identical(x, rep(0, length(FAMT.DEF.PWR[[1]])))])
FAMT.DEF.PWR <- Filter(length, FAMT.DEF.PWR)

## FAMT Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="FAMT DEF ROC Curves", asp=1)

for (i in seq(length(FAMT.DEF.TPR))){
  
  lines(FAMT.DEF.FPR[[i]], FAMT.DEF.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "FAMT"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


FAMT.DEF.TPR.avmat <- FAMT.DEF.TPR[[1]]

for (i in (2:length(FAMT.DEF.TPR))){
  
  FAMT.DEF.TPR.avmat <- cbind(FAMT.DEF.TPR.avmat, FAMT.DEF.TPR[[i]])
  
}

FAMT.DEF.TPR.avg <- rowMeans(FAMT.DEF.TPR.avmat, na.rm=T)


## FPR Averages


FAMT.DEF.FPR.avmat <- FAMT.DEF.FPR[[1]]

for (i in (2:length(FAMT.DEF.FPR))){
  
  FAMT.DEF.FPR.avmat <- cbind(FAMT.DEF.FPR.avmat, FAMT.DEF.FPR[[i]])
  
}

FAMT.DEF.FPR.avg <- rowMeans(FAMT.DEF.FPR.avmat, na.rm=T)

## PPV Averages


FAMT.DEF.PPV.avmat <- FAMT.DEF.PPV[[1]]

for (i in (2:length(FAMT.DEF.PPV))){
  
  FAMT.DEF.PPV.avmat <- cbind(FAMT.DEF.PPV.avmat, FAMT.DEF.PPV[[i]])
  
}

FAMT.DEF.PPV.avg <- rowMeans(FAMT.DEF.PPV.avmat, na.rm=T)


## FNR Averages


FAMT.DEF.FNR.avmat <- FAMT.DEF.FNR[[1]]

for (i in (2:length(FAMT.DEF.FNR))){
  
  FAMT.DEF.FNR.avmat <- cbind(FAMT.DEF.FNR.avmat, FAMT.DEF.FNR[[i]])
  
}

FAMT.DEF.FNR.avg <- rowMeans(FAMT.DEF.FNR.avmat, na.rm=T)



## null Averages (By Rank)


FAMT.DEF.null.avmat <- FAMT.DEF.null[[1]][order(FAMT.DEF.null[[1]], decreasing=T)]

for (i in (2:length(FAMT.DEF.null))){
  
  FAMT.DEF.null.avmat <- cbind(FAMT.DEF.null.avmat, FAMT.DEF.null[[i]][order(FAMT.DEF.null[[i]],
                                                                             decreasing=T)])
  
}

FAMT.DEF.null.avg <- rowMeans(FAMT.DEF.null.avmat, na.rm=T)


#--------------------------#
# UNSUPERVISED SVA + LIMMA #
#--------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])


UNSUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  batch <- sva(set, mod, mod0)
  
  UNSUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}



UNSUPSVA_limma.tstats <- as.list(rep(NA, nsims))

UNSUPSVA_limma.null   <- as.list(rep(NA, nsims))

UNSUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

UNSUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
UNSUPSVA_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups


for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- UNSUPSVA.mods[[i]]                                     # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  UNSUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  UNSUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in UNSUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  UNSUPSVA_limma.TPR[[i]] <- TPR.vect
  UNSUPSVA_limma.FPR[[i]] <- FPR.vect
  
  UNSUPSVA_limma.PPV[[i]] <- PPV.vect
  UNSUPSVA_limma.FNR[[i]] <- FNR.vect
  
  UNSUPSVA_limma.FDR[[i]] <- FDR.vect
  UNSUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## UNSUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="UNSUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(UNSUPSVA_limma.TPR))){
  
  lines(UNSUPSVA_limma.FPR[[i]], UNSUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "UNSUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


UNSUPSVA_limma.TPR.avmat <- UNSUPSVA_limma.TPR[[1]]

for (i in (2:length(UNSUPSVA_limma.TPR))){
  
  UNSUPSVA_limma.TPR.avmat <- cbind(UNSUPSVA_limma.TPR.avmat, UNSUPSVA_limma.TPR[[i]])
  
}

UNSUPSVA_limma.TPR.avg <- rowMeans(UNSUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


UNSUPSVA_limma.FPR.avmat <- UNSUPSVA_limma.FPR[[1]]

for (i in (2:length(UNSUPSVA_limma.FPR))){
  
  UNSUPSVA_limma.FPR.avmat <- cbind(UNSUPSVA_limma.FPR.avmat, UNSUPSVA_limma.FPR[[i]])
  
}

UNSUPSVA_limma.FPR.avg <- rowMeans(UNSUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


UNSUPSVA_limma.PPV.avmat <- UNSUPSVA_limma.PPV[[1]]

for (i in (2:length(UNSUPSVA_limma.PPV))){
  
  UNSUPSVA_limma.PPV.avmat <- cbind(UNSUPSVA_limma.PPV.avmat, UNSUPSVA_limma.PPV[[i]])
  
}

UNSUPSVA_limma.PPV.avg <- rowMeans(UNSUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


UNSUPSVA_limma.FNR.avmat <- UNSUPSVA_limma.FNR[[1]]

for (i in (2:length(UNSUPSVA_limma.FNR))){
  
  UNSUPSVA_limma.FNR.avmat <- cbind(UNSUPSVA_limma.FNR.avmat, UNSUPSVA_limma.FNR[[i]])
  
}

UNSUPSVA_limma.FNR.avg <- rowMeans(UNSUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


UNSUPSVA_limma.null.avmat <- UNSUPSVA_limma.null[[1]][order(UNSUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(UNSUPSVA_limma.null))){
  
  UNSUPSVA_limma.null.avmat <- cbind(UNSUPSVA_limma.null.avmat, 
                                     UNSUPSVA_limma.null[[i]][order(UNSUPSVA_limma.null[[i]], 
                                                                    decreasing=T)])
  
}

UNSUPSVA_limma.null.avg <- rowMeans(UNSUPSVA_limma.null.avmat, na.rm=T)


#------------------------#
# SUPERVISED SVA + LIMMA #
#------------------------#


library(sva)
library(limma)


mod <- model.matrix(~simulations$Treatment.Groups)
mod0 <- cbind(mod[,1])

SUPSVA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- t(as.matrix(simulations$Simulated.Data[[i]]))
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, nrow(set))
  
  controls[QCs] <- 1
  
  batch <- sva(set, mod, mod0, controls=controls, method="supervised")
  
  SUPSVA.mods[[i]] = cbind(mod, batch$sv)
  
}


SUPSVA_limma.tstats <- as.list(rep(NA, nsims))

SUPSVA_limma.null   <- as.list(rep(NA, nsims))

SUPSVA_limma.TPR    <- as.list(rep(NA, nsims))
SUPSVA_limma.FPR    <- as.list(rep(NA, nsims))

SUPSVA_limma.PPV    <- as.list(rep(NA, nsims))
SUPSVA_limma.FNR    <- as.list(rep(NA, nsims))

SUPSVA_limma.FDR    <- as.list(rep(NA, nsims))
SUPSVA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- SUPSVA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  SUPSVA_limma.tstats[[i]] <- ebayes$t[,2]
  
  SUPSVA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in SUPSVA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  SUPSVA_limma.TPR[[i]] <- TPR.vect
  SUPSVA_limma.FPR[[i]] <- FPR.vect
  
  SUPSVA_limma.PPV[[i]] <- PPV.vect
  SUPSVA_limma.FNR[[i]] <- FNR.vect
  
  SUPSVA_limma.FDR[[i]] <- FDR.vect
  SUPSVA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## SUPSVA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="SUPSVA + LIMMA ROC Curves", asp=1)

for (i in seq(length(SUPSVA_limma.TPR))){
  
  lines(SUPSVA_limma.FPR[[i]], SUPSVA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "SUPSVA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


SUPSVA_limma.TPR.avmat <- SUPSVA_limma.TPR[[1]]

for (i in (2:length(SUPSVA_limma.TPR))){
  
  SUPSVA_limma.TPR.avmat <- cbind(SUPSVA_limma.TPR.avmat, SUPSVA_limma.TPR[[i]])
  
}

SUPSVA_limma.TPR.avg <- rowMeans(SUPSVA_limma.TPR.avmat, na.rm=T)


## FPR Averages


SUPSVA_limma.FPR.avmat <- SUPSVA_limma.FPR[[1]]

for (i in (2:length(SUPSVA_limma.FPR))){
  
  SUPSVA_limma.FPR.avmat <- cbind(SUPSVA_limma.FPR.avmat, SUPSVA_limma.FPR[[i]])
  
}

SUPSVA_limma.FPR.avg <- rowMeans(SUPSVA_limma.FPR.avmat, na.rm=T)

## PPV Averages


SUPSVA_limma.PPV.avmat <- SUPSVA_limma.PPV[[1]]

for (i in (2:length(SUPSVA_limma.PPV))){
  
  SUPSVA_limma.PPV.avmat <- cbind(SUPSVA_limma.PPV.avmat, SUPSVA_limma.PPV[[i]])
  
}

SUPSVA_limma.PPV.avg <- rowMeans(SUPSVA_limma.PPV.avmat, na.rm=T)


## FNR Averages


SUPSVA_limma.FNR.avmat <- SUPSVA_limma.FNR[[1]]

for (i in (2:length(SUPSVA_limma.FNR))){
  
  SUPSVA_limma.FNR.avmat <- cbind(SUPSVA_limma.FNR.avmat, SUPSVA_limma.FNR[[i]])
  
}

SUPSVA_limma.FNR.avg <- rowMeans(SUPSVA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


SUPSVA_limma.null.avmat <- SUPSVA_limma.null[[1]][order(SUPSVA_limma.null[[1]], decreasing=T)]

for (i in (2:length(SUPSVA_limma.null))){
  
  SUPSVA_limma.null.avmat <- cbind(SUPSVA_limma.null.avmat, 
                                   SUPSVA_limma.null[[i]][order(SUPSVA_limma.null[[i]],
                                                                decreasing=T)])
  
}

SUPSVA_limma.null.avg <- rowMeans(SUPSVA_limma.null.avmat, na.rm=T)


#-------------#
# PCA + LIMMA #
#-------------#


library(limma)


PCA.mods <- as.list(rep(NA, length(simulations$Simulated.Data)))


for (i in 1:length(simulations$Simulated.Data)){
  
  set <- as.matrix(simulations$Simulated.Data[[i]])

  batches <- svd(t(set) - rowMeans(t(set)))$v[,1]
  
  PCA.mods[[i]] <- model.matrix(~simulations$Treatment.Groups+batches)
  
}


PCA_limma.tstats <- as.list(rep(NA, nsims))

PCA_limma.null   <- as.list(rep(NA, nsims))

PCA_limma.TPR    <- as.list(rep(NA, nsims))
PCA_limma.FPR    <- as.list(rep(NA, nsims))

PCA_limma.PPV    <- as.list(rep(NA, nsims))
PCA_limma.FNR    <- as.list(rep(NA, nsims))

PCA_limma.FDR    <- as.list(rep(NA, nsims))
PCA_limma.PWR    <- as.list(rep(NA, nsims))

for (i in 1:length(simulations$Simulated.Data)){
  
  set    <- simulations$Simulated.Data[[i]]                        # Extract Dataset
  mod    <- PCA.mods[[i]]                                       # Extract Design Matrix
  fit    <- lmFit(t(set), mod)                                     # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  PCA_limma.tstats[[i]] <- ebayes$t[,2]
  
  PCA_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in PCA_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR
    
    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR
    
    j <- j + 1
    
  }
  
  PCA_limma.TPR[[i]] <- TPR.vect
  PCA_limma.FPR[[i]] <- FPR.vect
  
  PCA_limma.PPV[[i]] <- PPV.vect
  PCA_limma.FNR[[i]] <- FNR.vect
  
  PCA_limma.FDR[[i]] <- FDR.vect
  PCA_limma.PWR[[i]] <- PWR.vect
  
  
  i <- i + 1
  
}


## PCA + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="PCA + LIMMA ROC Curves", asp=1)

for (i in seq(length(PCA_limma.TPR))){
  
  lines(PCA_limma.FPR[[i]], PCA_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "PCA + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


PCA_limma.TPR.avmat <- PCA_limma.TPR[[1]]

for (i in (2:length(PCA_limma.TPR))){
  
  PCA_limma.TPR.avmat <- cbind(PCA_limma.TPR.avmat, PCA_limma.TPR[[i]])
  
}

PCA_limma.TPR.avg <- rowMeans(PCA_limma.TPR.avmat, na.rm=T)


## FPR Averages


PCA_limma.FPR.avmat <- PCA_limma.FPR[[1]]

for (i in (2:length(PCA_limma.FPR))){
  
  PCA_limma.FPR.avmat <- cbind(PCA_limma.FPR.avmat, PCA_limma.FPR[[i]])
  
}

PCA_limma.FPR.avg <- rowMeans(PCA_limma.FPR.avmat, na.rm=T)

## PPV Averages


PCA_limma.PPV.avmat <- PCA_limma.PPV[[1]]

for (i in (2:length(PCA_limma.PPV))){
  
  PCA_limma.PPV.avmat <- cbind(PCA_limma.PPV.avmat, PCA_limma.PPV[[i]])
  
}

PCA_limma.PPV.avg <- rowMeans(PCA_limma.PPV.avmat, na.rm=T)


## FNR Averages


PCA_limma.FNR.avmat <- PCA_limma.FNR[[1]]

for (i in (2:length(PCA_limma.FNR))){
  
  PCA_limma.FNR.avmat <- cbind(PCA_limma.FNR.avmat, PCA_limma.FNR[[i]])
  
}

PCA_limma.FNR.avg <- rowMeans(PCA_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


PCA_limma.null.avmat <- PCA_limma.null[[1]][order(PCA_limma.null[[1]], decreasing=T)]

for (i in (2:length(PCA_limma.null))){
  
  PCA_limma.null.avmat <- cbind(PCA_limma.null.avmat, PCA_limma.null[[i]][order(PCA_limma.null[[i]],
                                                                                decreasing=T)])
  
}

PCA_limma.null.avg <- rowMeans(PCA_limma.null.avmat, na.rm=T)


#------------------------------------------#
# RUV WITH NEGATIVE CONTROLS KNOWN + LIMMA #
#------------------------------------------#


library(MetNorm)
library(limma)


RUV.sets <- as.list(rep(NA, length(simulations$Simulated.Data)))

for (i in 1:length(simulations$Simulated.Data)){
  
  QCs <- simulations$Quality.Controls[[i]]
  
  controls <- rep(0, ncol(set))
  
  controls[QCs] <- 1
  
  controls <- as.logical(controls)
  
  set <- as.matrix(simulations$Simulated.Data[[i]])
  
  RUV.sets[[i]] <- NormalizeRUVRand(set, k=2, ctl=controls)$newY
  
}
  

RUV_limma.tstats <- as.list(rep(NA, nsims))

RUV_limma.null   <- as.list(rep(NA, nsims))

RUV_limma.TPR    <- as.list(rep(NA, nsims))
RUV_limma.FPR    <- as.list(rep(NA, nsims))

RUV_limma.PPV    <- as.list(rep(NA, nsims))
RUV_limma.FNR    <- as.list(rep(NA, nsims))

RUV_limma.FDR    <- as.list(rep(NA, nsims))
RUV_limma.PWR    <- as.list(rep(NA, nsims))


trmt.ind  <- simulations$Treatment.Groups                          # Set Treatment Groups

i <- 1

for (set in RUV.sets){
  
  design <- model.matrix(~factor(trmt.ind))                        # Create Design Matrix
  fit    <- lmFit(t(set),design)                                   # Fit Model
  ebayes <- eBayes(fit)                                            # Get Statistics
  
  RUV_limma.tstats[[i]] <- ebayes$t[,2]
  
  RUV_limma.null[[i]]   <- 1 - exp(ebayes$lods[,2]) / ( 1 + exp(ebayes$lods[,2]))
  
  i <- i + 1  
  
}  

t.crits <- seq(0.0, 50.0, by=0.1)

diff.genes <- simulations$Differential.Compounds

i <- 1

for (tstats in RUV_limma.tstats){
  
  TPR.vect <- rep(NA, length(t.crits))
  FPR.vect <- rep(NA, length(t.crits))
  
  PPV.vect <- rep(NA, length(t.crits))
  FNR.vect <- rep(NA, length(t.crits)) 
  
  FDR.vect <- rep(NA, length(t.crits)) 
  PWR.vect <- rep(NA, length(t.crits))
  
  j <- 1
  
  for (t.crit in t.crits){
    
    diff <- which(diff.genes[[i]] == 1)
    pos  <- which((tstats > t.crit) | (tstats < -t.crit))
    
    TP   <- length(intersect(pos, diff))                                        # TP
    TN   <- length(intersect(setdiff((1:G), diff), setdiff((1:G), pos)))        # TN
    FP   <- length(setdiff(pos, diff))                                          # FP
    FN   <- length(setdiff(diff, pos))                                          # FN
    
    TPR  <- TP/(TP+FN)
    FPR  <- FP/(FP+TN)
    
    TPR.vect[j] <- TPR
    FPR.vect[j] <- FPR
    
    PPV  <- TP/(TP+FP)
    
    PPV.vect[j] <- PPV
    
    FNR <- FN/(TP+FN)
    
    FNR.vect[j] <- FNR

    FDR  <- FP/(TP+FP)
    PWR  <- TP/(TP+FN)
    
    FDR.vect[j] <- FDR
    PWR.vect[j] <- PWR

    j <- j + 1
    
  }
  
  RUV_limma.TPR[[i]] <- TPR.vect
  RUV_limma.FPR[[i]] <- FPR.vect
  
  RUV_limma.PPV[[i]] <- PPV.vect
  RUV_limma.FNR[[i]] <- FNR.vect
 
  RUV_limma.FDR[[i]] <- FDR.vect
  RUV_limma.PWR[[i]] <- PWR.vect

  
  i <- i + 1
  
}


## RUV + LIMMA Simulated ROC Curves

par(pty="s")
plot(c(0,1), c(0,1),
     type="l",
     lty=2,
     col="gray",
     xlab="False Positive Rate (1-Specificity)",
     ylab="True Positive Rate (Sensitivity)",
     main="RUV + LIMMA ROC Curves", asp=1)

for (i in seq(length(RUV_limma.TPR))){
  
  lines(RUV_limma.FPR[[i]], RUV_limma.TPR[[i]])
  
}

legend("bottomright", 
       legend = c("Random Guessing", "RUV + LIMMA"), 
       col    = c("gray",1),
       lty    = c(2,1), 
       lwd    = c(1,1),
       cex    = 0.6)


## TPR Averages


RUV_limma.TPR.avmat <- RUV_limma.TPR[[1]]

for (i in (2:length(RUV_limma.TPR))){
  
  RUV_limma.TPR.avmat <- cbind(RUV_limma.TPR.avmat, RUV_limma.TPR[[i]])
  
}

RUV_limma.TPR.avg <- rowMeans(RUV_limma.TPR.avmat, na.rm=T)


## FPR Averages


RUV_limma.FPR.avmat <- RUV_limma.FPR[[1]]

for (i in (2:length(RUV_limma.FPR))){
  
  RUV_limma.FPR.avmat <- cbind(RUV_limma.FPR.avmat, RUV_limma.FPR[[i]])
  
}

RUV_limma.FPR.avg <- rowMeans(RUV_limma.FPR.avmat, na.rm=T)

## PPV Averages


RUV_limma.PPV.avmat <- RUV_limma.PPV[[1]]

for (i in (2:length(RUV_limma.PPV))){
  
  RUV_limma.PPV.avmat <- cbind(RUV_limma.PPV.avmat, RUV_limma.PPV[[i]])
  
}

RUV_limma.PPV.avg <- rowMeans(RUV_limma.PPV.avmat, na.rm=T)


## FNR Averages


RUV_limma.FNR.avmat <- RUV_limma.FNR[[1]]

for (i in (2:length(RUV_limma.FNR))){
  
  RUV_limma.FNR.avmat <- cbind(RUV_limma.FNR.avmat, RUV_limma.FNR[[i]])
  
}

RUV_limma.FNR.avg <- rowMeans(RUV_limma.FNR.avmat, na.rm=T)

## null Averages (By Rank)


RUV_limma.null.avmat <- RUV_limma.null[[1]][order(RUV_limma.null[[1]], decreasing=T)]

for (i in (2:length(RUV_limma.null))){
  
  RUV_limma.null.avmat <- cbind(RUV_limma.null.avmat, RUV_limma.null[[i]][order(RUV_limma.null[[i]],
                                                                                decreasing=T)])
  
}

RUV_limma.null.avg <- rowMeans(RUV_limma.null.avmat, na.rm=T)


#-----#
# AUC #
#-----#


# t-Test

height = (ttest.TPR.avg[-1]+ttest.TPR.avg[-length(ttest.TPR.avg)])/2
width = -diff(ttest.FPR.avg)
ttest.auc <- sum(height*width)
ttest.auc

# LIMMA

height = (limma.TPR.avg[-1]+limma.TPR.avg[-length(limma.TPR.avg)])/2
width = -diff(limma.FPR.avg)
limma.auc <- sum(height*width)
limma.auc

# RRmix

height = (RRmix.TPR.avg[-1]+RRmix.TPR.avg[-length(RRmix.TPR.avg)])/2
width = -diff(RRmix.FPR.avg)
RRmix.auc <- sum(height*width)
RRmix.auc 

# FAMT NBF

height = (FAMT.NBF.TPR.avg[-1]+FAMT.NBF.TPR.avg[-length(FAMT.NBF.TPR.avg)])/2
width = -diff(FAMT.NBF.FPR.avg)
FAMT.NBF.auc <- sum(height*width)
FAMT.NBF.auc

# FAMT DEF

height = (FAMT.DEF.TPR.avg[-1]+FAMT.DEF.TPR.avg[-length(FAMT.DEF.TPR.avg)])/2
width = -diff(FAMT.DEF.FPR.avg)
FAMT.DEF.auc <- sum(height*width)
FAMT.DEF.auc

# UNSUPSVA + LIMMA

height = (UNSUPSVA_limma.TPR.avg[-1]+UNSUPSVA_limma.TPR.avg[-length(UNSUPSVA_limma.TPR.avg)])/2
width = -diff(UNSUPSVA_limma.FPR.avg)
UNSUPSVA_limma.auc <- sum(height*width)
UNSUPSVA_limma.auc

# SUPSVA + LIMMA

height = (SUPSVA_limma.TPR.avg[-1]+SUPSVA_limma.TPR.avg[-length(SUPSVA_limma.TPR.avg)])/2
width = -diff(SUPSVA_limma.FPR.avg)
SUPSVA_limma.auc <- sum(height*width)
SUPSVA_limma.auc

# PCA + LIMMA

height = (PCA_limma.TPR.avg[-1]+PCA_limma.TPR.avg[-length(PCA_limma.TPR.avg)])/2
width = -diff(PCA_limma.FPR.avg)
PCA_limma.auc <- sum(height*width)
PCA_limma.auc

# RUV + LIMMA

height = (RUV_limma.TPR.avg[-1]+RUV_limma.TPR.avg[-length(RUV_limma.TPR.avg)])/2
width = -diff(RUV_limma.FPR.avg)
RUV_limma.auc <- sum(height*width)
RUV_limma.auc

auc.vect <- c(ttest.auc, limma.auc, RRmix.auc, FAMT.NBF.auc, FAMT.DEF.auc, UNSUPSVA_limma.auc,
              SUPSVA_limma.auc, PCA_limma.auc, RUV_limma.auc)


#-------------------------------------#
# Average Model ROC Curve Comparisons #
#-------------------------------------#


dat.FPR <- rbind(as.matrix(ttest.FPR.avg, ncol=1), as.matrix(limma.FPR.avg, ncol=1), 
                 as.matrix(RRmix.FPR.avg, ncol=1), as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.FPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.FPR.avg, ncol=1), as.matrix(PCA_limma.FPR.avg, ncol=1), 
                 as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.TPR <- rbind(as.matrix(ttest.TPR.avg, ncol=1), as.matrix(limma.TPR.avg, ncol=1), 
                 as.matrix(RRmix.TPR.avg, ncol=1), as.matrix(FAMT.NBF.TPR.avg, ncol=1), 
                 as.matrix(FAMT.DEF.TPR.avg, ncol=1), as.matrix(UNSUPSVA_limma.TPR.avg, ncol=1), 
                 as.matrix(SUPSVA_limma.TPR.avg, ncol=1), as.matrix(PCA_limma.TPR.avg, ncol=1), 
                 as.matrix(RUV_limma.TPR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       TPR = dat.TPR)


dat.refline <- data.frame(x=seq(0,1, by=0.1), y=seq(0,1, by=0.1))

methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


p200x500.4F.L <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
         coord_fixed() +
         ggtitle("200 x 500 - 4 Factors") +
         labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
         theme_bw() +
         theme(title=element_text(size=50, margin=margin(t=5)),
               axis.title.x=element_text(size=30, margin=margin(t=20)),
               axis.title.y=element_text(size=30, margin=margin(r=20)),
               legend.text=element_text(size=30)) +
         geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
         scale_size_manual(values=rep(1.1, 9)) +
         scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
         scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
         geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
         guides(colour = guide_legend(override.aes = list(size=5))) +       
         annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
                  label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))


p200x500.4F <- ggplot(data=dat.plot, aes(x=FPR, y=TPR, color=Methods)) + 
  coord_fixed() +
  ggtitle("200 x 500 - 4 Factors") +
  labs(x="False Positive Rate (1-Specificity)", y="True Positive Rate (Sensitivity)") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=TPR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2) +
  annotate("text", x=0.75, y=seq(0, 0.48, by=0.06), size=7,
           label=paste0(methods, " AUC: ", round(auc.vect, 3), "\n"))

p200x500.4F

png(filename="200x500_4F_ROC.png")

p200x500.4F

dev.off()


#-----------------#
# DET Curve Plots #
#-----------------#


dat.refline <- data.frame(x=seq(1,0, by=-0.1), y=seq(0,1, by=0.1))

dat.FPR  <- rbind(as.matrix(ttest.FPR.avg, ncol=1), 
                  as.matrix(limma.FPR.avg, ncol=1), 
                  as.matrix(RRmix.FPR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FPR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FPR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FPR.avg, ncol=1), 
                  as.matrix(PCA_limma.FPR.avg, ncol=1), 
                  as.matrix(RUV_limma.FPR.avg, ncol=1))


dat.FNR  <- rbind(as.matrix(ttest.FNR.avg, ncol=1), 
                  as.matrix(limma.FNR.avg, ncol=1), 
                  as.matrix(RRmix.FNR.avg, ncol=1), 
                  as.matrix(FAMT.NBF.FNR.avg, ncol=1), 
                  as.matrix(FAMT.DEF.FNR.avg, ncol=1), 
                  as.matrix(UNSUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(SUPSVA_limma.FNR.avg, ncol=1), 
                  as.matrix(PCA_limma.FNR.avg, ncol=1), 
                  as.matrix(RUV_limma.FNR.avg, ncol=1))


dat.plot <- data.frame(Methods=c(rep("t-test", length(ttest.FPR.avg)),
                                 rep("LIMMA", length(limma.FPR.avg)),
                                 rep("RRmix", length(RRmix.FPR.avg)),
                                 rep("FAMT NBF", length(FAMT.NBF.FPR.avg)),                                 
                                 rep("FAMT DEF", length(FAMT.DEF.FPR.avg)),                                 
                                 rep("UNSUPSVA + LIMMA", length(UNSUPSVA_limma.FPR.avg)),
                                 rep("SUPSVA + LIMMA", length(SUPSVA_limma.FPR.avg)),
                                 rep("PCA + LIMMA", length(PCA_limma.FPR.avg)),
                                 rep("RUV + LIMMA", length(RUV_limma.FPR.avg))),
                       FPR = dat.FPR,
                       FNR = dat.FNR)


methods <- c("t-test", "LIMMA", "RRmix", "FAMT NBF", "FAMT DEF", "UNSUPSVA", "SUPSVA", "PCA", "RUV")


d200x500.4F <- ggplot(data=dat.plot, aes(x=FPR, y=FNR, color=Methods)) + 
  coord_fixed(ratio = 1/3) +
  xlim(0, 0.2) +
  ylim(min(dat.plot$FNR[which(dat.plot$FPR <= 0.2)]), 
       max(dat.plot$FNR[which(dat.plot$FPR <= 0.2)])) + 
  ggtitle("200 x 500 - 4 Factors") +
  labs(x="FPR", y="FNR") +
  theme_bw() +
  theme(title=element_text(size=50, margin=margin(t=5)),
        axis.title.x=element_text(size=30, margin=margin(t=20)),
        axis.title.y=element_text(size=30, margin=margin(r=20)),
        legend.text=element_text(size=30),
        legend.position = "none") +
  geom_line(aes(x=FPR, y=FNR, color=Methods, linetype=Methods, size=Methods)) +
  scale_size_manual(values=rep(1.1, 9)) +
  scale_linetype_manual(values=c(5,4,2,8,3,9,7,1,6)) +  
  scale_color_manual(values=c(5,4,2,"darkorange",3,"darkcyan",7,1,6)) +
  geom_line(aes(x, y, color="Random Guessing"), data=dat.refline, color="gray", linetype=2)


d200x500.4F

png(filename="200x500_4F_DET.png", width=1000, height=500)

d200x500.4F

dev.off()



#--------#
# Legend #
#--------#


library(gridExtra)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

Legend <- g_legend(p200x500.4F.L)


#-----------------#
# Manuscript Plot #
#-----------------#


library(Rmisc)


## ROC Curve Plots


save(p6x265.0F, p12x265.4F, p50x265.0F, p100x265.4F, p100x265.0F, p200x265.4F, 
     p100x500.0F, p200x500.4F, Legend, file="ROC_Plots.RDATA")

lay <- rbind(c(1,2,3,4,5),
             c(6,7,8,9,5))

grid.arrange(p6x265.0F, p50x265.0F, p100x265.0F, p100x500.0F, Legend, 
             p12x265.4F, p100x265.4F, p200x265.4F, p200x500.4F,
             layout_matrix = lay)


png(filename="ROC_Plots.png", width = 3500, height = 2000)

grid.arrange(p6x265.0F, p50x265.0F, p100x265.0F, p100x500.0F, Legend, 
             p12x265.4F, p100x265.4F, p200x265.4F, p200x500.4F,
             layout_matrix = lay)

dev.off()


## DET Curve Plots


save(d6x265.0F, d12x265.4F, d50x265.0F, d100x265.4F, d100x265.0F, d200x265.4F, 
     d100x500.0F, d200x500.4F, Legend, file="DET_Plots.RDATA")

lay <- rbind(c(1,2,3,4,5),
             c(6,7,8,9,5))

grid.arrange(d6x265.0F, d50x265.0F, d100x265.0F, d100x500.0F, Legend, 
             d12x265.4F, d100x265.4F, d200x265.4F, d200x500.4F,
             layout_matrix = lay)

png(filename="DET_Plots.png", width = 3500, height = 2000)

grid.arrange(d6x265.0F, d50x265.0F, d100x265.0F, d100x500.0F, Legend, 
             d12x265.4F, d100x265.4F, d200x265.4F, d200x500.4F,
             layout_matrix = lay)

dev.off()


#--------------#
# SESSION INFO #
#--------------#

sessionInfo()