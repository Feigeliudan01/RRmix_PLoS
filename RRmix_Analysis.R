#------------------------------------------------------------------#
# Mass-Spectrometry Based Metabolomic Data - Batch Effect Analysis #
# Stephen Salerno (ss2658)                                         #
# From Data File 'filtered-II-5000-CV-50.csv'                      #
# Last Update: July 18, 2016                                       #
#------------------------------------------------------------------#


#----------------#
# Initialization #
#----------------#

set.seed(1212)

library(knitr); library(limma); library(FAMT)

data.full <- read.csv("filtered-II-5000-CV-50.csv",              # Read-In Data
                      row.names = 1,                             # Set Row Names to First Column
                      col.names = rep(c("AC","AT",               # Manually Set Column Names
                                        "XC","XT",
                                        "DC","DT",
                                        "ZC","ZT"),each=3))

data.full <- log(as.matrix(data.full))                           # Log-Transform Data


#-----------------#
# CLUSTERING PLOT #
#-----------------#

library(ggplot2)

expr.full <- t(data.full)

Batch <- rep(c("AC","AT","XC","XT","DC","DT","ZC","ZT"),each=3)     # Batch Status      

Size <- 0.1


## SVD


batch.SVD <- svd(t(expr.full) - rowMeans(t(expr.full)))$v           # SVD atch Estimation

v.SVD <- data.frame(batch.SVD)


SVD_V1xV2 <- ggplot(data=v.SVD, aes(x=v.SVD[,1], y=v.SVD[,2], color=Batch)) + 
  ggtitle("SVD: V1xV2") +
  labs(x="Right Singular Vector 1", y="Right Singular Vector 2") +
  theme_bw() +
  xlim(-0.75, 0.75) +
  ylim(-0.75, 0.75) +
  theme(title=element_text(size=20, margin=margin(t=5)),
        axis.title.x=element_text(size=15, margin=margin(t=20)),
        axis.title.y=element_text(size=15, margin=margin(r=20))) +
  geom_point(aes(x=v.SVD[,1], y=v.SVD[,2], color=Batch, size=Batch)) +
  scale_size_manual(values = rep(3,8)) +
  coord_fixed()
  
  
SVD_V2xV3 <- ggplot(data=v.SVD, aes(x=v.SVD[,2], y=v.SVD[,3], color=Batch)) + 
  ggtitle("SVD: V2xV3") +
  labs(x="Right Singular Vector 2", y="Right Singular Vector 3") +
  theme_bw() +
  xlim(-0.75, 0.75) +
  ylim(-0.75, 0.75) +
  theme(title=element_text(size=20, margin=margin(t=5)),
        axis.title.x=element_text(size=15, margin=margin(t=20)),
        axis.title.y=element_text(size=15, margin=margin(r=20))) +
  geom_point(aes(x=v.SVD[,2], y=v.SVD[,3], color=Batch, size=Batch)) +
  scale_size_manual(values = rep(3,8)) +
  coord_fixed()


SVD_V3xV4 <- ggplot(data=v.SVD, aes(x=v.SVD[,3], y=v.SVD[,4], color=Batch)) + 
  ggtitle("SVD: V3xV4") +
  labs(x="Right Singular Vector 3", y="Right Singular Vector 4") +
  theme_bw() +
  xlim(-0.75, 0.75) +
  ylim(-0.75, 0.75) +
  theme(title=element_text(size=20, margin=margin(t=5)),
        axis.title.x=element_text(size=15, margin=margin(t=20)),
        axis.title.y=element_text(size=15, margin=margin(r=20))) +
  geom_point(aes(x=v.SVD[,3], y=v.SVD[,4], color=Batch, size=Batch)) +
  scale_size_manual(values = rep(3,8)) +
  coord_fixed()


## RRmix


trmt.full  <- rep(c(0,1,0,1,0,1,0,1), each=3)                       # Treatment Status

G     <- ncol(expr.full)                                            # Number of Metabolites
n     <- nrow(expr.full)                                            # Number of Observations
Xc    <- matrix(nrow=0, ncol=0)                                     # Covariate Matrix
mu.0  <- 1/G * as.matrix(expr.full) %*% rep(1,G)                    # Initialize mu
eta.0 <- matrix(0, 2+ncol(Xc), G)                                   # Initialize eta

betac.0 <- matrix(nrow=0, ncol=0)                                   # Initialize beta_c
sig20.0 <- 1                                                        # Initialize sig^2_0
sig21.0 <- 0.1                                                      # Initialize sig^2_1

source('HEFT-RRmix-fast4-NoCovar-TB.R')                             # Source RRmix Script

t0                <- proc.time()                                    # Start Time for Script

result.full.RRmix <- runHEFTmix(G.in=G,                             # Run RRmix Model
                             n.in=n, 
                             Xc.in=Xc, 
                             Y.in=as.matrix(expr.full), 
                             SNP.in=trmt.full,
                             mu.0=mu.0, 
                             betac.0=betac.0, 
                             sig20.0=sig20.0, 
                             sig21.0=sig21.0, 
                             p.0=0.05, 
                             er_tol.in=10^(-3), 
                             q.in=4)

t.run <- proc.time()-t0                                             # Runtime for Script

save(result.full.RRmix, t.run, file = 'full-RRmix.RData')           # Save Results and Runtime


## Plot RRmix Loadings


load('full-RRmix.RData')                                           # Load Saved File
Lambda <- result.full.RRmix[['Lam']]                               # Obtain Hidden Factor Loading Matrix

Lambda <- data.frame(Lambda)


RRmix_L1xL2 <- ggplot(data=Lambda, aes(x=Lambda[,1], y=Lambda[,1], color=Batch)) + 
  ggtitle("RRmix: L1xL2") +
  labs(x="Loading 1", y="Loading 2") +
  theme_bw() +
  xlim(-0.75, 0.75) +
  ylim(-0.75, 0.75) +
  theme(title=element_text(size=20, margin=margin(t=5)),
        axis.title.x=element_text(size=15, margin=margin(t=20)),
        axis.title.y=element_text(size=15, margin=margin(r=20))) +
  geom_point(aes(x=Lambda[,1], y=Lambda[,2], color=Batch, size=Batch)) +
  scale_size_manual(values = rep(3,8)) +
  coord_fixed()


RRmix_L2xL3 <- ggplot(data=Lambda, aes(x=Lambda[,2], y=Lambda[,3], color=Batch)) + 
  ggtitle("RRmix: L2xL3") +
  labs(x="Loading 2", y="Loading 3") +
  theme_bw() +
  xlim(-0.75, 0.75) +
  ylim(-0.75, 0.75) +
  theme(title=element_text(size=20, margin=margin(t=5)),
        axis.title.x=element_text(size=15, margin=margin(t=20)),
        axis.title.y=element_text(size=15, margin=margin(r=20))) +
  geom_point(aes(x=Lambda[,2], y=Lambda[,3], color=Batch, size=Batch)) +
  scale_size_manual(values = rep(3,8)) +
  coord_fixed()


RRmix_L3xL4 <- ggplot(data=Lambda, aes(x=Lambda[,3], y=Lambda[,4], color=Batch)) + 
  ggtitle("RRmix: L3xL4") +
  labs(x="Loading 3", y="Loading 4") +
  theme_bw() +
  xlim(-0.75, 0.75) +
  ylim(-0.75, 0.75) +
  theme(title=element_text(size=20, margin=margin(t=5)),
        axis.title.x=element_text(size=15, margin=margin(t=20)),
        axis.title.y=element_text(size=15, margin=margin(r=20))) +
  geom_point(aes(x=Lambda[,3], y=Lambda[,4], color=Batch, size=Batch)) +
  scale_size_manual(values = rep(3,8)) +
  coord_fixed()


library(Rmisc)


save(SVD_V1xV2, SVD_V2xV3, SVD_V3xV4,
     RRmix_L1xL2, RRmix_L2xL3, RRmix_L3xL4, file="clustplot2.RDATA")


multiplot(SVD_V1xV2, RRmix_L1xL2, SVD_V2xV3, RRmix_L2xL3, 
          SVD_V3xV4, RRmix_L3xL4, cols=3)


png(filename="clustplot2.png", width = 1500, height = 1000)

multiplot(SVD_V1xV2, RRmix_L1xL2, SVD_V2xV3, RRmix_L2xL3, 
          SVD_V3xV4, RRmix_L3xL4, cols=3)

dev.off()



#------------#
# Operator A #
#------------#


data.A <- data.full[,c(1:6)]                                     # Subset Data for Operator A


## Number of Latent Factors


expr.A <- t(data.A)                                              # Transpose 'data.A' Matrix

G <- ncol(expr.A)                                                # Number of Metabolites
Y <- as.matrix(expr.A)                                           # Metabolite Abundance Matrix
m <- 1/(G-1)*tcrossprod(Y-tcrossprod(rowMeans(Y),rep(1,G)))      # Variance-Covariance Matrix

fit     <- princomp(covmat=m, cor=TRUE)                          # PCA on Correlation Matrix
var.exp <- fit$sdev^2/sum(fit$sdev^2)                            # Var. Explained by Prin. Comp.

## Scree Plot

plot(var.exp,                                                    # Plot Variance Explained
     type = "b",                                                 # Specify Lines and Points
     col  = 4,                                                   # Point Color
     pch  = 16,                                                  # Point Style
     xlab = "Component Index",                                   # X-Axis Label
     ylab = "% Var. Exp.",                                       # Y-Axis Label
     main = "PCA on Correlation Matrix - Subset A",              # Plot Title
     ylim = c(0,1))                                              # Y-Axis Range


text({var.exp}[1],                                               # Label Points by Value
     labels = round(as.vector(var.exp),2)[1],
     pos    = 4)

kable(t(cumsum(fit$sdev^2/sum(fit$sdev^2))[1:6]))                # Cumulative Variance Explained


## Individual t-tests


### Bonferroni Correction


G <- nrow(data.A); 
tstat <- rep(NA, G)
pval <- rep(NA, G)
df <- rep(NA, G)

for (i in seq(G)){
  tstat[i] <- t.test(data.A[i,1:3], data.A[i,4:6])$statistic
  pval[i]  <- t.test(data.A[i,1:3], data.A[i,4:6])$p.value
  df[i]    <- t.test(data.A[i,1:3], data.A[i,4:6])$parameter
}

metabolites <- rep(NA, G)
for(i in seq(G)){ metabolites[i] <- substr(rownames(data.A)[i], start=1, stop=60) }

results.A.ttest <- data.frame(metabolites=metabolites[order(pval)],
                              t.statistic=tstat[order(pval)],
                              p.value=pval[order(pval)], 
                              row.names=NULL)

sig.A.ttest.BF <- length(which(pval < 0.05/G))

cat(sig.A.ttest.BF, 'Bonferroni-corrected significant metabolites.')

kable(results.A.ttest[1:sig.A.ttest.BF,], format='markdown')


### Benjamini-Hochberg Procedure (FDR = 0.1)


p.sorted <- sort(pval)
FDR <- 0.1
BH <- FDR*(1:length(p.sorted))/length(p.sorted)
sig.A.ttest.BH <- length(which(p.sorted <= BH))

cat(sig.A.ttest.BH, 'significant metabolites at FDR = 0.1.')

mets.t.A <- as.vector(results.A.ttest[1:sig.A.ttest.BH, 1])

write.table(mets.t.A, "mets_t_A.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(results.A.ttest[1:sig.A.ttest.BH,], format='markdown')


## RRmix Model Analysis


trmt.A  <- rep(c(0,1), each=3)                                   # Treatment Status
batch.A <- rep(c(1,2), each=3)                                   # Batch Status

G     <- ncol(expr.A)                                            # Number of Metabolites
n     <- nrow(expr.A)                                            # Number of Observations
Xc    <- matrix(nrow=0, ncol=0)                                  # Covariate Matrix
mu.0  <- 1/G * as.matrix(expr.A) %*% rep(1,G)                    # Initialize mu
eta.0 <- matrix(0, 2+ncol(Xc), G)                                # Initialize eta

betac.0 <- matrix(nrow=0, ncol=0)                                # Initialize beta_c
sig20.0 <- 1                                                     # Initialize sig^2_0
sig21.0 <- 0.1                                                   # Initialize sig^2_1

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

t0             <- proc.time()                                    # Start Time for Script

result.A.RRmix <- runHEFTmix(G.in=G,                             # Run RRmix Model
                             n.in=n, 
                             Xc.in=Xc, 
                             Y.in=as.matrix(expr.A), 
                             SNP.in=trmt.A,
                             mu.0=mu.0, 
                             betac.0=betac.0, 
                             sig20.0=sig20.0, 
                             sig21.0=sig21.0, 
                             p.0=0.05, 
                             er_tol.in=10^(-3), 
                             q.in=2)

t.run <- proc.time()-t0                                          # Runtime for Script

save(result.A.RRmix, t.run, file = 'A-RRmix.RData')              # Save Results and Runtime


## RRmix Model Results


load('A-RRmix.RData')                                            # Load Saved File

t.run                                                            # Display Runtime for Script

b_g           <- result.A.RRmix[['b_g']]                         # Extract Posterior Probabilities
sig.RRmix.A.9 <- length(which(b_g > 0.9))                        # Number of Significant Metabolites (0.9)

cat(sig.RRmix.A.9,                                               # Display Number of Metabolites
    'significant metabolites with post. prob. = 0.9.')

toptable.RRmix.A <- data.frame(metabolite = metabolites[order(-b_g)], 
                               post.prob  = b_g[order(-b_g)], 
                               row.names=NULL)

mets.RRmix.A <- toptable.RRmix.A[1:sig.RRmix.A.9, 1]

write.table(mets.RRmix.A, "mets_RRmix_A.txt", sep="\t", row.names=FALSE, col.names=FALSE)


kable(toptable.RRmix.A[1:sig.RRmix.A.9,],                        # Display Table (0.9)
      format='markdown')     

sig.RRmix.A.8 <- length(which(b_g>0.8))                          # Number of Significant Metabolites (0.8)      

cat(sig.RRmix.A.8,                                               # Displat Number of Metabolites
    'significant metabolites with post. prob. = 0.8.')       

kable(toptable.RRmix.A[1:sig.RRmix.A.8,],                        # Display Table (0.8)
      format='markdown')     


## Plot RRmix Loadings


load('A-RRmix.RData')                                           # Load Saved File
Lambda <- result.A.RRmix[['Lam']]                               # Obtain Hidden Factor Loading Matrix

plot(Lambda[,1],                                                # Plot First Loading
     Lambda[,2],                                                # Plot Second Loading
     xlab="Loading-1",                                          # X-Axis Label
     ylab="Loading-2",                                          # Y-Axis Label
     main="RRmix Analysis - Operator A")                        # Plot Title

points(Lambda[which(batch.A==1),1],                             # Color Points from First Batch Red
       Lambda[which(batch.A==1),2],                             
       pch=16, 
       col="red")

points(Lambda[which(batch.A==2),1],                             # Color Points from Second Batch Blue
       Lambda[which(batch.A==2),2], 
       pch=16, 
       col="blue")

plot.new(); legend("center",                                    # Create Legend
                   paste(c("AC","AT")), 
                   pch=16, 
                   col=c("red","blue"))


## LIMMA Model Analysis


### Bonferroni Correction


expr.A  <- t(data.A)                                             # Set Abundance Matrix
trmt.A  <- rep(c(0,1),each=3)                                    # Set Treatment Groups
batch.A <- rep(c(1,2), each=3)                                   # Set Batch Groups

t0 <- proc.time()                                                # Start Time for Script

design   <- model.matrix(~factor(trmt.A))                        # Create Design Matrix
fit      <- lmFit(t(expr.A),design)                              # Fit Model
ebayes.A <- eBayes(fit)                                          # Get Statistics

t.run <- proc.time()-t0                                          # Runtime for Script

save(ebayes.A, t.run, file = 'A-LIMMA.RData')                    # Save Results and Runtime

load('A-LIMMA.RData')                                            # Load Saved File 

t.run                                                            # Display Runtime for Script

tab.A.BF <- topTable(ebayes.A,                                   # Bonferroni-Corrected Toptable
                     coef=2, 
                     number=G, 
                     adjust="bonferroni")

toptable.A.LIMMA.BF <- data.frame(metabolites = substr(rownames(tab.A.BF),1,50),
                                  pvalue = tab.A.BF[,'adj.P.Val'], 
                                  row.names = NULL)

kable(toptable.A.LIMMA.BF[which(tab.A.BF[,'adj.P.Val'] < 0.05), ])    

heatmap(t(expr.A)[which(tab.A.BF[,'adj.P.Val'] < 0.05),])        # Display Heatmap


### Benjamini-Hochberg Procedure: FDR=0.1


tab.A.BH <- topTable(ebayes.A,                                   # FDR Threshold Toptable
                     coef=2, 
                     number=G, 
                     adjust="BH")

toptable.A.LIMMA.BH <- data.frame(metabolites = substr(rownames(tab.A.BH),1,50),
                                  pvalue = tab.A.BH[,'adj.P.Val'], 
                                  row.names = NULL)

mets.limma.A <- toptable.A.LIMMA.BH[1:length(which(tab.A.BH[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.limma.A, "mets_limma_A.txt", sep="\t", row.names=FALSE, col.names=FALSE)


kable(toptable.A.LIMMA.BH[1:length(which(tab.A.BH[,'adj.P.Val'] < 0.1)),],
      format='markdown')

heatmap(t(expr.A)[which(tab.A.BH[,'adj.P.Val'] < 0.1),])         # Display Heatmap


## FAMT Model Analysis


expr.A.FAMT <- data.A                                            # Set Data Matrix

cov.A.FAMT  <- data.frame(id    = rownames(expr.A),              # Set Covariates Matrix
                          trmt  = as.factor(trmt.A), 
                          batch = as.factor(batch.A))

data.A.FAMT <- as.FAMTdata(expression = expr.A.FAMT,             # Make Data Structure
                           covariates = cov.A.FAMT, 
                           idcovar    = 1)

fit.A.FAMT  <- modelFAMT(data.A.FAMT,                            # Fit FAMT Model
                         x    = 2, 
                         test = 2, 
                         nbf  = 2)


toptable.A.FAMT <- data.frame(metabolites = metabolites[order(fit.A.FAMT$adjpval)],
                              pvalues     = fit.A.FAMT$adjpval[order(fit.A.FAMT$adjpval)],
                              row.names   = NULL)

sig.A.FAMT      <- length(which(fit.A.FAMT$adjpval < 0.1))       # Number of Significant Metabolites (BH)

cat(sig.A.FAMT,                                                  # Display Number of Metabolites (BH)
    ' significant metabolites at FDR = 0.1.')

mets.FAMT.A <- toptable.A.FAMT[1:sig.A.FAMT, 1]

write.table(mets.FAMT.A, "mets_FAMT_A.txt", sep="\t", row.names=FALSE, col.names=FALSE)


kable(toptable.A.FAMT[1:sig.A.FAMT,],                            # Display Significant Metabolites (BH)
      row.names=NA)


## PCA + LIMMA


batch.A.PCA <- svd(t(expr.A) - rowMeans(t(expr.A)))$v[,1]                       # Batch Estimation

mod.A.PCA <- model.matrix(~trmt.A+batch.A.PCA)                                  # Model Matrix

fit.A.PCA <- lmFit(t(expr.A), mod.A.PCA)                                        # Fit Model

ebayes.A.PCA <- eBayes(fit.A.PCA)                                               # Calculate Statistics
 
tab.A.PCA <- toptable(fit.A.PCA, coef=2,                                        # Toptable
                      number=ncol(expr.A), p.value=0.1)     

cat(toString(nrow(tab.A.PCA)), "significant metabolites found by PCA.")         # Print Statement

toptable.A.PCA.BH <- data.frame(metabolites = substr(rownames(tab.A.PCA),1,50), # Toptable
                                  pvalue = tab.A.PCA[,'adj.P.Val'], 
                                  row.names = NULL)

mets.PCA.A <- toptable.A.PCA.BH[1:length(which(tab.A.PCA[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.PCA.A, "mets_PCA_A.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.A.PCA.BH[1:length(which(tab.A.PCA[,'adj.P.Val'] < 0.1)),],       # Print Statement
      format='markdown')


## SVA + LIMMA


library(sva)

mod  <- model.matrix(~trmt.A)                                             # Initial Model Matrix

mod0 <- cbind(mod[,1])                                                    # Initial Model Matrix

batch.A.SVA <- sva(t(expr.A), mod, mod0)                                  # Batch Estimation

mod.A.SVA <- cbind(mod, batch.A.SVA$sv)                                   # SVA Model Matrix

fit.A.SVA <- lmFit(t(expr.A), mod.A.SVA)                                  # Fit Model

ebayes.A.SVA <- eBayes(fit.A.SVA)                                         # Calculate Statistics

tab.A.SVA <- toptable(fit.A.SVA, coef=2,                                  # Toptable
                      number=ncol(expr.A), p.value=0.1)      

cat(toString(nrow(tab.A.SVA)), "significant metabolites found by SVA.")   # Print Statement

toptable.A.SVA.BH <- data.frame(metabolites = substr(rownames(tab.A.SVA),1,50), # Toptable
                                pvalue = tab.A.SVA[,'adj.P.Val'], 
                                row.names = NULL)

mets.SVA.A <- toptable.A.SVA.BH[1:length(which(tab.A.SVA[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.SVA.A, "mets_SVA_A.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.A.SVA.BH[1:length(which(tab.A.SVA[,'adj.P.Val'] < 0.1)),],       # Print Statement
      format='markdown')


#---------------#
# Venn Diagrams #
#---------------#


discoveries.A <- list("t-tests"=mets.t.A, "LIMMA"=mets.limma.A, "SVA + LIMMA"=mets.SVA.A, 
                      "PCA + LIMMA"=mets.PCA.A, "RRmix"=mets.RRmix.A, "FAMT"=mets.FAMT.A)

Venn.A <- Venn(discoveries.A)

Venn.A <- compute.Venn(Venn.A, type="AWFE")

SetLabels.A <- VennGetSetLabels(Venn.A)

SetLabels.A[, "x"] <- c(-3.8, -4, 1, 2, -1.4, -2)
SetLabels.A[, "y"] <- c(3.6, 4.3, 2, 1.55, 2.7, 2)

Venn.A <- VennSetSetLabels(Venn.A, SetLabels.A)

gp.A <- VennThemes(Venn.A)

gp.A[["SetText"]][["Set1"]]$col <- 1
gp.A[["SetText"]][["Set2"]]$col <- 2
gp.A[["SetText"]][["Set3"]]$col <- 6
gp.A[["SetText"]][["Set4"]]$col <- "darkorange"
gp.A[["SetText"]][["Set5"]]$col <- 3
gp.A[["SetText"]][["Set6"]]$col <- 5

gp.A[["Set"]][["Set1"]]$col <- 1
gp.A[["Set"]][["Set2"]]$col <- 2
gp.A[["Set"]][["Set3"]]$col <- 6
gp.A[["Set"]][["Set4"]]$col <- "darkorange"
gp.A[["Set"]][["Set5"]]$col <- 3
gp.A[["Set"]][["Set6"]]$col <- 5

plot(Venn.A, gp=gp.A)

pdf("VennA.pdf")

plot(Venn.A, gp=gp.A)

dev.off()


discoveries.A.4 <- list("t-tests"=mets.t.A, "LIMMA"=mets.limma.A, 
                        "RRmix"=mets.RRmix.A, "FAMT"=mets.FAMT.A)

Venn.A.4 <- Venn(discoveries.A.4)

Venn.A.4 <- compute.Venn(Venn.A.4, type="ellipses", doWeights=FALSE)

SetLabels.A.4 <- VennGetSetLabels(Venn.A.4)

SetLabels.A.4[, "x"] <- c(9, -11, 7, -9)
SetLabels.A.4[, "y"] <- c(6, 6, 8, 8)

Venn.A.4 <- VennSetSetLabels(Venn.A.4, SetLabels.A.4)

gp.A.4 <- VennThemes(Venn.A.4)

gp.A.4[["SetText"]][["Set1"]]$col <- 1
gp.A.4[["SetText"]][["Set2"]]$col <- 2
gp.A.4[["SetText"]][["Set3"]]$col <- 3
gp.A.4[["SetText"]][["Set4"]]$col <- 5

gp.A.4[["Set"]][["Set1"]]$col <- 1
gp.A.4[["Set"]][["Set2"]]$col <- 2
gp.A.4[["Set"]][["Set3"]]$col <- 3
gp.A.4[["Set"]][["Set4"]]$col <- 5

plot(Venn.A.4, gp=gp.A.4)

pdf("VennA_4.pdf")

plot(Venn.A.4, gp=gp.A.4)

dev.off()


#------------#
# Operator X #
#------------#


data.X <- data.full[,c(7:12)]                                    # Subset Data for Operator X


## Number of Latent Factors


expr.X <- t(data.X)                                              # Transpose 'data.X' Matrix

G <- ncol(expr.X)                                                # Number of Metabolites
Y <- as.matrix(expr.X)                                           # Metabolite Abundance Matrix
m <- 1/(G-1)*tcrossprod(Y-tcrossprod(rowMeans(Y),rep(1,G)))      # Variance-Covariance Matrix

fit     <- princomp(covmat=m, cor=TRUE)                          # PCA on Correlation Matrix
var.exp <- fit$sdev^2/sum(fit$sdev^2)                            # Var. Explained by Prin. Comp.

## Scree Plot

plot(var.exp,                                                    # Plot Variance Explained
     type = "b",                                                 # Specify Lines and Points
     col  = 4,                                                   # Point Color
     pch  = 16,                                                  # Point Style
     xlab = "Component Index",                                   # X-Axis Label
     ylab = "% Var. Exp.",                                       # Y-Axis Label
     main = "PCA on Correlation Matrix - Subset X",              # Plot Title
     ylim = c(0,1))                                              # Y-Axis Range


text({var.exp}[1],                                               # Label Points by Value
     labels = round(as.vector(var.exp),2)[1],
     pos    = 4)

kable(t(cumsum(fit$sdev^2/sum(fit$sdev^2))[1:6]))                # Cumulative Variance Explained


## Individual t-tests


### Bonferroni Correction


G <- nrow(data.X); 
tstat <- rep(NA, G)
pval <- rep(NA, G)
df <- rep(NA, G)

for (i in seq(G)){
  tstat[i] <- t.test(data.X[i,1:3], data.X[i,4:6])$statistic
  pval[i]  <- t.test(data.X[i,1:3], data.X[i,4:6])$p.value
  df[i]    <- t.test(data.X[i,1:3], data.X[i,4:6])$parameter
}

metabolites <- rep(NA, G)
for(i in seq(G)){ metabolites[i] <- substr(rownames(data.X)[i], start=1, stop=60) }

results.X.ttest <- data.frame(metabolites=metabolites[order(pval)],
                              t.statistic=tstat[order(pval)],
                              p.value=pval[order(pval)], 
                              row.names=NULL)

sig.X.ttest.BF <- length(which(pval < 0.05/G))

cat(sig.X.ttest.BF, 'Bonferroni-corrected significant metabolites.')

kable(results.X.ttest[1:sig.X.ttest.BF,], format='markdown')


### Benjamini-Hochberg Procedure (FDR = 0.1)


p.sorted <- sort(pval)
FDR <- 0.1
BH <- FDR*(1:length(p.sorted))/length(p.sorted)
sig.X.ttest.BH <- length(which(p.sorted <= BH))

cat(sig.X.ttest.BH, 'significant metabolites at FDR = 0.1.')

mets.t.X <- as.vector(results.X.ttest[1:sig.X.ttest.BH, 1])

write.table(mets.t.X, "mets_t_X.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(results.X.ttest[1:sig.X.ttest.BH,], format='markdown')


## RRmix Model Analysis


trmt.X  <- rep(c(0,1), each=3)                                   # Treatment Status
batch.X <- rep(c(1,2), each=3)                                   # Batch Status

G     <- ncol(expr.X)                                            # Number of Metabolites
n     <- nrow(expr.X)                                            # Number of Observations
Xc    <- matrix(nrow=0, ncol=0)                                  # Covariate Matrix
mu.0  <- 1/G * as.matrix(expr.X) %*% rep(1,G)                    # Initialize mu
eta.0 <- matrix(0, 2+ncol(Xc), G)                                # Initialize eta

betac.0 <- matrix(nrow=0, ncol=0)                                # Initialize beta_c
sig20.0 <- 1                                                     # Initialize sig^2_0
sig21.0 <- 0.1                                                   # Initialize sig^2_1

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

t0             <- proc.time()                                    # Start Time for Script

result.X.RRmix <- runHEFTmix(G.in=G,                             # Run RRmix Model
                             n.in=n, 
                             Xc.in=Xc, 
                             Y.in=as.matrix(expr.X), 
                             SNP.in=trmt.X,
                             mu.0=mu.0, 
                             betac.0=betac.0, 
                             sig20.0=sig20.0, 
                             sig21.0=sig21.0, 
                             p.0=0.05, 
                             er_tol.in=10^(-3), 
                             q.in=2)

t.run <- proc.time()-t0                                          # Runtime for Script

save(result.X.RRmix, t.run, file = 'X-RRmix.RData')              # Save Results and Runtime


## RRmix Model Results


load('X-RRmix.RData')                                            # Load Saved File

t.run                                                            # Display Runtime for Script

b_g           <- result.X.RRmix[['b_g']]                         # Extract Posterior Probabilities
sig.RRmix.X.9 <- length(which(b_g > 0.9))                        # Number of Significant Metabolites (0.9)

cat(sig.RRmix.X.9,                                               # Display Number of Metabolites
    'significant metabolites with post. prob. = 0.9.')

toptable.RRmix.X <- data.frame(metabolite = metabolites[order(-b_g)], 
                               post.prob  = b_g[order(-b_g)], 
                               row.names=NULL)

mets.RRmix.X <- toptable.RRmix.X[1:sig.RRmix.X.9, 1]

write.table(mets.RRmix.X, "mets_RRmix_X.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.RRmix.X[1:sig.RRmix.X.9,],                        # Display Table (0.9)
      format='markdown')     

sig.RRmix.X.8 <- length(which(b_g>0.8))                          # Number of Significant Metabolites (0.8)      

cat(sig.RRmix.X.8,                                               # Displat Number of Metabolites
    'significant metabolites with post. prob. = 0.8.')       

kable(toptable.RRmix.X[1:sig.RRmix.X.8,],                        # Display Table (0.8)
      format='markdown')     


## Plot RRmix Loadings


load('X-RRmix.RData')                                           # Load Saved File
Lambda <- result.X.RRmix[['Lam']]                               # Obtain Hidden Factor Loading Matrix

plot(Lambda[,1],                                                # Plot First Loading
     Lambda[,2],                                                # Plot Second Loading
     xlab="Loading-1",                                          # X-Axis Label
     ylab="Loading-2",                                          # Y-Axis Label
     main="RRmix Analysis - Operator X")                        # Plot Title

points(Lambda[which(batch.X==1),1],                             # Color Points from First Batch Red
       Lambda[which(batch.X==1),2],                             
       pch=16, 
       col="red")

points(Lambda[which(batch.X==2),1],                             # Color Points from Second Batch Blue
       Lambda[which(batch.X==2),2], 
       pch=16, 
       col="blue")

plot.new(); legend("center",                                    # Create Legend
                   paste(c("XC","XT")), 
                   pch=16, 
                   col=c("red","blue"))


## LIMMA Model Analysis


### Bonferroni Correction


expr.X  <- t(data.X)                                             # Set Abundance Matrix
trmt.X  <- rep(c(0,1),each=3)                                    # Set Treatment Groups
batch.X <- rep(c(1,2), each=3)                                   # Set Batch Groups

t0 <- proc.time()                                                # Start Time for Script

design   <- model.matrix(~factor(trmt.X))                        # Create Design Matrix
fit      <- lmFit(t(expr.X),design)                              # Fit Model
ebayes.X <- eBayes(fit)                                          # Get Statistics

t.run <- proc.time()-t0                                          # Runtime for Script

save(ebayes.X, t.run, file = 'X-LIMMA.RData')                    # Save Results and Runtime

load('X-LIMMA.RData')                                            # Load Saved File 

t.run                                                            # Display Runtime for Script

tab.X.BF <- topTable(ebayes.X,                                   # Bonferroni-Corrected Toptable
                     coef=2, 
                     number=G, 
                     adjust="bonferroni")

toptable.X.LIMMA.BF <- data.frame(metabolites = substr(rownames(tab.X.BF),1,50),
                                  pvalue = tab.X.BF[,'adj.P.Val'], 
                                  row.names = NULL)

kable(toptable.X.LIMMA.BF[which(tab.X.BF[,'adj.P.Val'] < 0.05), ])  

heatmap(t(expr.X)[which(tab.X.BF[,'adj.P.Val'] < 0.05),])        # Display Heatmap


### Benjamini-Hochberg Procedure: FDR=0.1


tab.X.BH <- topTable(ebayes.X,                                   # FDR Threshold Toptable
                     coef=2, 
                     number=G, 
                     adjust="BH")

toptable.X.LIMMA.BH <- data.frame(metabolites = substr(rownames(tab.X.BH),1,50),
                                  pvalue = tab.X.BH[,'adj.P.Val'], 
                                  row.names = NULL)

mets.limma.X <- toptable.X.LIMMA.BH[1:length(which(tab.X.BH[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.limma.X, "mets_limma_X.txt", sep="\t", row.names=FALSE, col.names=FALSE)


kable(toptable.X.LIMMA.BH[1:length(which(tab.X.BH[,'adj.P.Val'] < 0.1)),],
      format='markdown')

heatmap(t(expr.X)[which(tab.X.BH[,'adj.P.Val'] < 0.1),])         # Display Heatmap


## FAMT Model Analysis


expr.X.FAMT <- data.X                                            # Set Data Matrix

cov.X.FAMT  <- data.frame(id    = rownames(expr.X),              # Set Covariates Matrix
                          trmt  = as.factor(trmt.X), 
                          batch = as.factor(batch.X))

data.X.FAMT <- as.FAMTdata(expression = expr.X.FAMT,             # Make Data Structure
                           covariates = cov.X.FAMT, 
                           idcovar    = 1)

fit.X.FAMT  <- modelFAMT(data.X.FAMT,                            # Fit FAMT Model
                         x    = 2, 
                         test = 2, 
                         nbf  = 2)


toptable.X.FAMT <- data.frame(metabolites = metabolites[order(fit.X.FAMT$adjpval)],
                              pvalues     = fit.X.FAMT$adjpval[order(fit.X.FAMT$adjpval)],
                              row.names   = NULL)

sig.X.FAMT      <- length(which(fit.X.FAMT$adjpval < 0.1))       # Number of Significant Metabolites (BH)

cat(sig.X.FAMT,                                                  # Display Number of Metabolites (BH)
    ' significant metabolites at FDR = 0.1.')

mets.FAMT.X <- toptable.X.FAMT[1:sig.X.FAMT, 1]

write.table(mets.FAMT.X, "mets_FAMT_X.txt", sep="\t", row.names=FALSE, col.names=FALSE)


kable(toptable.X.FAMT[1:sig.X.FAMT,],                            # Display Significant Metabolites (BH)
      row.names=NA)


## PCA + LIMMA


batch.X.PCA <- svd(t(expr.X) - rowMeans(t(expr.X)))$v[,1]                       # Batch Estimation

mod.X.PCA <- model.matrix(~trmt.X+batch.X.PCA)                                  # Model Matrix

fit.X.PCA <- lmFit(t(expr.X), mod.X.PCA)                                        # Fit Model

ebayes.X.PCA <- eBayes(fit.X.PCA)                                               # Calculate Statistics

tab.X.PCA <- toptable(fit.X.PCA, coef=2,                                        # Toptable
                      number=ncol(expr.X), p.value=0.1)     

cat(toString(nrow(tab.X.PCA)), "significant metabolites found by PCA.")         # Print Statement

toptable.X.PCA.BH <- data.frame(metabolites = substr(rownames(tab.X.PCA),1,50), # Toptable
                                pvalue = tab.X.PCA[,'adj.P.Val'], 
                                row.names = NULL)

mets.PCA.X <- toptable.X.PCA.BH[1:length(which(tab.X.PCA[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.PCA.X, "mets_PCA_X.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.X.PCA.BH[1:length(which(tab.X.PCA[,'adj.P.Val'] < 0.1)),],       # Print Statement
      format='markdown')


## SVA + LIMMA


library(sva)

mod  <- model.matrix(~trmt.X)                                             # Initial Model Matrix

mod0 <- cbind(mod[,1])                                                    # Initial Model Matrix

batch.X.SVA <- sva(t(expr.X), mod, mod0)                                  # Batch Estimation

mod.X.SVA <- cbind(mod, batch.X.SVA$sv)                                   # SVA Model Matrix

fit.X.SVA <- lmFit(t(expr.X), mod.X.SVA)                                  # Fit Model

ebayes.X.SVA <- eBayes(fit.X.SVA)                                         # Calculate Statistics

tab.X.SVA <- toptable(fit.X.SVA, coef=2,                                  # Toptable
                      number=ncol(expr.X), p.value=0.1)      

cat(toString(nrow(tab.X.SVA)), "significant metabolites found by SVA.")   # Print Statement

toptable.X.SVA.BH <- data.frame(metabolites = substr(rownames(tab.X.SVA),1,50), # Toptable
                                pvalue = tab.X.SVA[,'adj.P.Val'], 
                                row.names = NULL)

mets.SVA.X <- toptable.X.SVA.BH[1:length(which(tab.X.SVA[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.SVA.X, "mets_SVA_X.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.X.SVA.BH[1:length(which(tab.X.SVA[,'adj.P.Val'] < 0.1)),],       # Print Statement
      format='markdown')


#---------------#
# Venn Diagrams # 
#---------------#


discoveries.X <- list("t-tests"=mets.t.X, "LIMMA"=mets.limma.X, "SVA + LIMMA"=mets.SVA.X, 
                      "PCA + LIMMA"=mets.PCA.X, "RRmix"=mets.RRmix.X, "FAMT"=mets.FAMT.X)

Venn.X <- Venn(discoveries.X)

Venn.X <- compute.Venn(Venn.X, type="AWFE")

SetLabels.X <- VennGetSetLabels(Venn.X)

SetLabels.X[, "x"] <- c(-3.8, -4, 1, 2, -1.4, -2)
SetLabels.X[, "y"] <- c(3.6, 4.3, 2, 1.55, 2.7, 2)

Venn.X <- VennSetSetLabels(Venn.X, SetLabels.X)

gp.X <- VennThemes(Venn.X)

gp.X[["SetText"]][["Set1"]]$col <- 1
gp.X[["SetText"]][["Set2"]]$col <- 2
gp.X[["SetText"]][["Set3"]]$col <- 6
gp.X[["SetText"]][["Set4"]]$col <- "darkorange"
gp.X[["SetText"]][["Set5"]]$col <- 3
gp.X[["SetText"]][["Set6"]]$col <- 5

gp.X[["Set"]][["Set1"]]$col <- 1
gp.X[["Set"]][["Set2"]]$col <- 2
gp.X[["Set"]][["Set3"]]$col <- 6
gp.X[["Set"]][["Set4"]]$col <- "darkorange"
gp.X[["Set"]][["Set5"]]$col <- 3
gp.X[["Set"]][["Set6"]]$col <- 5

plot(Venn.X, gp=gp.X)

pdf("VennX.pdf")

plot(Venn.X, gp=gp.X)

dev.off()


discoveries.X.4 <- list("t-tests"=mets.t.X, "LIMMA"=mets.limma.X, 
                         "RRmix"=mets.RRmix.X, "FAMT"=mets.FAMT.X)

Venn.X.4 <- Venn(discoveries.X.4)

Venn.X.4 <- compute.Venn(Venn.X.4, type="ellipses", doWeights=FALSE)

SetLabels.X.4 <- VennGetSetLabels(Venn.X.4)

SetLabels.X.4[, "x"] <- c(9, -11, 7, -9)
SetLabels.X.4[, "y"] <- c(6, 6, 8, 8)

Venn.X.4 <- VennSetSetLabels(Venn.X.4, SetLabels.X.4)

gp.X.4 <- VennThemes(Venn.X.4)

gp.X.4[["SetText"]][["Set1"]]$col <- 1
gp.X.4[["SetText"]][["Set2"]]$col <- 2
gp.X.4[["SetText"]][["Set3"]]$col <- 3
gp.X.4[["SetText"]][["Set4"]]$col <- 5

gp.X.4[["Set"]][["Set1"]]$col <- 1
gp.X.4[["Set"]][["Set2"]]$col <- 2
gp.X.4[["Set"]][["Set3"]]$col <- 3
gp.X.4[["Set"]][["Set4"]]$col <- 5

plot(Venn.X.4, gp=gp.X.4)

pdf("VennX_4.pdf")

plot(Venn.X.4, gp=gp.X.4)

dev.off()



#-----------------#
# Operators A & X #
#-----------------#


data.AX <- data.full[,c(1:12)]


## Number of Latent Factors


expr.AX <- t(data.AX)                                            # Transpose 'data.AX' Matrix

G <- ncol(expr.AX)                                               # Number of Metabolites
Y <- as.matrix(expr.AX)                                          # Metabolite Abundance Matrix
m <- 1/(G-1)*tcrossprod(Y-tcrossprod(rowMeans(Y),rep(1,G)))      # Variance-Covariance Matrix

fit     <- princomp(covmat=m, cor=TRUE)                          # PCA on Correlation Matrix
var.exp <- fit$sdev^2/sum(fit$sdev^2)                            # Var. Explained by Prin. Comp.

## Scree Plot

plot(var.exp,                                                    # Plot Variance Explained
     type = "b",                                                 # Specify Lines and Points
     col  = 4,                                                   # Point Color
     pch  = 16,                                                  # Point Style
     xlab = "Component Index",                                   # X-Axis Label
     ylab = "% Var. Exp.",                                       # Y-Axis Label
     main = "PCA on Correlation Matrix - Subset AX",             # Plot Title
     ylim = c(0,1))                                              # Y-Axis Range


text({var.exp}[1],                                               # Label Points by Value
     labels = round(as.vector(var.exp),2)[1],
     pos    = 4)

kable(t(cumsum(fit$sdev^2/sum(fit$sdev^2))[1:6]))                # Cumulative Variance Explained


## Individual t-tests


### Bonferroni Correction


G <- nrow(data.AX); 
tstat <- rep(NA, G)
pval <- rep(NA, G)
df <- rep(NA, G)

for (i in seq(G)){
  tstat[i] <- t.test(data.AX[i,c(1:3,7:9)], data.AX[i,c(4:6,10:12)])$statistic
  pval[i]  <- t.test(data.AX[i,c(1:3,7:9)], data.AX[i,c(4:6,10:12)])$p.value
  df[i]    <- t.test(data.AX[i,c(1:3,7:9)], data.AX[i,c(4:6,10:12)])$parameter
}

metabolites <- rep(NA, G)
for(i in seq(G)){ metabolites[i] <- substr(rownames(data.AX)[i], start=1, stop=60) }

results.AX.ttest <- data.frame(metabolites=metabolites[order(pval)],
                               t.statistic=tstat[order(pval)],
                               p.value=pval[order(pval)], 
                               row.names=NULL)

sig.AX.ttest.BF <- length(which(pval < 0.05/G))

cat(sig.AX.ttest.BF, 'Bonferroni-corrected significant metabolites.')

kable(results.AX.ttest[1:sig.AX.ttest.BF,], format='markdown')


### Benjamini-Hochberg Procedure (FDR = 0.1)


p.sorted <- sort(pval)
FDR <- 0.1
BH <- FDR*(1:length(p.sorted))/length(p.sorted)
sig.AX.ttest.BH <- length(which(p.sorted <= BH))

cat(sig.AX.ttest.BH, 'significant metabolites at FDR = 0.1.')

mets.t.AX <- as.vector(results.AX.ttest[1:sig.AX.ttest.BH, 1])

write.table(mets.t.AX, "mets_t_AX.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(results.AX.ttest[1:sig.AX.ttest.BH,], format='markdown')


## RRmix Model Analysis


trmt.AX  <- rep(c(0,1,0,1), each=3)                              # Treatment Status
batch.AX <- rep(c(1,2,3,4), each=3)                              # Batch Status

G     <- ncol(expr.AX)                                           # Number of Metabolites
n     <- nrow(expr.AX)                                           # Number of Observations
Xc    <- matrix(nrow=0, ncol=0)                                  # Covariate Matrix
mu.0  <- 1/G * as.matrix(expr.AX) %*% rep(1,G)                   # Initialize mu
eta.0 <- matrix(0, 2+ncol(Xc), G)                                # Initialize eta

betac.0 <- matrix(nrow=0, ncol=0)                                # Initialize beta_c
sig20.0 <- 1                                                     # Initialize sig^2_0
sig21.0 <- 0.1                                                   # Initialize sig^2_1

source('HEFT-RRmix-fast4-NoCovar-TB.R')                          # Source RRmix Script

t0             <- proc.time()                                    # Start Time for Script

result.AX.RRmix <- runHEFTmix(G.in=G,                            # Run RRmix Model
                              n.in=n, 
                              Xc.in=Xc, 
                              Y.in=as.matrix(expr.AX), 
                              SNP.in=trmt.AX,
                              mu.0=mu.0, 
                              betac.0=betac.0, 
                              sig20.0=sig20.0, 
                              sig21.0=sig21.0, 
                              p.0=0.05, 
                              er_tol.in=10^(-3), 
                              q.in=1)

t.run <- proc.time()-t0                                          # Runtime for Script

save(result.AX.RRmix, t.run, file = 'AX-RRmix.RData')            # Save Results and Runtime


## RRmix Model Results


load('AX-RRmix.RData')                                           # Load Saved File

t.run                                                            # Display Runtime for Script

b_g           <- result.AX.RRmix[['b_g']]                        # Extract Posterior Probabilities
sig.RRmix.AX.9 <- length(which(b_g > 0.9))                       # Number of Significant Metabolites (0.9)

cat(sig.RRmix.AX.9,                                              # Display Number of Metabolites
    'significant metabolites with post. prob. = 0.9.')

toptable.RRmix.AX <- data.frame(metabolite = metabolites[order(-b_g)], 
                                post.prob  = b_g[order(-b_g)], 
                                row.names=NULL)


mets.RRmix.AX <- toptable.RRmix.AX[1:sig.RRmix.AX.9, 1]

write.table(mets.RRmix.AX, "mets_RRmix_AX.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.RRmix.AX[1:sig.RRmix.AX.9,],                      # Display Table (0.9)
      format='markdown')     

sig.RRmix.AX.8 <- length(which(b_g>0.8))                         # Number of Significant Metabolites (0.8)      

cat(sig.RRmix.AX.8,                                              # Displat Number of Metabolites
    'significant metabolites with post. prob. = 0.8.')       

kable(toptable.RRmix.AX[1:sig.RRmix.AX.8,],                      # Display Table (0.8)
      format='markdown')     


## LIMMA Model Analysis


### Bonferroni Correction


expr.AX  <- t(data.AX)                                           # Set Abundance Matrix

trmt.AX  <- rep(c(0,1,0,1),each=3)                               # Set Treatment Groups
batch.AX <- rep(c(1,2,3,4), each=3)                              # Set Batch Groups

t0 <- proc.time()                                                # Start Time for Script

design   <- model.matrix(~factor(trmt.AX))                       # Create Design Matrix
fit      <- lmFit(t(expr.AX),design)                             # Fit Model
ebayes.AX <- eBayes(fit)                                         # Get Statistics

t.run <- proc.time()-t0                                          # Runtime for Script

save(ebayes.AX, t.run, file = 'AX-LIMMA.RData')                  # Save Results and Runtime

load('AX-LIMMA.RData')                                           # Load Saved File 

t.run                                                            # Display Runtime for Script

tab.AX.BF <- topTable(ebayes.AX,                                 # Bonferroni-Corrected Toptable
                      coef=2, 
                      number=G, 
                      adjust="bonferroni")

toptable.AX.LIMMA.BF <- data.frame(metabolites = substr(rownames(tab.AX.BF),1,50),
                                   pvalue = tab.AX.BF[,'adj.P.Val'], 
                                   row.names = NULL)

kable(toptable.AX.LIMMA.BF[which(tab.AX.BF[,'adj.P.Val'] < 0.05), ])  

heatmap(t(expr.AX)[which(tab.AX.BF[,'adj.P.Val'] < 0.05),])      # Display Heatmap


### Benjamini-Hochberg Procedure: FDR=0.1


tab.AX.BH <- topTable(ebayes.AX,                                 # FDR Threshold Toptable
                      coef=2, 
                      number=G, 
                      adjust="BH")

toptable.AX.LIMMA.BH <- data.frame(metabolites = substr(rownames(tab.AX.BH),1,50),
                                   pvalue = tab.AX.BH[,'adj.P.Val'], 
                                   row.names = NULL)


mets.limma.AX <- toptable.AX.LIMMA.BH[1:length(which(tab.AX.BH[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.limma.AX, "mets_limma_AX.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.AX.LIMMA.BH[1:length(which(tab.AX.BH[,'adj.P.Val'] < 0.1)),],
      format='markdown')

heatmap(t(expr.AX)[which(tab.AX.BH[,'adj.P.Val'] < 0.1),])       # Display Heatmap


## FAMT Model Analysis


expr.AX.FAMT <- data.AX                                          # Set Data Matrix

cov.AX.FAMT  <- data.frame(id    = rownames(expr.AX),            # Set Covariates Matrix
                           trmt  = as.factor(trmt.AX), 
                           batch = as.factor(batch.AX))

data.AX.FAMT <- as.FAMTdata(expression = expr.AX.FAMT,           # Make Data Structure
                            covariates = cov.AX.FAMT, 
                            idcovar    = 1)

fit.AX.FAMT  <- modelFAMT(data.AX.FAMT,                          # Fit FAMT Model
                          x    = 2, 
                          test = 2, 
                          nbf  = 2)


toptable.AX.FAMT <- data.frame(metabolites = metabolites[order(fit.AX.FAMT$adjpval)],
                               pvalues     = fit.AX.FAMT$adjpval[order(fit.AX.FAMT$adjpval)],
                               row.names   = NULL)

sig.AX.FAMT      <- length(which(fit.AX.FAMT$adjpval < 0.1))     # Number of Significant Metabolites (BH)

cat(sig.AX.FAMT,                                                 # Display Number of Metabolites (BH)
    ' significant metabolites at FDR = 0.1.')


mets.FAMT.AX <- toptable.AX.FAMT[1:sig.AX.FAMT, 1]

write.table(mets.FAMT.AX, "mets_FAMT_AX.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.AX.FAMT[1:sig.AX.FAMT,],                          # Display Significant Metabolites (BH)
      row.names=NA)


## PCA + LIMMA


batch.AX.PCA <- svd(t(expr.AX) - rowMeans(t(expr.AX)))$v[,1]                       # Batch Estimation

mod.AX.PCA <- model.matrix(~trmt.AX+batch.AX.PCA)                                  # Model Matrix

fit.AX.PCA <- lmFit(t(expr.AX), mod.AX.PCA)                                        # Fit Model

ebayes.AX.PCA <- eBayes(fit.AX.PCA)                                               # Calculate Statistics

tab.AX.PCA <- toptable(fit.AX.PCA, coef=2,                                        # Toptable
                      number=ncol(expr.AX), p.value=0.1)     

cat(toString(nrow(tab.AX.PCA)), "significant metabolites found by PCA.")         # Print Statement

toptable.AX.PCA.BH <- data.frame(metabolites = substr(rownames(tab.AX.PCA),1,50), # Toptable
                                pvalue = tab.AX.PCA[,'adj.P.Val'], 
                                row.names = NULL)

mets.PCA.AX <- toptable.AX.PCA.BH[1:length(which(tab.AX.PCA[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.PCA.AX, "mets_PCA_AX.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.AX.PCA.BH[1:length(which(tab.AX.PCA[,'adj.P.Val'] < 0.1)),],       # Print Statement
      format='markdown')


## SVA + LIMMA


library(sva)

mod  <- model.matrix(~trmt.AX)                                             # Initial Model Matrix

mod0 <- cbind(mod[,1])                                                    # Initial Model Matrix

batch.AX.SVA <- sva(t(expr.AX), mod, mod0)                                  # Batch Estimation

mod.AX.SVA <- cbind(mod, batch.AX.SVA$sv)                                   # SVA Model Matrix

fit.AX.SVA <- lmFit(t(expr.AX), mod.AX.SVA)                                  # Fit Model

ebayes.AX.SVA <- eBayes(fit.AX.SVA)                                         # Calculate Statistics

tab.AX.SVA <- toptable(fit.AX.SVA, coef=2,                                  # Toptable
                      number=ncol(expr.AX), p.value=0.1)      

cat(toString(nrow(tab.AX.SVA)), "significant metabolites found by SVA.")   # Print Statement

toptable.AX.SVA.BH <- data.frame(metabolites = substr(rownames(tab.AX.SVA),1,50), # Toptable
                                pvalue = tab.AX.SVA[,'adj.P.Val'], 
                                row.names = NULL)

mets.SVA.AX <- toptable.AX.SVA.BH[1:length(which(tab.AX.SVA[,'adj.P.Val'] < 0.1)), 1]

write.table(mets.SVA.AX, "mets_SVA_AX.txt", sep="\t", row.names=FALSE, col.names=FALSE)

kable(toptable.AX.SVA.BH[1:length(which(tab.AX.SVA[,'adj.P.Val'] < 0.1)),],       # Print Statement
      format='markdown')


#---------------#
# Venn Diagrams #
#---------------#


discoveries.AX <- list("t-tests"=mets.t.AX, "LIMMA"=mets.limma.AX, "SVA + LIMMA"=mets.SVA.AX, 
                       "PCA + LIMMA"=mets.PCA.AX, "RRmix"=mets.RRmix.AX, "FAMT"=mets.FAMT.AX)

Venn.AX <- Venn(discoveries.AX)

Venn.AX <- compute.Venn(Venn.AX, type="AWFE")

SetLabels.AX <- VennGetSetLabels(Venn.AX)

SetLabels.AX[, "x"] <- c(-3.8, -4, 1, 2, -1.4, -2)
SetLabels.AX[, "y"] <- c(3.6, 4.3, 2, 1.55, 2.7, 2)

Venn.AX <- VennSetSetLabels(Venn.AX, SetLabels.AX)

gp.AX <- VennThemes(Venn.AX)

gp.AX[["SetText"]][["Set1"]]$col <- 1
gp.AX[["SetText"]][["Set2"]]$col <- 2
gp.AX[["SetText"]][["Set3"]]$col <- 6
gp.AX[["SetText"]][["Set4"]]$col <- "darkorange"
gp.AX[["SetText"]][["Set5"]]$col <- 3
gp.AX[["SetText"]][["Set6"]]$col <- 5

gp.AX[["Set"]][["Set1"]]$col <- 1
gp.AX[["Set"]][["Set2"]]$col <- 2
gp.AX[["Set"]][["Set3"]]$col <- 6
gp.AX[["Set"]][["Set4"]]$col <- "darkorange"
gp.AX[["Set"]][["Set5"]]$col <- 3
gp.AX[["Set"]][["Set6"]]$col <- 5

plot(Venn.AX, gp=gp.AX)

pdf("VennAX.pdf")

plot(Venn.AX, gp=gp.AX)

dev.off()


discoveries.AX.4 <- list("t-tests"=mets.t.AX, "LIMMA"=mets.limma.AX, 
                         "RRmix"=mets.RRmix.AX, "FAMT"=mets.FAMT.AX)

Venn.AX.4 <- Venn(discoveries.AX.4)

Venn.AX.4 <- compute.Venn(Venn.AX.4, type="ellipses", doWeights=FALSE)

SetLabels.AX.4 <- VennGetSetLabels(Venn.AX.4)

SetLabels.AX.4[, "x"] <- c(9, -11, 7, -9)
SetLabels.AX.4[, "y"] <- c(6, 6, 8, 8)

Venn.AX.4 <- VennSetSetLabels(Venn.AX.4, SetLabels.AX.4)

gp.AX.4 <- VennThemes(Venn.AX.4)

gp.AX.4[["SetText"]][["Set1"]]$col <- 1
gp.AX.4[["SetText"]][["Set2"]]$col <- 2
gp.AX.4[["SetText"]][["Set3"]]$col <- 3
gp.AX.4[["SetText"]][["Set4"]]$col <- 5

gp.AX.4[["Set"]][["Set1"]]$col <- 1
gp.AX.4[["Set"]][["Set2"]]$col <- 2
gp.AX.4[["Set"]][["Set3"]]$col <- 3
gp.AX.4[["Set"]][["Set4"]]$col <- 5

plot(Venn.AX.4, gp=gp.AX.4)

pdf("VennAX_4.pdf")

plot(Venn.AX.4, gp=gp.AX.4)

dev.off()



#--------------#
# Session Info #
#--------------#


sessionInfo()