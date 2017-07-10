# ########################################################################################
# FUNCTION to test alphaRao, defined as the functional or phylogenetic divergence 
# of dominant values within communities.This function calculates Rao's quadratic entropy (Rao 1982)
# within communities using functional trait distances (alphaFD) or phylogenetic distances (alphaPD).
# Observed values are compared against a null model shuffling abundances among species wthin communities. 
# Both FD and PD can be calculated at the same time if trait data and a phylogeny are provided.
# 
# INPUTS :
#   X :     a numeric community data matrix with samples (or sites) as rows and species as columns
#           must include numeric species abundances within sites for the null model to make sense
# 
#   df :    (only for FD) a data frame containing trait data (as numeric columns) for the species listed in X
# 
#   var :   (only for FD) a character vector of column names corresponding to the traits contained in df
#           If only one trait nae is provided, then a univariate index of functional diversity 
#           will be calculated.
#           If a vector of several names is provided, then a multivariate index will be tested.
#   
#    phy :  (only for PD) object of class phylo, a phylogenetic tree of the species contained in the community matrix
#   
#   nreps : numeric, the number of repetition of the null model
#   
#   Jost :  boolean. Optional correction on the functional distance matrix, 
#           according to recommendations by Lou Jost (Ecology, 10:88, 2007), 
#           and as implemented by De Bello et al.(Journal of Vegetation Science, 21:5, 2010).
#   
# 
# OUTPUTS :
#   a data frame with community sites/samples as rows, and diversity results as columns :
#   - "obs.FD" :         the observed Rao alpha diversity index in each site 
#   - "mean.null.FD" :   the mean of the distribution of null Rao values generated independently for each site  
#   - "sd.null.FD" :     the standard deviation of the null rao distribution for each site
#   - "prop.lower.FD" :  the proportion of null rao values which were lower than the observed Rao
#   - "alphaFD.ES" :     the calculated effect size for the alpha Rao, as described in Bernard-Verdier et al. 2013. 
# 
#   + same columns for PD if a phylogeny is provided.
#
# NOTE : the more classic "standardized effect sizes" defined as (observed-mean(null))/sd(null)
# can also be calculated from the output if needed.
# 
# This Rcode was written by Maud Bernard-Verdier, as used in :
# Bernard-Verdier, M., Flores, O., Navas, M.-L., & Garnier, E. (2013). Partitioning phylogenetic and functional diversity into alpha and beta 
# components along an environmental gradient in a Mediterranean rangeland. (F. de Bello, Ed.)Journal of Vegetation Science, 24(5), 877-889. doi:10.1111/jvs.12048
#
# #######################################################################################################

alphaRao.test=function(X, df, var, phy, Jost=T, nreps=9999) {

  require (ape)
  require(picante)

  # For functional diversity
  if (!is.null(df)) {
   
    # if no trait names are provided in 'var', then the whole df is used
    if (is.null(var)) {
      t <- na.omit(df)
      t[row.names(t)%in%names(X),]
      if (is.null(t))       stop("error : trait values are missing for all species")
      if (is.numeric(t)==F) stop("error : non numeric trait data")
      X <- X[,rownames(t)]              # matching the order of species names with the trait data
      X <- X/rowSums(X)                 # re-calculating relative abundances
    }
    
    # if only one trait:
    if (length(var)==1) {
      t <- na.omit(df2vec(df,var)) 
      if (is.null(t)) stop("error : trait values are missing for all species")
      t <- t[names(t)%in%names(X)]  # keeping only species names included in the community matrix
      X <- X[,names(t)]             # matching the order of species names with the trait data
      X <- X/rowSums(X)             # re-calculating relative abundances
    }
  
    # if more than one trait :
    if (length(var)>1) {
      t <- na.omit(df[,var]) # selects the trait data and only species without missing values
      if (is.null(t)) stop("error : some trait values are missing for all species")
      
      t <- t[rownames(t)%in%names(X),]  # keeping only species names included in the community matrix
      X <- X[,rownames(t)]              # matching the order of species names with the trait data
      X <- X/rowSums(X)                 # re-calculating relative abundances
    }  
  
    # Matrix of functional distances among species :
    dfunc <- dist(t, method="euclidean")
    # OPTIONAL :
    if (Jost) dfunc=dfunc/max(dfunc, na.rm=T)  
    dfunc=as.matrix(dfunc)
  } 
  
######## For phylogenetic diversity:
  if (!is.null(phy)) {
    
    # removing any species name in the phylogeny which is absent from the community data matrix :
    if (length(which(!phy$tip.label%in%names(X)))!=0)  { 
      phy=drop.tip(phy,phy$tip.label[!phy$tip.label%in%names(X)],trim.internal=T)
    }
    
      X <- X[,phy$tip.label]
      X <- X/rowSums(X)
   
  # phylogenetic distance matrix
  dphyl=as.dist(cophenetic(phylo))/2
  # OPTIONAL:
  if (Jost) dphyl=dphyl/max(dphyl, na.rm=T)    
  dphyl=as.matrix(dphyl) # transform it in a matrix for easier manipulation
  dphyl=dphyl[names(X),names(X)] 
  }
  
########## Final matching of species names if both phylogenetic and functional informations are provided :
  if( !is.null(phy) & !is.null(df)) {
    dphyl=dphyl[names(X),names(X)] 
    dfunc=dfunc[names(X),names(X)] 
  }
  
   
###### Initiating result objects
  FD=NULL
  obs.FDs <- NULL
  mean.null.FDs <- NULL
  sd.null.FDs <- NULL
  prop.lower.FD <- NULL
  alphaFD.ES <- NULL
  
  PD=NULL
  obs.PDs <- NULL
  mean.null.PDs <- NULL
  sd.null.PDs <- NULL
  prop.lower.PD <- NULL
  alphaPD.ES <- NULL

###### Looping on the n sites
  for (a in 1:nrow(X)) {
    com.obs <- X[a,X[a,]>0]
    
    if (!is.null(df)) {
    Fdist.obs <- as.matrix(dfunc)[colnames(com.obs),colnames(com.obs)]  
    obs.FD <- (as.matrix(com.obs) %*% Fdist.obs %*% as.matrix(t(com.obs)))[1,1]
    if (Jost) obs.FD=1/(1-obs.FD) # OPTIONAL Jost correction
    null_FD_array <- NULL
    }
    
    if (!is.null(phy)) {
    Pdist.obs <- as.matrix(dphyl)[colnames(com.obs),colnames(com.obs)]  
    obs.PD <- (as.matrix(com.obs) %*% Pdist.obs %*% as.matrix(t(com.obs)))[1,1]
    # OPTIONAL :
    if (Jost) obs.PD=1/(1-obs.PD)
    null_PD_array <- NULL
    }
    
    
    # looping on the nreps null samplings
      for (k in 1:nreps) {
        s.com <- sample(com.obs) # shuffling abundances within the community species pool
        
        if (!is.null(df)) {
          null_FD_array[k] <- (as.matrix(s.com) %*% Fdist.obs %*% as.matrix(t(s.com)))[1,1]
          if (Jost) null_FD_array[k]=1/(1-null_FD_array[k])        # OPTIONAL Jost correction
        }
        
        if (!is.null(phy)) {
          null_PD_array[k] <- (as.matrix(s.com) %*% Pdist.obs %*% as.matrix(t(s.com)))[1,1]
          if (Jost) null_PD_array[k]=1/(1-null_PD_array[k])        # OPTIONAL Jost correction
        }
      }
    
######### Comparing to null model :
    
   if (!is.null(df)) {
    null_FD_array <- c(null_FD_array,obs.FD) # the null distribution of FD values
    
    # Proportion of null values which are lower than the observed FD (= the corresponding quntile of the null distribution)
    prop.lower.FD[a] <- (sum(as.numeric(null_FD_array<obs.FD))+sum(as.numeric(null_FD_array==obs.FD))/2)/(nreps+1)
    
    # calculating the effect size (ranges between -1 and 1) as described in 
    # Bernard-Verdier et al. (Journal of Ecology, 2012)
    alphaFD.ES[a] <- 2*(prop.lower.FD[a] - 0.5)
    
    # saving additional results for site a :
    obs.FDs[a] <- obs.FD
    mean.null.FDs[a] <- mean(null_FD_array)     
    sd.null.FDs[a] <- sd(null_FD_array)
     }
        
  if (!is.null(phy)) {
    null_PD_array <- c(null_PD_array,obs.PD) # the null distribution of PD values
    
    # Proportion of null values which are lower than the observed PD (= the corresponding quntile of the null distribution)
    prop.lower.PD[a] <- (sum(as.numeric(null_PD_array<obs.PD))+sum(as.numeric(null_PD_array==obs.PD))/2)/(nreps+1)
    
    # calculating the effect size (ranges between -1 and 1) as described in 
    # Bernard-Verdier et al. (Journal of Ecology, 2012)
    alphaPD.ES[a] <- 2*(prop.lower.PD[a] - 0.5)
    
    # saving additional results for site a :
    obs.PDs[a] <- obs.PD
    mean.null.PDs[a] <- mean(null_PD_array)     
    sd.null.PDs[a] <- sd(null_PD_array)
 
  }
}      
  
  # results data frame :
 if (!is.null(df)) {
   FD <- data.frame(obs.FD=obs.FDs,mean.null.FD=mean.null.FDs,sd.null.FD=sd.null.FDs,
                   prop.lower.FD=prop.lower.FD,alphaFD.ES=alphaFD.ES)
   output=FD
 }
 if (!is.null(phy)) {
   PD <- data.frame(obs.PD=obs.PDs,mean.null.PD=mean.null.PDs,sd.null.PD=sd.null.PDs,
                   prop.lower.PD=prop.lower.PD,alphaPD.ES=alphaPD.ES)
   output=PD
 }
 if( !is.null(phy) & !is.null(df))  output=cbind(FD,PD)
  row.names(output) <- row.names(X)
  return(output)
}

