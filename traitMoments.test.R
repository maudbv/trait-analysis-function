######################################################################################################
#
# Testing the different MOMENTS of the local community trait distribution 
# against random abundance permutations 
# 
# Null model: within each community independently, observed species abundances are randomly 
# permutated among the species of the community.
# This effectively breaks the link between trait values and abundance values, and thus tests if 
# traits are distributed randomly with respect to abundances within the community.
# 
# We calculate the 4 moments of a trait distribution : mean, variance, skewness, kurtosis
# These four moments are weighted by species abundances within each community.
# We test each observed moment against the null distribution created by *nreps* independent
# random permutations. 
# 
# Significance is tested as the proportion of null values lower than the observed value 
# (i.e. the lower quantile). For instance, significantly lower variance than expected 
# (e.g. in the lower 5th quantile) indicates that abundant trait values tend to converge 
# more than expected within the community,given the local pool of coexisting species
# and the evenness of the abundance distribution.
#
#
#### INPUTS
# comm :    a community matrix with samples (or communities) as rows, and species as columns, 
#           which must include local species abundances (or else the null model is meaningless)
# trait:    the name of the trait under study, has to be the name of the column in df)
# df:       a matrix or data frame with species as rows, and traits as columns, must have species names as rownames.
# nreps :   number of repetitions of the null sampling
#
# #### OUTPUT
#   a data frame listing for each sample/community and each of the four moments: the observed value, 
#   the mean and standard deviation of null values, and the proportion of null values found to be
#   lower than the observed value.
# 
#######################################################################################################
#
# Code written by Maud Bernard-Verdier, as used in the following paper :
# Bernard-Verdier, M., Navas, M., Vellend, M., Violle, C., Fayolle, A., & Garnier, E. (2012).
# Community assembly along a soil depth gradient : contrasting patterns of plant trait convergence 
# and divergence in a Mediterranean rangeland. Journal of Ecology, 100(6), 1422-1433.
# doi:10.1111/1365-2745.12003
#
#  Copyright GPL-2 2012 Maud Bernar-Verdier
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.


######################################################################################################

null.mom=function(comm,trait,df=traits,nreps=9999) {
  
  require(picante)
  
  # trait vector:
  df <- as.data.frame(df)
  df[,trait]=as.numeric( df[,trait])
  t=na.omit(df2vec(df,trait))
  t=t[names(t)%in%names(comm)]
  
  # Community matrix reduced to species for which the trait is known
  comm <- comm[,names(t)]
  comm <- comm/rowSums(comm)  # transform abundance data in relative abundance data
  
  
  # Number of community samples
  n_sites <- nrow(comm)
  
  # Initiating result objects 
  
  # community weighted means
  obs.cwms=NULL
  mean.null.cwms=NULL
  sd.null.cwms=NULL
  P.cwm.lower=NULL
  
  # community weighted variance
  obs.cwvs=NULL
  mean.null.cwvs=NULL
  sd.null.cwvs=NULL
  P.cwv.lower=NULL
  
  # community weighted skewness
  obs.cwss=NULL
  mean.null.cwss=NULL
  sd.null.cwss=NULL
  P.cws.lower=NULL
 
  # community weighted kurtosis
  obs.cwks=NULL
  mean.null.cwks=NULL
  sd.null.cwks=NULL
  P.cwk.lower=NULL
  
  # Looping on the n sites
  for (a in 1:n_sites) {
    com.obs=comm[a,comm[a,]>0]  # We only keep the observed species in community i
    t.obs=t[names(com.obs)]
    
    # Observed moment values in community i :
    obs.cwm=wtd.mean(t.obs,com.obs)   # Weighted mean
    obs.cwv=sum(com.obs[com.obs>0]*(t.obs[com.obs>0]-obs.cwm)^2)  # Weighted variance
    obs.cws=sum(com.obs[com.obs>0]*(t.obs[com.obs>0]-obs.cwm)^3)/obs.cwv^(3/2)    # Weighted skewness
    obs.cwk=sum(com.obs[com.obs>0]*(t.obs[com.obs>0]-obs.cwm)^4)/obs.cwv^2 - 3    # Weighted kurtosis
    
    null_cwm_array=NULL
    null_cwv_array=NULL
    null_cws_array=NULL
    null_cwk_array=NULL
    
    # looping on the nreps null samplings
    for (k in 1:nreps)
    {
      s.com=sample(com.obs) # shuffles abundance values among the observed species in the community
      
      # Null values of the four moments for community i
      null_cwm_array[k] =  wtd.mean(t.obs,s.com)  # Weighted mean
      null_cwv_array[k]=sum(s.com[s.com>0]*(t.obs[s.com>0]-null_cwm_array[k])^2) # Weighted variance
      null_cws_array[k]=sum(s.com[s.com>0]*(t.obs[s.com>0]-null_cwm_array[k])^3)/null_cwv_array[k]^(3/2)  # Weighted skewness
      null_cwk_array[k]=sum(s.com[s.com>0]*(t.obs[s.com>0]-null_cwm_array[k])^4)/null_cwv_array[k]^2 - 3  # Weighted kurtosis
      
    }
    
    # Distribution of null values, to which the observed value is added, for each of the four moments :
    null_cwm_array=c(null_cwm_array,obs.cwm)
    null_cwv_array=c(null_cwv_array,obs.cwv)
    null_cwk_array=c(null_cwk_array,obs.cwk)
    null_cws_array=c(null_cws_array,obs.cws)
    
    ######## Test of the observed value against the null distribution
    
    # Community weighted mean :
    
    # Proportion of null values lower than the observed value (one sided test)
    P.cwm.lower[a]=(sum(as.numeric(null_cwm_array<obs.cwm))+sum(as.numeric(null_cwm_array==obs.cwm))/2)/(nreps+1)
  
    # additional results:
    obs.cwms[a]=obs.cwm
    mean.null.cwms[a]= mean(null_cwm_array)     
    sd.null.cwms[a]= sd(null_cwm_array)
    
    # CW variance
    P.cwv.lower[a]=(sum(as.numeric(null_cwv_array<obs.cwv))+sum(as.numeric(null_cwv_array==obs.cwv))/2)/(nreps+1)
    obs.cwvs[a]=obs.cwv
    mean.null.cwvs[a]= mean(null_cwv_array)     
    sd.null.cwvs[a]= sd(null_cwv_array)
    
    # CW Kurtosis
    P.cwk.lower[a]=(sum(as.numeric(null_cwk_array<obs.cwk))+sum(as.numeric(null_cwk_array==obs.cwk))/2)/(nreps+1)
    obs.cwks[a]=obs.cwk
    mean.null.cwks[a]= mean(null_cwk_array)     
    sd.null.cwks[a]= sd(null_cwk_array)
    
    # CW skewness
    P.cws.lower[a]=(sum(as.numeric(null_cws_array<obs.cws))+sum(as.numeric(null_cws_array==obs.cws))/2)/(nreps+1)
    obs.cwss[a]=obs.cws
    mean.null.cwss[a]= mean(null_cws_array)     
    sd.null.cwss[a]= sd(null_cws_array)
    
  }
  
  # results data frame
  results=data.frame(round(obs.cwms,3),round(mean.null.cwms,3),round(sd.null.cwms,3),round(P.cwm.lower,5),
                     round(obs.cwvs,3),round(mean.null.cwvs,3),round(sd.null.cwvs,3),round(P.cwv.lower,5),
                     round(obs.cwks,3),round(mean.null.cwks,3),round(sd.null.cwks,3),round(P.cwk.lower,5),
                     round(obs.cwss,3),round(mean.null.cwss,3),round(sd.null.cwss,3),round(P.cws.lower,5))
  
  names(results)=c("obs.cwm","mean.null.cwm","sd.null.cwm","P.lower.cwm",
                   "obs.cwv","mean.null.cwv","sd.null.cwv","P.lower.cwv",
                   "obs.cwk","mean.null.cwk","sd.null.cwk","P.lower.cwk",
                   "obs.cws","mean.null.cws","sd.null.cws","P.lower.cws")
  
  row.names(results)=row.names(comm)
  return(results)
}


