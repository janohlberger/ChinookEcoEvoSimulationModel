##########################################################################
##                                                                      ##
## Parameters for simulating trend in Chinook salmon age-size structure ##
##                                                                      ##
##########################################################################
parameters<-c("s",s<-1,  					## scenario number 
"nR",nR<-1, 								## number of stochastic runs 
"bin.size",bin.size<-10,					## width size and trait bins
##------------------------------------------## time periods
"nY.equi",nY.equi<-100,		 				## equilibrate, only evolution
"nY.trans",nY.trans<-100,					## transition periods (1&2)
"nY.obs",nY.obs<-50, 						## 'observation' period 
##------------------------------------------## von Bertalanffy growth
"Linf",Linf<-1200,						 	## average asymptotic length
"k",k<-0.3,									## growth rate
"Linf.vs.k",Linf.vs.k<--0.0005, 			## slope of Linf vs k
"sd.g",sd.g<-0.2,							## variance in growth
"sd.L",sd.L<-120,							## variance in Linf 
"SaS",SaS<-100,								## size at smolt stage
##------------------------------------------## correlated traits
"Linf.vs.pmrn",Linf.vs.pmrn<-0.3,			## slope Linf-vs-PMRN intercept
##------------------------------------------## fecundity (egg production)
"a_fec",a_fec<-0.3,  						## fecundity constant
"b_fec",b_fec<-1.5, 						## fecundity exponent
##------------------------------------------## recruitment
"a_rec",a_rec<-0.1,							## egg-to-smolt productivity
"b_rec",b_rec<-250000,						## smolt capacity
"sd.rec",sd.rec<-0.2,						## variance in recruitment
##------------------------------------------## background survival
"surv_z",surv_z<-0.1,						## survival to ocean-1
"surv_a",surv_a<-0.9,						## annual survival ocean 1-4
##------------------------------------------## maturation reaction norm
"a_mat",a_mat<-850,							## PMRN intercept
"b_mat",b_mat<--50,							## PMRN slope
"sd_mat",sd_mat<-100,						## PMRN standard deviation
##------------------------------------------## predation removals
"a_pred",a_pred<-800,						## predation size at midpoint 
"b_pred",b_pred<-0.025,						## predation steepness
"PR_max",PR_max<-0.4,						## maximum predation rate
##------------------------------------------## harvest removals
"a_harv",a_harv<-725,						## size at maximum selectivity
"sd_harv",sd_harv<-0.5,  					## harvest selectivity variance
"HR_max",HR_max<-0.6)						## maximum harvest rate
##########################################################################