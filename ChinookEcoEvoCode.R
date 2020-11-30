##########################################################################
##                                                                      ##
##   Simualtion model of changes in Chinook salmon age-size structure   ##
##                                                                      ##
##########################################################################
ChinookEcoEvoModel<-function(parameters) {
  
  ############################################### variables and functions
  
  ############################################## years and ages
  nY<-nY.equi+nY.trans+nY.obs					## total number years 
  years<-seq(nY)  							## index of all years
  ages<-seq(5)								## age-classes (ocean age)
  nA<-length(ages)							## number of ages 
  
  ############################################## discretized sizes
  s.breaks<-seq(0,1500,bin.size) 				## vector of size bins
  nS<-length(s.breaks)-1						## number of bins
  size.bins<-(s.breaks+bin.size/2)[-(nS+1)] 	## size mid-points
  all.size.bins<-data.frame(size.bin=size.bins)
  ##---------------------------------------## function to make size bins
  SBIN<-function(size,bin=bin.size) { bin*floor(size/bin)+bin/2 }
  
  ############################################ fecundity (egg production)
  FEC<-function(x,a=a_fec,b=b_fec) { a*x^b }
  
  ########################################### recruitment (1:1 sex ratio)
  REC<-function(sizes,alpha=a_rec,beta=b_rec) { sum(FEC(sizes)*0.5)/((1/alpha)+(1/beta)*sum(FEC(sizes)*0.5)) } 
  
  ########################################### vonBertalanffy growth model
  ##---------------------------------------## based on age for first year
  sVBGM_a<-function(age,Linf_vB=Linf,k_vB=k,sd_vB=sd.g) { Linf_vB*(1-exp(-k_vB*age))*exp(rnorm(length(Linf_vB),-0.5*sd.g^2,sd.g)) }
  ##---------------------------------------## based on size (stochastic)
  sVBGM_s<-function(size,Linf_vB=Linf,k_vB=k,sd_vB=sd.g) { size+k_vB*(Linf_vB-size)*exp(rnorm(length(Linf_vB),-0.5*sd.g^2,sd.g)) } 
  
  ########################################### trade-off growth parameters
  a_vonB<-k-Linf*Linf.vs.k					## intercept
  TOFF<-function(L.inf,a=a_vonB,b=Linf.vs.k) { a+L.inf*b } 	
  
  ############################################## discretized traits
  t.breaks<-seq(500,2500,bin.size) 	 		## vector of trait bins
  nT<-length(t.breaks)-1						## number of bins 
  trait.bins<-(t.breaks+bin.size/2)[-(nT+1)]	## trait mid-points
  all.trait.bins<-data.frame(trait.bin=trait.bins)
  ##---------------------------------------## function to make trait bins
  TBIN<-function(trait,bin=bin.size) { bin*floor(trait/bin)+bin/2 }
  
  ############################################## maturation reaction norm
  PMRN<-function(age,size,a=a_mat,b=b_mat,sd=sd_mat) { 1/(1+exp(-(size-(a+b*age))/sd)) }
  
  ########################################## correlated maturation trait
  a_tcor<-a_mat-Linf*Linf.vs.pmrn             ## PMRN intercept
  AMAT<-function(L.inf,a=a_tcor,b=Linf.vs.pmrn) { a+L.inf*b } 
  
  ############################################## predation removals
  ##------------------------------------------## predation selectivity 
  PRED<-function(size,a=a_pred,b=b_pred) { 1/(1+exp(-b*(size-a))) }
  ##------------------------------------------## annual predation rate
  PR_y<-c(rep(0,nY.equi),seq(0,PR_max*0.25,length=nY.trans*0.5), rep(PR_max*0.25,nY.trans*0.5),seq(PR_max*0.25,PR_max,length=nY.obs)) 
  
  ############################################## harvest removals
  ##------------------------------------------## harvest selectivity 
  HARV<-function(size,a=a_harv,sd=sd_harv) { exp(-(log(size)-log(a))^2/(2*sd^2)) }
  ##------------------------------------------## annual harvest rate
  HR_y<-c(rep(0,nY.equi),seq(0,HR_max,length=nY.trans*0.5), rep(HR_max,nY.trans*0.5),seq(HR_max,HR_max*0.5,length=nY.obs))		
  
  ################################################### arrays for results
  runs<-seq(nR)
  adims<-c(nY,nA,nS,nT,nR)						
  dnames<-list(years,ages,size.bins,trait.bins,runs) 
  IMM_yastr<-array(dim=adims,dimnames=dnames) ## immature
  MAT_yastr<-array(dim=adims,dimnames=dnames) ## mature
  PRED_yastr<-array(dim=adims,dimnames=dnames)## predation
  RET_yastr<-array(dim=adims,dimnames=dnames)	## return
  HARV_yastr<-array(dim=adims,dimnames=dnames)## harvested
  ESC_yastr<-array(dim=adims,dimnames=dnames)	## escapement
  
  ###################################################### start simulation
  ZERO<-function(x) { (abs(x)+x)/2 } ## zero negative values
  
  ################################################## loop stochastic runs
  for(m in 1:nR) {
    print(paste0("scenario=",s,"_run=",m)) ## print scenario and run 
    set.seed(m) ## set only when used for reproducibility 
    
    ######################################################### loop years
    PE.rec<-rnorm(nY,0,sd.rec) ## lognormal recruitment error
    for (y in 1:nY) {
      print(paste0("year=",y)) ## print year
      
      ########################################################### year-1
      if(y==1) {
        N_a_y1<-round(b_rec*surv_z*(c(1,0.75,0.5,0.15,0.02))) 
        totNs_y1<-sum(N_a_y1);cumNs_y1<-cumsum(N_a_y1)
        ##-------------------------------------------------## all traits 
        Linfs<-rtnorm(totNs_y1,Linf,sd.L,lower=0) 
        ks<-TOFF(Linfs)
        ##+++++++++++++++++++++++++++++++++++++++++++++++++++## age loop
        for (a in 1:nA) { 
          if(a==1) index<-1:cumNs_y1[a]
          if(a!=1) index<-(cumNs_y1[a-1]+1):cumNs_y1[a]
          ##===================================================## growth
          SaSs<-rlnorm(length(Linf[index]),log(SaS),0.1)
          age.tab<-data.frame(cbind(size.bin=SBIN(SaSs+sVBGM_a(a,Linfs[index], ks[index])),trait.bin=TBIN(Linfs[index]))) ## based on age in 1st year
          ##-------------------------------------## trait and size bins
          counts<-ddply(age.tab,.(size.bin=age.tab$size.bin, trait.bin=age.tab$trait.bin),nrow)
          new.tab<-reshape(counts,idvar="size.bin",timevar="trait.bin", direction="wide")
          ##---------------------------------## expand to all size bins
          tab.all<-merge(all.size.bins,new.tab,by="size.bin",all.x=T)
          tab.all<-data.frame(t(tab.all))
          names(tab.all)<-tab.all[1,]
          tab.all<-tab.all[-1,]
          ##--------------------------------## expand to all trait bins
          tab.all$trait.bin<-gsub("V1.","",rownames(tab.all))
          tab.all<-merge(all.trait.bins,tab.all,by="trait.bin",all.x=T)
          rownames(tab.all)<-tab.all[,1]
          tab.all<-t(tab.all[,-1])
          tab.all[is.na(tab.all)]<-0	
          ##===============================================## maturation
          grid<-expand.grid(size.bin=size.bins,trait.bin=trait.bins) 
          grid$PMRN<-round(PMRN(age=a,size=grid$size.bin,a=AMAT(grid$trait.bin)),4)
          probM_ast<-reshape(grid,idvar="size.bin",timevar="trait.bin", direction="wide")[,-1]
          rownames(probM_ast)<-size.bins;names(probM_ast)<-trait.bins
          MAT_yastr[y,a,,,m]<-as.matrix(round(tab.all*probM_ast))
          IMM_yastr[y,a,,,m]<-round(tab.all-MAT_yastr[y,a,,,m])		
        } ##++++++++++++++++++++++++++++++++++++++++++++## end age loop
      } ##------------------------------------## end statement if(y==1)
      
      ##################################################### years 2:nY 
      if(y!=1) {	
        ##++++++++++++++++++++++++++++++++++++++++++++++++++## age loop
        for (a in 1:5) {
          if(a==1) {
            ##===================================## number of offspring
            N.spawn_st<-apply(ESC_yastr[y-1,,,,m],c(2,3),sum) ##spawners by trait/size
            N.eggs_st<-round(N.spawn_st*FEC(size.bins)) ## times fecundity
            N.eggs_t<-apply(N.eggs_st,2,sum) ## number eggs by trait value
            Egg.traits<-as.numeric(rep(names(N.eggs_t),N.eggs_t))
            N.spawn_s<-rowSums(N.spawn_st) ## number of spawners by size
            Spawner.sizes<-as.numeric(rep(names(N.spawn_s),N.spawn_s)) ##spawner sizes
            ##-------------------------------## density-dependent # smolts
            N.smolt.det<-round(REC(Spawner.sizes)) 
            N.smolt<-round(N.smolt.det*exp(PE.rec[y]))	
            ##--------------------------------------## number of ocean-1s
            N1<-round(N.smolt*surv_z) ## independent of size, trait, and density
            ##===========================## offspring trait distribution 
            ## mean parental trait value and fixed variance
            Linfs.off<-rtnorm(N1,mean(Egg.traits),sd.L,lower=0)
            ks.off<-TOFF(Linfs.off)
            ##==================================================## growth
            ##-------------------------------------------## size-at-age-1
            SaSs<-rlnorm(length(Linfs.off),log(SaS),0.1)
            SaA1<-sVBGM_s(SaSs,Linfs.off,ks.off)
            ##------------------------------## make size and trait matrix
            age.mat<-data.frame(cbind(size.bin=SBIN(SaA1),trait.bin=TBIN(Linfs.off)))
          } ##---------------------------------## end if(a==1) statement
          ##++++++++++++++++++++++++++++++++++++++++++++++++++## ages 2-4 
          if(a!=1) {
            ##------------------------------## immature previous year/age
            if(a %in% c(2,3,4)) { prev.tab<-data.frame(IMM_yastr[y-1,a-1,,,m]) }
            ##-----------------------------## plus group (previous 4s+5s)
            if(a==5) { prev.tab<-data.frame(IMM_yastr[y-1,a-1,,,m]+IMM_yastr[y-1,a,,,m]) }
            ##--------------------------------## long format for expansion
            names(prev.tab)<-as.numeric(as.character(gsub("X","",names(prev.tab))))
            prev.long<-reshape(prev.tab,idvar="size.bin",ids=row.names(prev.tab), times=names(prev.tab),timevar="trait.bin",varying=list(names(prev.tab)), direction="long");names(prev.long)[2]<-"count"
            ##------------------------------## all sizes and traits by bin
            prev.all<-data.frame(cbind(size.bin=as.numeric(rep(prev.long$size.bin,prev.long$count)),trait.bin=as.numeric(rep(prev.long$trait.bin,prev.long$count))))
            ##----------------------------------## calculate Linfs and ks
            Linfs.prev<-prev.all$trait.bin
            ks.prev<-TOFF(Linfs.prev)
            ##-------------------------------## matrix with new size bins
            age.mat<-data.frame(cbind(size.bin=SBIN(sVBGM_s(prev.all$size.bin, Linfs.prev,ks.prev)),trait.bin=Linfs.prev))
          } ##--------------------------------## end if(a!=1) statement
          ##--------------------------------------## trait and size bins
          counts<-ddply(age.mat,.(size.bin=age.mat$size.bin, trait.bin=age.mat$trait.bin),nrow)
          ##----------------------------------## reshape to wide format
          new.tab<-reshape(counts,idvar="size.bin",timevar="trait.bin", direction="wide")
          ##----------------------------------## expand to all size bins
          tab.all<-merge(all.size.bins,new.tab,by="size.bin",all.x=T)
          tab.all<-data.frame(t(tab.all))
          names(tab.all)<-tab.all[1,];tab.all<-tab.all[-1,]
          ##----------------------------------## expand to all trait bins
          tab.all$trait.bin<-gsub("V1.","",rownames(tab.all))
          tab.all<-merge(all.trait.bins,tab.all,by="trait.bin",all.x=T)
          rownames(tab.all)<-tab.all[,1];tab.all<-tab.all[,-1]
          ##-------------------------------## NAs to zeros and transpose
          tab.all[is.na(tab.all)]<-0
          tab.all<-t(tab.all)
          ##--------------------------------------## background survival
          tab.all<-round(tab.all*surv_a)		
          ##===============================================## maturation
          grid<-expand.grid(size.bin=size.bins,trait.bin=trait.bins) 
          grid$PMRN<-round(PMRN(age=a,size=grid$size.bin,a=AMAT(grid$trait.bin)),4)
          probM_ast<-reshape(grid,idvar="size.bin",timevar="trait.bin", direction="wide")[,-1]
          rownames(probM_ast)<-size.bins;names(probM_ast)<-trait.bins
          MAT_yastr[y,a,,,m]<-as.matrix(round(tab.all*probM_ast))
          IMM_yastr[y,a,,,m]<-round(tab.all-MAT_yastr[y,a,,,m])
        } ##+++++++++++++++++++++++++++++++++++++++++++++## end age loop
      } ##-------------------------------------## end statement if(y!=1)
      
      ##====================================================## predation
      counts_ast<-apply(MAT_yastr[y,,,,m],c(1,2,3),sum)
      props_ast<-prop.table(counts_ast,2) ## props mature by age-trait each size
      props_ast[is.na(props_ast)]<-0
      props_s<-apply(MAT_yastr[y,,,,m],2,sum)*PRED(size.bins)/sum(apply(MAT_yastr[y,,,,m],2,sum)*PRED(size.bins)) ## props by size 
      pred_s<-props_s*PR_y[y]*sum(MAT_yastr[y,,,,m]) ## numbers by size to remove
      pred_ast<-array(dim=c(nA,nS,nT)) ## numbers by age-size-trait to remove
      for(bin in 1:nS) { pred_ast[,bin,]<-pred_s[bin]*props_ast[,bin,] }
      pred_ast[is.na(pred_ast)]<-0
      PRED_yastr[y,,,,m]<-apply(pred_ast,c(1,2,3),function(x) floor(x)+ rbinom(n=1,size=1,prob=x-floor(x))) ## binomial draw for fractions of fish
      ##--------------------------------------## mature-predation=return
      RET_yastr[y,,,,m]<-round(MAT_yastr[y,,,,m]-PRED_yastr[y,,,,m])
      RET_yastr[y,,,,m]<-ZERO(RET_yastr[y,,,,m]) ## no negative values allowed
      PRED_yastr[y,,,,m]<-MAT_yastr[y,,,,m]-RET_yastr[y,,,,m] ## actual predation
      
      ##====================================================## harvesting 
      counts_ast<-apply(RET_yastr[y,,,,m],c(1,2,3),sum)
      props_ast<-prop.table(counts_ast,2) ## props return by age-trait each size
      props_ast[is.na(props_ast)]<-0
      props_s<-apply(RET_yastr[y,,,,m],2,sum)*HARV(size.bins)/sum(apply(RET_yastr[y,,,,m],2,sum)*HARV(size.bins)) ## props by size
      harv_s<-props_s*HR_y[y]*sum(RET_yastr[y,,,,m]) ## numbers by size to remove
      harv_ast<-array(dim=c(nA,nS,nT)) ## numbers by age-size-trait to remove
      for(bin in 1:nS) { harv_ast[,bin,]<-harv_s[bin]*props_ast[,bin,] }
      harv_ast[is.na(harv_ast)]<-0
      HARV_yastr[y,,,,m]<-apply(harv_ast,c(1,2,3),function(x) floor(x)+ rbinom(n=1,size=1,prob=x-floor(x))) ## binomial draw for fractions of fish
      ##-------------------------------------## return-harvest=escapement
      ESC_yastr[y,,,,m]<-round(RET_yastr[y,,,,m]-HARV_yastr[y,,,,m])	
      ESC_yastr[y,,,,m]<-ZERO(ESC_yastr[y,,,,m]) ## no negative values allowed
      HARV_yastr[y,,,,m]<-RET_yastr[y,,,,m]-ESC_yastr[y,,,,m] ## actual harvest
      
    } ################################################### end year loop
    
  } ##################################### end loop over stochastic runs
  
  save(list=ls(),file="Rworkspace.RData") ## save output
  
} ###################################################### close function

#######################################################################