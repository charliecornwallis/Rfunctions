# prior=list(R=list(R1=list(V=diag(6), nu=5.002)),
#            G=list(G1=list(V=diag(6), nu=5.002),
#                   G2=list(V=diag(6), nu=5.002),
#                   G3=list(V=diag(6), nu=5.002)))
# treetoburn=0
# data=dat
# samples=10
# trees=trees
# tiplabel=dat$animal
# fixed=as.formula(cbind(Ztemp,Zsqrt_precip,Zsqrt_temp_sd,Zsqrt_precip_sd,Zsqrt_temp_sd_byrs,Zsqrt_precip_sd_byrs)~
#                    at.level(fam_non,'Pair'):trait+
#                    at.level(fam_non,'Family'):trait+
#                    at.level(fam_non,'Nonfamily'):trait-1)
# random=as.formula(~ us(at.level(fam_non,'Pair'):trait):animal+
#                     us(at.level(fam_non,'Family'):trait):animal+
#                     us(at.level(fam_non,'Nonfamily'):trait):animal)
# residual=as.formula(~us(trait):units)
# Gnumber=3
# Rnumber=1
# Geffects=c((6*6),(6*6),(6*6))
# Reffects=(6*6)
# family=c("gaussian","gaussian","gaussian","gaussian","gaussian","gaussian")
# nitts_tree=1
# thins_tree=1
# burns_tree=0

#***************************************
#Function for sampling of across trees using MCMCglmm models#
#***************************************

MCMCglmmTrees<-function(prior, treetoburn,data, mev=0,trees,tiplabel,fixed,random,residual,Gnumber,Rnumber,Geffects,Reffects,family="gaussian",samples=1000,nitts_tree,thins_tree,burns_tree,pr=FALSE){
# Explanation
# prior #model prior
# treetoburn #number of trees to discard before saving iterations
# data #data to be used
# mev=0 #is there a mev effect in the model mev=0 (no), mev=1 (yes)
# trees #trees to sample across
# tiplabel #name of tiplabel in dataset
# fixed #fixed effect formula
# random #G side random effects formula
# residual #R side random effects formula
# Gnumber #number of G side random effects (Number of G in prior list)
# Rnumber #number of R side random effects (Number of R in prior list)
# Geffects #number of parameters per G term (e.g. 2x2 covar per G = 4) given as a vector [e.g c(1,4,4) ] with length = Gnumber [e.g. 3: G1 = single variance, G2=2x2, G3=2x2]
# Reffects #
# family="gaussian" # vector of families e.g. c("gaussian","gaussian") for 2 trait multi-response model
# samples=1000
# nitts_tree
# thins_tree
# burns_tree

#1. Create object to write model output to ----
  tree1<-inverseA(trees[[1]])$Ainv
  data$tiplabel<-tiplabel
  
  model<-MCMCglmm(fixed=fixed, random=random, rcov = residual,
                  ginverse=list(tiplabel=tree1),family =family, data = data,prior=prior, 
                  nitt=samples, burnin=0, thin=1,pr=pr,verbose = F,pl=T)
  modeltmp<-model
  
  #2. Run model on each tree. Estimates of last tree provide start values for next tree. Use first 500 trees as a burnin. ----
  trees<-trees[1:(samples+treetoburn)]#make sure list of trees is not longer than burnin + samples to be saved
  
  for(i in 1:length(trees)){
    i=1
    tree<-trees[[i]]
    invtree<-inverseA(tree)$Ainv
    
    #Setup starting values
    #**Liability design ----
    liab_start <- list(Liab=modeltmp$Liab[1,])
    
    #**G design ----
    G_start=list()
    for(j in 1:Gnumber){
      k=j-1
      a=ifelse(j==1,1,sum(Geffects[1:k]))
      b=sum(Geffects[1:j])
      G=list(modeltmp$VCV[1,a:b])
      names(G)<-paste("G",j,sep="")
      G_start=c(G_start,G)
    }
    names(G_start)<-"G"
    
    #**R design ----
    R_start=list()
    for(j in 1:Rnumber){
      j=1
      k=j-1
      glength=length(G_start$G)
      a=ifelse(j==1,glength+1,(glength+sum(Reffects[1:k])))
      b=glength+sum(Reffects[1:j])
      R=list(modeltmp$VCV[1,a:b])
      names(R)<-paste("R",j,sep="")
      R_start=c(R_start,R)
    }
    names(R_start)<-"R"
    
    #combine
    start=list(liab_start,G_start,R_start)
    
    #3. Sample across trees ----
    modeltmp<-MCMCglmm(fixed=fixed, random=random, rcov = residual,
                       ginverse=list(tiplabel=invtree),family =family, data = data,prior=prior, 
                       nitt=nitts_tree, thin=thins_tree,burnin=burns_tree, start=start, pr=pr,verbose = F,pl=T)
    #saving after tree burnin
    if(i>treetoburn){
      model$VCV[i-treetoburn,]<-modeltmp$VCV[1,]
      model$Sol[i-treetoburn,]<-modeltmp$Sol[1,]
      model$Liab[i-treetoburn,]<-modeltmp$Liab[1,]
    }
    #print(i)
  }
  return(model)
}    
