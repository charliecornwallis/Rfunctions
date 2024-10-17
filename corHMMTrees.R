#******************************************************************************************
#HMM models run on list of trees
#******************************************************************************************
hmm_trees<-function(trees,data,rate_cats=1,node.states="joint"){
  results <- vector("list", length(trees))
  
  for(i in 1:length(trees)) {
    mod<-corHMM(phy = trees[[i]], data = data, rate.cat = rate_cats,node.states=node.states)
    results[[i]]<-mod
  }
  
  return(results) 
} 
