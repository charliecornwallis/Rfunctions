#******************************************************************************************
#MCMCglmm processing to extract values from models for ancestral nodes - can take list of trees (first part of function) or single trees (second part of function)
#******************************************************************************************

mcmcglmm_pred_states<-function(mr="No",trees,phy_name="animal",model,dat,trait1,trait2=NULL,binomial="No",state1=NULL,state2=NULL,cutoff=0.5){
  #multiple responses or not
  trait_sub<-ifelse(mr=="No",paste(phy_name,".",sep=""),paste0("trait",trait1,".",paste(phy_name,".",sep="")))
  int<-ifelse(mr=="No","(Intercept)",paste0("trait",trait1))
  trait1<-ifelse(mr=="No","",paste0("\\",trait1))
  if(is.null(trees$tip.label)) {
    #Multiple trees
    #Check trees have node labels
    if(is.null(trees[[1]]$node.label)) {
      trees<-lapply(trees,makeNodeLabel)
    } else  {trees<-trees
    }
    
    #Create a dataframe that gives estimate for each node & tip for each tree - assumes model = tree are in same order
    results <- vector("list", length(trees))
    
    for(i in 1:length(trees)) {
      #extract model estimates
      tree<-trees[[i]]
      if(binomial == "Yes") {
        pred_states<-data.frame(tree=i,species=c(tree$node.label,tree$tip.label),estimate=boot::inv.logit(model$Sol[i,grepl(trait1,colnames(model$Sol))]))
        pred_states$pred_state<-ifelse(pred_states$estimate>cutoff,state2,ifelse(pred_states$estimate<(1-cutoff),state1,"Unknown"))
      } else  {
        pred_states<-data.frame(tree=i,species=c(tree$node.label,tree$tip.label),estimate=model$Sol[i,grepl(trait1,colnames(model$Sol))])
        pred_states$estimate<-ifelse(pred_states$species == "Node1",pred_states$estimate,pred_states$estimate + pred_states$estimate[pred_states$species == "Node1"])
      }
      results[[i]]<-pred_states
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    return(results)
    } else  {
    
    #Function for single tree
    #Check tree has node labels
      if(is.null(trees$node.label)) {
        trees<-makeNodeLabel(trees)
      } else  {trees<-trees
      }
      tree<-trees
      if(binomial == "Yes") {
        mod_est=data.frame(species=gsub(trait_sub,"",colnames(model$Sol)[grep(trait1,colnames(model$Sol))]),estimate=boot::inv.logit(posterior.mode(model$Sol[,grepl(trait1,colnames(model$Sol))])))
        mod_est$species[mod_est$species == int]<-"Node1"
        pred_states<-data.frame(species=c(tree$node.label,tree$tip.label))
        pred_states$estimate<-mod_est$estimate[match(pred_states$species,mod_est$species)]
        pred_states$pred_state<-ifelse(pred_states$estimate>cutoff,state2,ifelse(pred_states$estimate<(1-cutoff),state1,"Unknown"))
      } else  {
        pred_states=data.frame(species=gsub(trait_sub,"",colnames(model$Sol)[grep(trait1,colnames(model$Sol))]),estimate=posterior.mode(model$Sol[,grepl(trait1,colnames(model$Sol))]))
        pred_states$species[pred_states$species == int]<-"Node1"
        pred_states$estimate<-ifelse(pred_states$species == "Node1",pred_states$estimate,pred_states$estimate + pred_states$estimate[pred_states$species == "Node1"])
      }
    rownames(pred_states)<-NULL
    return(pred_states)
    }
}      
 

#******************************************************************************************
#MCMCglmm processing to transition dataset - can take list of trees (first part of function) or single trees (second part of function)
#******************************************************************************************

mcmcglmm_trans<-function(mr="No",phy_name="animal",trees,model,dat,trait_raw,trait1,trait2=NULL,state1,state2,cutoff=0.5){
  #multiple responses or not
  trait_sub<-ifelse(mr=="No",paste(phy_name,".",sep=""),paste0("trait",trait1,".",paste(phy_name,".",sep="")))
  int<-ifelse(mr=="No","(Intercept)",paste0("trait",trait1))
  trait1<-ifelse(mr=="No","",paste0("\\",trait1))
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
    #Multiple trees
    
    #Create a dataframe that gives estimate for each node & tip for each tree - assumes model = tree are in same order
    results <- vector("list", length(trees))
    
    for(i in 1:length(trees)) {
      #extract model estimates
      tree<-trees[[i]]
        pred_states<-data.frame(tree=i,species=c(tree$node.label,tree$tip.label),prob=boot::inv.logit(model$Sol[i,grepl(trait1,colnames(model$Sol))]))
        pred_states$pred_state<-ifelse(pred_states$prob>cutoff,state2,ifelse(pred_states$prob<(1-cutoff),state1,"Unknown"))
      
      #make transition dataset
      tree_df<-tidytree::as_tibble(tree)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_prob=pred_states$prob[match(ancestor,pred_states$species)],
               anc_state=pred_states$pred_state[match(ancestor,pred_states$species)],
               des_prob=pred_states$prob[match(descendant,pred_states$species)],
               des_state=pred_states$pred_state[match(descendant,pred_states$species)]) 
      
      #work out descendant nodes
      trans <- tree_df %>% dplyr::arrange(ancestor) %>%
        dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
        tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
        dplyr::mutate(des12_state=paste(des1,des2,sep=".")) 
      
      #Code different combinations
      comb1<-paste("Unknown",state1,sep=".")
      comb2<-paste(state1,"Unknown",sep=".")
      comb3<-paste("Unknown",state1,sep=".")
      comb4<-paste(state1,"Unknown",sep=".")
      comb5<-paste(state1,state1,sep=".")
      comb6<-paste(state2,state2,sep=".")
      comb7<-paste(state1,state2,sep=".")
      comb8<-paste(state2,state1,sep=".")
      
      #Recode descendant states
      trans$des12_state[trans$des12_state == "Unknown.Unknown"]<-"Unknown"
      trans$des12_state[trans$des12_state == comb1]<-state1
      trans$des12_state[trans$des12_state == comb2]<-state1
      trans$des12_state[trans$des12_state == comb3]<-state2
      trans$des12_state[trans$des12_state == comb4]<-state2
      trans$des12_state[trans$des12_state == comb5]<-state1
      trans$des12_state[trans$des12_state == comb6]<-state2
      trans$des12_state[trans$des12_state == comb7]<-"Both"
      trans$des12_state[trans$des12_state == comb8]<-"Both"
      
      #Construct transitions
      trans <- trans %>% dplyr::mutate(des_comb=paste(anc_state,des12_state,sep="."))%>%
        dplyr::mutate(des_comb=ifelse(anc_state == state2 & des12_state == "Both",comb8,des_comb),
               des_comb=ifelse(anc_state == state1 & des12_state == "Both",comb7,des_comb)) %>%
        dplyr::mutate(des_comb=ifelse(grepl("Unknown",des_comb),"Unknown",des_comb))
      
      tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
      tree_df$raw_state<-dat[,trait_raw][match(tree_df$descendant,dat[,phy_name])]
      tree_df$tree=i
      results[[i]]<-tree_df
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    results$des_comb<-ifelse(grepl("NA",results$des_comb),"NA",results$des_comb)
    return(results)
    
    } else  {
    
    #Function for single tree
    
    #Check tree has node labels
    if(is.null(trees$node.label)) {
      trees<-makeNodeLabel(trees)
    } else  {trees<-trees
    }
      tree<-trees
      mod_est=data.frame(species=gsub(trait_sub,"",colnames(model$Sol)),prob=boot::inv.logit(posterior.mode(model$Sol[,grepl(trait1,colnames(model$Sol))])))
      mod_est$species[mod_est$species == int]<-"Node1"
      pred_states<-data.frame(species=c(tree$node.label,tree$tip.label))
      pred_states$prob<-mod_est$prob[match(pred_states$species,mod_est$species)]
      pred_states$pred_state<-ifelse(pred_states$prob>cutoff,state2,ifelse(pred_states$prob<(1-cutoff),state1,"Unknown"))
      
    #add raw trait data
    pred_states$raw_state<-dat[,trait_raw][match(pred_states$species,dat[,phy_name])]
    
    #Create a dataset with ancestors and descendants
    tree_df<-tidytree::as_tibble(trees)
    data <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                       descendant=tree_df$label,
                       label_no=tree_df$node,
                       branch.length=tree_df$branch.length)
    data <- data %>% dplyr::filter(ancestor != descendant) %>%
      dplyr::mutate(anc_prob=pred_states$prob[match(ancestor,pred_states$species)],
             anc_state=pred_states$pred_state[match(ancestor,pred_states$species)],
             des_prob=pred_states$prob[match(descendant,pred_states$species)],
             des_state=pred_states$pred_state[match(descendant,pred_states$species)]) 
    data$anc_prob[data$ancestor == "Node1"]<-boot::inv.logit(mean(model$Sol[,"(Intercept)"]))
    data$anc_state[data$ancestor == "Node1"]<-ifelse(boot::inv.logit(mean(model$Sol[,"(Intercept)"]))>cutoff,state2,state1)#Root = intercept
    
    #Work out descendent states of each ancestor
    trans <- data %>% dplyr::arrange(ancestor) %>%
      dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
      tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
      dplyr::mutate(des12_state=paste(des1,des2,sep=".")) 
    
    #Code different combinations
    comb1<-paste("Unknown",state1,sep=".")
    comb2<-paste(state1,"Unknown",sep=".")
    comb3<-paste("Unknown",state1,sep=".")
    comb4<-paste(state1,"Unknown",sep=".")
    comb5<-paste(state1,state1,sep=".")
    comb6<-paste(state2,state2,sep=".")
    comb7<-paste(state1,state2,sep=".")
    comb8<-paste(state2,state1,sep=".")
    
    #Recode descendant states
    trans$des12_state[trans$des12_state == "Unknown.Unknown"]<-"Unknown"
    trans$des12_state[trans$des12_state == comb1]<-state1
    trans$des12_state[trans$des12_state == comb2]<-state1
    trans$des12_state[trans$des12_state == comb3]<-state2
    trans$des12_state[trans$des12_state == comb4]<-state2
    trans$des12_state[trans$des12_state == comb5]<-state1
    trans$des12_state[trans$des12_state == comb6]<-state2
    trans$des12_state[trans$des12_state == comb7]<-"Both"
    trans$des12_state[trans$des12_state == comb8]<-"Both"
    
    #Construct transitions
    trans <- trans %>% dplyr::mutate(des_comb=paste(anc_state,des12_state,sep="."))%>%
      dplyr::mutate(des_comb=ifelse(anc_state == state2 & des12_state == "Both",comb8,des_comb),
             des_comb=ifelse(anc_state == state1 & des12_state == "Both",comb7,des_comb)) %>%
      dplyr::mutate(des_comb=ifelse(grepl("Unknown",des_comb),"Unknown",des_comb))
    
    data$des_comb<-trans$des_comb[match(data$ancestor,trans$ancestor)]
    data$raw_state<-pred_states$raw_state[match(data$descendant,pred_states$species)]
    rownames(data)<-NULL
    data$des_comb<-ifelse(grepl("NA",data$des_comb),"NA",data$des_comb)
    return(data)
    }
}  


#******************************************************************************************
#SCM predicted ancestral states
#******************************************************************************************

scm_pred_states<-function(trees,scm_model,dat,species,trait_raw){
  #if trees$tip.label is null then must be list of trees
  
  if(is.null(trees$tip.label)) {
    
    #Multiple trees
    results <- vector("list", length(scm_model))
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    maps_per_tree = length(scm_model)/length(trees)
    
    for(i in 1:length(scm_model)) {
      pred_states<-data.frame(tree=ceiling(i/maps_per_tree),map=i,species=c(scm_model[[i]]$node.label,scm_model[[i]]$tip.label),states=c(nodes[,i],tips[,i]))
      results[[i]]<-pred_states
    }
    results = dplyr::bind_rows(results)
    
    results = results %>% dplyr::group_by(tree,species,states) %>% dplyr::summarise(n=n()) 
    results = results %>% tidyr::pivot_wider(names_from=states,values_from = n) %>% 
      dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
    
    #Convert numbers to proportion of maps
    results[,2:length(results)] = results[,2:length(results)]/length(scm_model)
    
    #Create column with most likely state
    results$state <- apply(results[, -1], 1, function(row) {
      colnames(results)[-1][which.max(row)]
    })
    
    #Create column with probability of most likely state 
    results$prob_state <- apply(results[, 2:(length(results)-1)], 1, function(row) {
      max(row)
    })
    
    
  } else  {
    results <- vector("list", length(scm_model))
    
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    for(i in 1:length(scm_model)) {
      pred_states<-data.frame(map=i,species=c(scm_model[[i]]$node.label,scm_model[[i]]$tip.label),states=c(nodes[,i],tips[,i]))
      results[[i]]<-pred_states
    }
    
    results<-dplyr::bind_rows(results)
    
    results = results %>% dplyr::group_by(species,states) %>% dplyr::summarise(n=n()) 
    results = results %>% tidyr::pivot_wider(names_from=states,values_from = n) %>% 
      dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
    
    #Convert numbers to proportion of maps
    results[,2:length(results)] = results[,2:length(results)]/length(scm_model)
    
    #Create column with most likely state
    results$state <- apply(results[, -1], 1, function(row) {
      colnames(results)[-1][which.max(row)]
    })
    
    #Create column with probability of most likely state 
    results$prob_state <- apply(results[, 2:(length(results)-1)], 1, function(row) {
      max(row)
    })
  }
  results = as.data.frame(results)
  dat = as.data.frame(dat)
  results$raw_state = dat[,trait_raw][match(results$species,dat[,species])]
  rownames(results)=NULL
  return(results)
} 

#******************************************************************************************
#SCM processing to transition dataset 
#******************************************************************************************

scm_pred_trans<-function(trees,scm_model,dat,species,trait_raw){
  dat = as.data.frame(dat)
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
    #Multiple trees
    
    #Create a dataframe that gives estimate for each node & tip for each tree - assumes model = tree are in same order
    results <- vector("list", length(scm_model))
    
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    maps_per_tree = length(scm_model)/length(trees)
    
    for(i in 1:length(scm_model)) {
      pred_states<-data.frame(species=c(scm_model[[i]]$node.label,scm_model[[i]]$tip.label),states=c(nodes[,i],tips[,i]))
      
      #make transition dataset
      tree<-scm_model[[i]]
      tree_df<-tidytree::as_tibble(tree)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_state=pred_states$states[match(ancestor,pred_states$species)],
               des_state=pred_states$states[match(descendant,pred_states$species)]) 
      
      #work out descendant nodes
      trans <- tree_df %>% dplyr::arrange(ancestor) %>%
        dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
        tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
        dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 

      tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
      tree_df$tree=ceiling(i/maps_per_tree)
      tree_df$map=i
      results[[i]]<-tree_df
    }
  } else  {
    #Single tree, multiple maps
    
    #Create a dataframe that gives estimate for each node & tip for each tree - assumes model = tree are in same order
    results <- vector("list", length(scm_model))
    
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    
    for(i in 1:length(scm_model)) {
      pred_states<-data.frame(species=c(scm_model[[i]]$node.label,scm_model[[i]]$tip.label),states=c(nodes[,i],tips[,i]))
      
      #make transition dataset
      tree<-scm_model[[i]]
      tree_df<-tidytree::as_tibble(tree)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_state=pred_states$states[match(ancestor,pred_states$species)],
                      des_state=pred_states$states[match(descendant,pred_states$species)]) 
      
      #work out descendant nodes
      trans <- tree_df %>% dplyr::arrange(ancestor) %>%
        dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
        tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
        dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
      
      tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
      tree_df$map=i
      results[[i]]<-tree_df
    }
  }
  
  results<-dplyr::bind_rows(results)
  
  results = results %>% dplyr::group_by(ancestor,descendant,label_no,branch.length,des_comb) %>% dplyr::summarise(n=n()) 
  results = results %>% tidyr::pivot_wider(names_from=des_comb,values_from = n) %>% 
    dplyr::mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  #Convert numbers to proportion of maps
  results[,5:length(results)] = results[,5:length(results)]/length(scm_model)
  
  #Create column with most likely state
  results$des_comb <- apply(results[, -c(1:5)], 1, function(row) {
    colnames(results)[-c(1:5)][which.max(row)]
  })
  
  #Create column with probability of most likely state 
  results$prob_state <- apply(results[, 5:(length(results)-1)], 1, function(row) {
    max(row)
  })
  results = as.data.frame(results)
  dat = as.data.frame(dat)
  results$raw_state<-dat[,trait_raw][match(results$descendant,dat[,species])]
  rownames(results)=NULL
  return(results)
}  

#******************************************************************************************
#HMM predicted ancestral states
#******************************************************************************************
hmm_pred_states<-function(trees,hmm_model,dat,species,trait_raw,node.states="marginal"){
  rate_mat = getStateMat4Dat(dat) # to retrieve state names
  
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
#Multiple trees
    results <- vector("list", length(hmm_model))#use length(hmm_model) here and  length(trees) to ensure they are same length (will throw and error if not)
    
    if(node.states=="marginal"){

#Multiple trees with probability of each node / species
      for(i in 1:length(trees)) {
      states = round(rbind(hmm_model[[i]]$states,hmm_model[[i]]$tip.states),4)
      pred_states = data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states)
      colnames(pred_states) = c("tree","species",rate_mat$legend)
      results[[i]]<-pred_states 
      }
      results = dplyr::bind_rows(results)
      
      #Create column with most likely state
      results$state <- apply(results[, -c(1:2)], 1, function(row) {
        colnames(results)[-c(1:2)][which.max(row)]
      })
      
      #Create column with probability of most likely state 
      results$prob_state <- apply(results[, 3:(length(results)-1)], 1, function(row) {
        max(row)
      }) 

#Multiple trees with mosty likely node / species     
      } else  {  
      
      for(i in 1:length(trees)) {  
      pred_states<-data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states=c(hmm_model[[i]]$states,hmm_model[[i]]$tip.states))
      results[[i]]<-pred_states
      }
      results = dplyr::bind_rows(results)
      #replace numerical states with trait states
      results$state = rate_mat$legend[match(results$state,names(rate_mat$legend))]
    }
    
  } else  {

#Single tree with probability of each node / species    
    #If marginal then probabilities of each state given, joint = most likely state so needs different processing
    if(node.states=="marginal"){
      #predicted ancestral states
      states = round(rbind(hmm_model$states,hmm_model$tip.states),4)
      results = data.frame(species=c(trees$node.label,trees$tip.label),states)
      colnames(results) = c("species",rate_mat$legend)
    
      #Create column with most likely state
      results$state <- apply(results[, -1], 1, function(row) {
        colnames(results)[-1][which.max(row)]
      })
      
      #Create column with probability of most likely state 
      results$prob_state <- apply(results[, 2:(length(results)-1)], 1, function(row) {
        max(row)
      }) 
      
      } else {

#Single tree with most likely node / species       
    #predicted ancestral states
    results = data.frame(species=c(trees$node.label,trees$tip.label),states=c(hmm_model$states,hmm_model$tip.states))
    #replace numerical states with trait states
    results$state = rate_mat$legend[match(results$state,names(rate_mat$legend))]
    }
  }

  results$raw_state = dat[,trait_raw][match(results$species,dat[,species])]
  rownames(results)=NULL
  return(results)
} 

#******************************************************************************************
#HMM processing to transition dataset 
#******************************************************************************************
# trees = trees_no_outgroups[1:2]
# hmm_model = multHMM_1
# dat = hmm_multi
# species = "tip"
# trait_raw= "multi_cat"
# node.states="marginal"

hmm_pred_trans<-function(trees,hmm_model,dat,species,trait_raw,node.states="marginal"){
  dat = as.data.frame(dat)
  rate_mat = getStateMat4Dat(dat)
  
#if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
#Multiple trees  
    results <- vector("list", length(hmm_model))#use length(hmm_model) here and  length(trees) to ensure they are same length (will throw and error if not)
    
#Multiple trees with probability of each state per node/species
    if(node.states=="marginal"){
    
    for(i in 1:length(trees)) {
        states = round(rbind(hmm_model[[i]]$states,hmm_model[[i]]$tip.states),4)
        pred_states = data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states)
        colnames(pred_states) = c("tree","species",rate_mat$legend)
        
        #Create column with most likely state
        pred_states$states <- apply(pred_states[, -c(1:2)], 1, function(row) {
          colnames(pred_states)[-c(1:2)][which.max(row)]
        })
        
        #Create column with probability of most likely state 
        pred_states$prob_state <- apply(pred_states[, 3:(length(pred_states)-1)], 1, function(row) {
          max(row)
        }) 
        #make transition dataset
        tree<-trees[[i]]
        tree_df<-tidytree::as_tibble(tree)
        tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                              descendant=tree_df$label,
                              label_no=tree_df$node,
                              branch.length=tree_df$branch.length)
        tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
          dplyr::mutate(anc_state=pred_states$states[match(ancestor,pred_states$species)],
                        anc_prob=pred_states$prob_state[match(ancestor,pred_states$species)],
                        des_state=pred_states$states[match(descendant,pred_states$species)],
                        des_prob=pred_states$prob_state[match(descendant,pred_states$species)]) 
        
        #work out descendant nodes
        trans <- tree_df %>% dplyr::arrange(ancestor) %>%
          dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
          tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
          dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
        
        tree_df = data.frame(tree=i,tree_df)
        results[[i]]<-tree_df 
      }
    results<-dplyr::bind_rows(results)

#Multiple trees with most likely state per node/species  
    } else {     
    
      for(i in 1:length(trees)) {
      pred_states<-data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states=c(hmm_model[[i]]$states,hmm_model[[i]]$tip.states))
      pred_states$states = rate_mat$legend[match(pred_states$states,names(rate_mat$legend))]
      
      #make transition dataset
      tree<-trees[[i]]
      tree_df<-tidytree::as_tibble(tree)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_state=pred_states$states[match(ancestor,pred_states$species)],
                      des_state=pred_states$states[match(descendant,pred_states$species)]) 
      
      #work out descendant nodes
      trans <- tree_df %>% dplyr::arrange(ancestor) %>%
        dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
        tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
        dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
      
      tree_df = data.frame(tree=i,tree_df)
      results[[i]]<-tree_df
    }
    results<-dplyr::bind_rows(results)
    }
  } else  {
    
#Single tree with probability of each state per node/species
    if(node.states=="marginal"){
      #predicted ancestral states
      states = round(rbind(hmm_model$states,hmm_model$tip.states),4)
      pred_states = data.frame(species=c(trees$node.label,trees$tip.label),states)
      colnames(pred_states) = c("species",rate_mat$legend)
      
      #Create column with most likely state
      pred_states$states <- apply(pred_states[, -1], 1, function(row) {
        colnames(pred_states)[-1][which.max(row)]
      })
      
      #Create column with probability of most likely state 
      pred_states$prob_state <- apply(pred_states[, 2:(length(pred_states)-1)], 1, function(row) {
        max(row)
      }) 
      
      #make transition dataset
      tree_df<-tidytree::as_tibble(trees)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_state=pred_states$states[match(ancestor,pred_states$species)],
                      anc_prob=pred_states$prob_state[match(ancestor,pred_states$species)],
                      des_state=pred_states$states[match(descendant,pred_states$species)],
                      des_prob=pred_states$prob_state[match(descendant,pred_states$species)]) 
      
      #work out descendant nodes
      trans <- tree_df %>% dplyr::arrange(ancestor) %>%
        dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
        tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
        dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
      
      results<-tree_df
      
    } else  { 
      
#Single tree with most likely state per node/species
    
    #predicted ancestral states
    pred_states = data.frame(species=c(trees$node.label,trees$tip.label),states=c(hmm_model$states,hmm_model$tip.states))
    pred_states$states = rate_mat$legend[match(pred_states$states,names(rate_mat$legend))]
    
    #make transition dataset
    tree_df<-tidytree::as_tibble(trees)
    tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                          descendant=tree_df$label,
                          label_no=tree_df$node,
                          branch.length=tree_df$branch.length)
    tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
      dplyr::mutate(anc_state=pred_states$states[match(ancestor,pred_states$species)],
                    des_state=pred_states$states[match(descendant,pred_states$species)]) 
    
    #work out descendant nodes
    trans <- tree_df %>% dplyr::arrange(ancestor) %>%
      dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
      tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
      dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
    
    results<-tree_df
  }
}  
  #add in information on both descendents
  
  results$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
  
  results$raw_state = dat[,trait_raw][match(results$descendant,dat[,species])]
  rownames(results)=NULL
  return(results)
}  
