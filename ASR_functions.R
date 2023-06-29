#******************************************************************************************
#MCMCglmm processing to extract values from models for ancestral nodes over multiple trees
#******************************************************************************************

pred_states_mcmcglmm<-function(mr="No",trees,phy_name="animal",model,dat,trait1,trait2=NULL,binomial="No",state1=NULL,state2=NULL,cutoff=0.5){
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
#MCMCglmm processing to transition dataset - can take single trees (first part of function) or list of trees (second part of function)
#******************************************************************************************

trans_mcmcglmm<-function(mr="No",phy_name="animal",trees,model,dat,trait_raw,trait1,trait2=NULL,state1,state2,cutoff=0.5){
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
      trans <- trans %>% dplyr::mutate(anc_des=paste(anc_state,des12_state,sep="."))%>%
        dplyr::mutate(anc_des=ifelse(anc_state == state2 & des12_state == "Both",comb8,anc_des),
               anc_des=ifelse(anc_state == state1 & des12_state == "Both",comb7,anc_des)) %>%
        dplyr::mutate(anc_des=ifelse(grepl("Unknown",anc_des),"Unknown",anc_des))
      
      tree_df$anc_des<-trans$anc_des[match(tree_df$ancestor,trans$ancestor)]
      tree_df$raw_state<-dat[,trait_raw][match(tree_df$descendant,dat[,phy_name])]
      tree_df$tree=i
      results[[i]]<-tree_df
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    results$anc_des<-ifelse(grepl("NA",results$anc_des),"NA",results$anc_des)
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
    trans <- trans %>% dplyr::mutate(anc_des=paste(anc_state,des12_state,sep="."))%>%
      dplyr::mutate(anc_des=ifelse(anc_state == state2 & des12_state == "Both",comb8,anc_des),
             anc_des=ifelse(anc_state == state1 & des12_state == "Both",comb7,anc_des)) %>%
      dplyr::mutate(anc_des=ifelse(grepl("Unknown",anc_des),"Unknown",anc_des))
    
    data$anc_des<-trans$anc_des[match(data$ancestor,trans$ancestor)]
    data$raw_state<-pred_states$raw_state[match(data$descendant,pred_states$species)]
    rownames(data)<-NULL
    data$anc_des<-ifelse(grepl("NA",data$anc_des),"NA",data$anc_des)
    return(data)
    }
}  


#******************************************************************************************
#SCM predicted ancestral states
#******************************************************************************************
pred_states_scm<-function(trees,scm_model,state1,state2,cutoff){
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    #Multiple trees
    results <- vector("list", length(scm_model))
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    for(i in 1:length(scm_model)) {
      pred_states<-data.frame(tree=i,species=c(scm_model[[i]]$node.label,scm_model[[i]]$tip.label),states=c(nodes[,i],tips[,i]))
      results[[i]]<-pred_states
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    return(results) 
  } else  {
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    
    tip_tots<-data.frame(species=scm_model[[1]]$tip.label)
    tip_tots$state1 <- rowSums(tips == state1)
    tip_tots$state2 <- rowSums(tips == state2)
    
    nodes_tots<-data.frame(species=scm_model[[1]]$node.label)
    nodes_tots$state1 <- rowSums(nodes == state1)
    nodes_tots$state2 <- rowSums(nodes == state2)
    
    pred_states<-rbind(nodes_tots,tip_tots)
    pred_states$prob<-pred_states$state1/(pred_states$state1+pred_states$state2)
    pred_states$state<-ifelse(pred_states$prob>cutoff,state2,ifelse(pred_states$prob<(1-cutoff),state1,"Unknown"))
    rownames(pred_states)<-NULL
    return(pred_states)
  }
} 

#******************************************************************************************
#SCM processing to transition dataset 
#******************************************************************************************

trans_scm<-function(trees,scm_model,dat,trait_raw,state1,state2,cutoff = 0.5){
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
    #Multiple trees
    
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
      trans <- trans %>% dplyr::mutate(anc_des=paste(anc_state,des12_state,sep="."))%>%
        dplyr::mutate(anc_des=ifelse(anc_state == state2 & des12_state == "Both",comb8,anc_des),
               anc_des=ifelse(anc_state == state1 & des12_state == "Both",comb7,anc_des)) %>%
        dplyr::mutate(anc_des=ifelse(grepl("Unknown",anc_des),"Unknown",anc_des))
      
      tree_df$anc_des<-trans$anc_des[match(tree_df$ancestor,trans$ancestor)]
      tree_df$raw_state<-dat[,trait_raw][match(tree_df$descendant,dat[,phy_name])]
      tree_df$tree=i
      results[[i]]<-tree_df
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    return(results)
    
  } else  {
    
    #predicted ancestral states
    nodes<-data.frame(getStates(scm_model,type=c("nodes")))
    tips<-data.frame(getStates(scm_model,type=c("tips")))
    
    tip_tots<-data.frame(species=scm_model[[1]]$tip.label)
    tip_tots$state1 <- rowSums(tips == state1)
    tip_tots$state2 <- rowSums(tips == state2)
    
    nodes_tots<-data.frame(species=scm_model[[1]]$node.label)
    nodes_tots$state1 <- rowSums(nodes == state1)
    nodes_tots$state2 <- rowSums(nodes == state2)
    
    pred_states<-rbind(nodes_tots,tip_tots)
    pred_states$prob<-pred_states$state1/(pred_states$state1+pred_states$state2)
    pred_states$state<-ifelse(pred_states$prob>cutoff,state2,ifelse(pred_states$prob<(1-cutoff),state1,"Unknown"))
    
    #Create a dataset with ancestors and descendants
    tree_df<-tidytree::as_tibble(trees)
    data <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                       descendant=tree_df$label,
                       label_no=tree_df$node,
                       branch.length=tree_df$branch.length)
    data <- data %>% dplyr::filter(ancestor != descendant) %>%
      dplyr::mutate(anc_prob=pred_states$prob[match(ancestor,pred_states$species)],
             anc_state=pred_states$state[match(ancestor,pred_states$species)],
             des_prob=pred_states$prob[match(descendant,pred_states$species)],
             des_state=pred_states$state[match(descendant,pred_states$species)]) 
    
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
    trans <- trans %>% dplyr::mutate(anc_des=paste(anc_state,des12_state,sep="."))%>%
      dplyr::mutate(anc_des=ifelse(anc_state == state2 & des12_state == "Both",comb8,anc_des),
             anc_des=ifelse(anc_state == state1 & des12_state == "Both",comb7,anc_des)) %>%
      dplyr::mutate(anc_des=ifelse(grepl("Unknown",anc_des),"Unknown",anc_des))
    
    data$anc_des<-trans$anc_des[match(data$ancestor,trans$ancestor)]
    data$raw_state<-pred_states$raw_state[match(data$descendant,pred_states$species)]
    rownames(data)<-NULL
    return(data)
  }
}  

#******************************************************************************************
#corHMM models run on list of trees
#******************************************************************************************
corHMM_trees<-function(trees,data,rate_cats=1,node.states="joint"){
    results <- vector("list", length(trees))

    for(i in 1:length(trees)) {
      mod<-corHMM(phy = trees[[i]], data = data, rate.cat = rate_cats,node.states=node.states)
      results[[i]]<-mod
    }

    return(results) 
  } 

#******************************************************************************************
#corHMM predicted ancestral states
#******************************************************************************************
#for 2 rate models for a binary trait
pred_states_corHMM<-function(trees,corHMM_models,state1,state2){
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    #Multiple trees
    results <- vector("list", length(trees))

    for(i in 1:length(trees)) {
      pred_states<-data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),state_2=c(corHMM_models[[i]]$states,corHMM_models[[i]]$tip.states))
      pred_states$state_2[pred_states$state_2 == 1]<-paste0(state1,"_1")
      pred_states$state_2[pred_states$state_2 == 2]<-paste0(state2,"_1")
      pred_states$state_2[pred_states$state_2 == 3]<-paste0(state1,"_2")
      pred_states$state_2[pred_states$state_2 == 4]<-paste0(state2,"_2")
      pred_states$state_2_sum=ifelse(grepl(state1,pred_states$state_2),state1,state2)
      
      results[[i]]<-pred_states
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    return(results) 
  } else  {
    #predicted ancestral states
    pred_states<-data.frame(species=c(trees$node.label,trees$tip.label),state_2=c(corHMM_models$states,corHMM_models$tip.states))
    pred_states$state_2[pred_states$state_2 == 1]<-paste0(state1,"_1")
    pred_states$state_2[pred_states$state_2 == 2]<-paste0(state2,"_1")
    pred_states$state_2[pred_states$state_2 == 3]<-paste0(state1,"_2")
    pred_states$state_2[pred_states$state_2 == 4]<-paste0(state2,"_2")
    pred_states$state_2_sum=ifelse(grepl(state1,pred_states$state_2),state1,state2)
    
    rownames(pred_states)<-NULL
    return(pred_states)
  }
} 


#******************************************************************************************
#corHMM processing to transition dataset 
#******************************************************************************************

trans_corHMM<-function(trees,corHMM_models,dat,trait_raw,state1,state2){
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
    #Multiple trees
    
    #Create a dataframe that gives estimate for each node & tip for each tree - assumes model = tree are in same order
    results <- vector("list", length(trees))

    for(i in 1:length(trees)) {
      pred_states<-data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),state_2=c(corHMM_models[[i]]$states,corHMM_models[[i]]$tip.states))
      pred_states$state_2[pred_states$state_2 == 1]<-paste0(state1,"_1")
      pred_states$state_2[pred_states$state_2 == 2]<-paste0(state2,"_1")
      pred_states$state_2[pred_states$state_2 == 3]<-paste0(state1,"_2")
      pred_states$state_2[pred_states$state_2 == 4]<-paste0(state2,"_2")
      pred_states$state_2_sum=ifelse(grepl(state1,pred_states$state_2),state1,state2)

      #make transition dataset
      tree<-trees[[i]]
      tree_df<-tidytree::as_tibble(tree)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_state=pred_states$state_2_sum[match(ancestor,pred_states$species)],
               des_state=pred_states$state_2_sum[match(descendant,pred_states$species)]) 
      
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
      trans <- trans %>% dplyr::mutate(anc_des=paste(anc_state,des12_state,sep="."))%>%
        dplyr::mutate(anc_des=ifelse(anc_state == state2 & des12_state == "Both",comb8,anc_des),
               anc_des=ifelse(anc_state == state1 & des12_state == "Both",comb7,anc_des)) %>%
        dplyr::mutate(anc_des=ifelse(grepl("Unknown",anc_des),"Unknown",anc_des))
      
      tree_df$anc_des<-trans$anc_des[match(tree_df$ancestor,trans$ancestor)]
      tree_df$raw_state<-dat[,trait_raw][match(tree_df$descendant,dat[,phy_name])]
      tree_df$tree=i
      results[[i]]<-tree_df
    }
    results<-dplyr::bind_rows(results)
    rownames(results)<-NULL
    return(results)
    
  } else  {
    
    #predicted ancestral states
    pred_states<-data.frame(species=c(trees$node.label,trees$tip.label),state_2=c(corHMM_models$states,corHMM_models$tip.states))
    pred_states$state_2[pred_states$state_2 == 1]<-paste0(state1,"_1")
    pred_states$state_2[pred_states$state_2 == 2]<-paste0(state2,"_1")
    pred_states$state_2[pred_states$state_2 == 3]<-paste0(state1,"_2")
    pred_states$state_2[pred_states$state_2 == 4]<-paste0(state2,"_2")
    pred_states$state_2_sum=ifelse(grepl(state1,pred_states$state_2),state1,state2)
    
    #Create a dataset with ancestors and descendants
    tree_df<-tidytree::as_tibble(trees)
    data <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                       descendant=tree_df$label,
                       label_no=tree_df$node,
                       branch.length=tree_df$branch.length)
    data <- data %>% dplyr::filter(ancestor != descendant) %>%
      dplyr::mutate(anc_state=pred_states$state_2_sum[match(ancestor,pred_states$species)],
             des_state=pred_states$state_2_sum[match(descendant,pred_states$species)]) 
    
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
    trans <- trans %>% dplyr::mutate(anc_des=paste(anc_state,des12_state,sep="."))%>%
      dplyr::mutate(anc_des=ifelse(anc_state == state2 & des12_state == "Both",comb8,anc_des),
             anc_des=ifelse(anc_state == state1 & des12_state == "Both",comb7,anc_des)) %>%
      dplyr::mutate(anc_des=ifelse(grepl("Unknown",anc_des),"Unknown",anc_des))
    
    data$anc_des<-trans$anc_des[match(data$ancestor,trans$ancestor)]
    data$raw_state<-pred_states$raw_state[match(data$descendant,pred_states$species)]
    rownames(data)<-NULL
    return(data)
  }
}  
