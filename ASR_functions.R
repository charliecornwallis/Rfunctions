#******************************************************************************************
#MCMCglmm processing to extract values from models for ancestral nodes - can take list of trees (first part of function) or single trees (second part of function)
#******************************************************************************************

mcmcglmm_pred_states<-function(mr="No",trees,phy_name="animal",model,dat,trait1,trait2=NULL,link="identity",state1=NULL,state2=NULL,cutoff=0.5){
  #remove formatting
  dat = as.data.frame(dat)
  
  #link function
  link_fun <- function(link) {
    if (link == "identity") {
      return(function(x) { 1 * x })
    } else if (link == "logit") {
      return(function(x) { boot::inv.logit(x) })
    } else if (link == "probit") {
      return(function(x) { VGAM::probitlink(x, inverse = T) })
    } else {
      stop("Invalid link function specified. Choose 'identity', 'logit', or 'probit'.")
    }
  }
  inv.link = link_fun(link)
  
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
      if(link != "identity") {
        pred_states<-data.frame(tree=i,species=c(tree$node.label,tree$tip.label),estimate=inv.link(model$Sol[i,grepl(trait1,colnames(model$Sol))]))
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
    if(link != "identity") {
      mod_est=data.frame(species=gsub(trait_sub,"",colnames(model$Sol)[grep(trait1,colnames(model$Sol))]),estimate=inv.link(posterior.mode(model$Sol[,grepl(trait1,colnames(model$Sol))])))
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

mcmcglmm_trans<-function(mr="No",phy_name="animal",trees,link="logit",model,dat,trait_raw,trait1,trait2=NULL,state1,state2,cutoff=0.5){
  #remove formatting
  dat = as.data.frame(dat)
  
  #link function
  link_fun <- function(link) {
    if (link == "identity") {
      return(function(x) { 1 * x })
    } else if (link == "logit") {
      return(function(x) { boot::inv.logit(x) })
    } else if (link == "probit") {
      return(function(x) { VGAM::probitlink(x, inverse = T) })
    } else {
      stop("Invalid link function specified. Choose 'identity', 'logit', or 'probit'.")
    }
  }
  inv.link = link_fun(link)
  
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
      pred_states<-data.frame(tree=i,species=c(tree$node.label,tree$tip.label),prob=boot::inv.link(model$Sol[i,grepl(trait1,colnames(model$Sol))]))
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
    mod_est=data.frame(species=gsub(trait_sub,"",colnames(model$Sol)),prob=boot::inv.link(posterior.mode(model$Sol[,grepl(trait1,colnames(model$Sol))])))
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
    data$anc_prob[data$ancestor == "Node1"]<-boot::inv.link(mean(model$Sol[,"(Intercept)"]))
    data$anc_state[data$ancestor == "Node1"]<-ifelse(boot::inv.link(mean(model$Sol[,"(Intercept)"]))>cutoff,state2,state1)#Root = intercept
    
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
    
    #Multiple trees with multiple maps per tree
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
    
    #Single tree with multiple maps
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
  
  dat = as.data.frame(dat) # remove any formatting
  
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
    #Multiple trees with multiple maps per tree
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
    
    #Single tree with multiple maps
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
  
  #***********************************************
  #create data frame with ancestor and descendents
  if(is.null(trees$tip.label)) {
    
    results2 <- vector("list", length(scm_model))
    
    #make transition dataset for each tree
    for(i in 1:length(trees)) {
      tree<-trees[[i]]
      tree_df<-tidytree::as_tibble(tree)
      tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                            descendant=tree_df$label,
                            label_no=tree_df$node,
                            branch.length=tree_df$branch.length)
      tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
        dplyr::mutate(anc_state=results$state[match(ancestor,results$species)],
                      anc_prob=results$prob_state[match(ancestor,results$species)],
                      des_state=results$state[match(descendant,results$species)],
                      des_prob=results$prob_state[match(descendant,results$species)]) 
      
      #work out descendant nodes
      trans <- tree_df %>% dplyr::arrange(ancestor) %>%
        dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
        tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
        dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
      
      #add in information on both descendents
      tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
      tree_df = data.frame(tree=i,tree_df)
      results2[[i]]<-tree_df 
    }  
    results2<-dplyr::bind_rows(results2)
    
  } else  {
    tree_df<-tidytree::as_tibble(trees)
    tree_df <- data.frame(ancestor=tree_df$label[match(tree_df$parent,tree_df$node)],
                          descendant=tree_df$label,
                          label_no=tree_df$node,
                          branch.length=tree_df$branch.length)
    tree_df <- tree_df %>% dplyr::filter(ancestor != descendant) %>%
      dplyr::mutate(anc_state=results$state[match(ancestor,results$species)],
                    anc_prob=results$prob_state[match(ancestor,results$species)],
                    des_state=results$state[match(descendant,results$species)],
                    des_prob=results$prob_state[match(descendant,results$species)]) 
    
    #work out descendant nodes
    trans <- tree_df %>% dplyr::arrange(ancestor) %>%
      dplyr::mutate(des=rep(c("des1","des2"),length(unique(ancestor)))) %>% 
      tidyr::pivot_wider(id_cols=c(ancestor,anc_state),names_from=c(des),values_from = c(des_state)) %>%
      dplyr::mutate(des_comb=paste(des1,des2,sep=".")) 
    
    #add in information on both descendents
    tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
    results2<-tree_df
  }
  
  results2$raw_state = dat[,trait_raw][match(results2$descendant,dat[,species])]
  rownames(results2)=NULL
  return(results2)
}

#******************************************************************************************
#HMM predicted ancestral states
#******************************************************************************************
hmm_pred_states<-function(trees,hmm_model,dat,species,trait_raw,node.states="marginal"){
  
  dat = as.data.frame(dat) # remove any formatting
  
  #Retrieve state names
  if(hmm_model$rate.cat >1) {
    
    state_names = character()
    for(i in 1:hmm_model$rate.cat) {
      tmp = paste(getStateMat4Dat(dat)$legend,"_rate",i,sep="")
      state_names = c(state_names,tmp)
    }
    names(state_names) = as.character(seq(hmm_model$rate.cat*length(getStateMat4Dat(dat)$legend))) 
  } else  {  
    state_names = getStateMat4Dat(dat)$legend 
  }
  
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    
    #Multiple trees
    results <- vector("list", length(hmm_model))#use length(hmm_model) here and  length(trees) to ensure they are same length (will throw and error if not)
    
    if(node.states=="marginal"){
      
      #Multiple trees with probability of each node / species
      for(i in 1:length(trees)) {
        states = round(rbind(hmm_model[[i]]$states,hmm_model[[i]]$tip.states),4)
        pred_states = data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states)
        colnames(pred_states) = c("tree","species",state_names)
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
      results$state = state_names[match(results$state,names(state_names))]
    }
    
  } else  {
    
    #Single tree with probability of each node / species    
    #If marginal then probabilities of each state given, joint = most likely state so needs different processing
    if(node.states=="marginal"){
      #predicted ancestral states
      states = round(rbind(hmm_model$states,hmm_model$tip.states),4)
      results = data.frame(species=c(trees$node.label,trees$tip.label),states)
      colnames(results) = c("species",state_names)
      
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
      results$state = state_names[match(results$state,names(state_names))]
    }
  }
  
  results$raw_state = dat[,trait_raw][match(results$species,dat[,species])]
  rownames(results)=NULL
  return(results)
} 

#******************************************************************************************
#HMM processing to transition dataset 
#******************************************************************************************
hmm_pred_trans<-function(trees,hmm_model,dat,species,trait_raw,node.states="marginal"){
  
  dat = as.data.frame(dat) # remove any formatting
  
  #Retrieve state names
  if(hmm_model$rate.cat >1) {
    
    state_names = character()
    for(i in 1:hmm_model$rate.cat) {
      tmp = paste(getStateMat4Dat(dat)$legend,"_rate",i,sep="")
      state_names = c(state_names,tmp)
    }
    names(state_names) = as.character(seq(hmm_model$rate.cat*length(getStateMat4Dat(dat)$legend))) 
  } else  {  
    state_names = getStateMat4Dat(dat)$legend 
  }
  
  #if trees$tip.label is null then must be list of trees
  if(is.null(trees$tip.label)) {
    #Multiple trees  
    results <- vector("list", length(hmm_model))#use length(hmm_model) here and  length(trees) to ensure they are same length (will throw and error if not)
    
    #Multiple trees with probability of each state per node/species
    if(node.states=="marginal"){
      
      for(i in 1:length(trees)) {
        states = round(rbind(hmm_model[[i]]$states,hmm_model[[i]]$tip.states),4)
        pred_states = data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states)
        colnames(pred_states) = c("tree","species",state_names)
        
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
        
        #add in information on both descendents
        tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
        tree_df = data.frame(tree=i,tree_df)
        results[[i]]<-tree_df 
      }
      results<-dplyr::bind_rows(results)
      
      #Multiple trees with most likely state per node/species  
    } else {     
      
      for(i in 1:length(trees)) {
        pred_states<-data.frame(tree=i,species=c(trees[[i]]$node.label,trees[[i]]$tip.label),states=c(hmm_model[[i]]$states,hmm_model[[i]]$tip.states))
        pred_states$states = state_names[match(pred_states$states,names(state_names))]
        
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
        
        #add in information on both descendents
        tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
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
      colnames(pred_states) = c("species",state_names)
      
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
      
      #add in information on both descendents
      tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
      results<-tree_df
      
    } else  { 
      
      #Single tree with most likely state per node/species
      
      #predicted ancestral states
      pred_states = data.frame(species=c(trees$node.label,trees$tip.label),states=c(hmm_model$states,hmm_model$tip.states))
      pred_states$states = state_names[match(pred_states$states,names(state_names))]
      
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
      
      #add in information on both descendents
      tree_df$des_comb<-trans$des_comb[match(tree_df$ancestor,trans$ancestor)]
      results<-tree_df
    }
  }  
  #add in information on both descendents
  
  results$raw_state = dat[,trait_raw][match(results$descendant,dat[,species])]
  rownames(results)=NULL
  return(results)
}  
