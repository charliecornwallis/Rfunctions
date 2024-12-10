#***************************************
#Function for processing MCMCglmm models#
#***************************************

#Trouble shooting tools----
# model=m
# responses=NULL
# dist_var=c(0)
# fixed_names=NULL
# Include_random = "yes"
# variances=
# covariances =
# randomvar_names=
# randomcovar_names =
# cor_diffs = 
# padding=3
# ginv="animal"
# fixed_diffinc="all"
# fixed_diff_diffs =NULL
# Include_random = "yes"
# padding=4
# dec_PM=2
# pvalues="exclude"
# S2var=0
# cor_diff=NULL

#Function ----

MCMCglmmProc<-function(model=NULL,responses=NULL,dist_var=NULL,ginv="animal",S2var=0,start_row=NULL,workbook=NULL, create_sheet="yes",sheet="sheet1",title="",fixed_names=NULL,fixed_diffinc="none",fixed_diff_diffs =NULL,variances=NULL,covariances=NULL,randomvar_names=NULL,randomcovar_names=NULL,Include_random = "no",padding=4,dec_PM=2,pvalues="include",cor_diffs=NULL)
{ 
  #Explanation of terms ----
  #model = MCMCglmm model
  #response = list of responses (e.g. c(trait1,trait2))
  #dist_var = distribution variance associated with link function: e.g."gaussian" = 0, "log" = log(1 + log(exp(intercept + 0.5*sumRE)), "logit" = pi^2/3, "probit" = 1. This is a complicated issue and should be given careful thought: see Nakagawa et al 2017
  #S2var = sampling variance if known - useful for meta-analyses
  #For Multi-response models can provide a list of dist_var and S2var corresponding to each response trait
  #start_row=starting row of workbook to add data to if NULL put data in first empty row 
  #workbook = adds data if specified, otherwise will make a new e.g. Results
  #create_sheet = should a new sheet be created e.g."yes" vs "no"
  #sheet= name of sheet "Analysis 1"
  #title = Title of table in e.g. "Table 1"
  #fixed_names = what you want fixed effects to be called e.g c("Intercept","Season length")
  #fixed_diffinc = differences between fixed effects to be included in output. Terms that come first in fixed_names have to come first in the comparison.
  #fixed_diff_diffs = calculates differences between differences e.g. c("effect1 vs effect2 - effect3 vs effect4"). Must exactly match names of fixed effects and be separate by " - "
  #variances = names of variance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed. If left as NULL then taken from the model object.
  #covariances = names of covariance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed. If left as NULL then they are not outputted.
  #randomvar_names = what you want random effect variances to be called. Note reserved term animal is renamed to phylogeny
  #randomcovar_names = what you want random effect covariances to be called. Note for correlations to be calculated this has to be same as variance1 : variance2
  #Include_random = should random effects be included in output
  #Padding = space between tables when outputting multiple models to same sheet
  #dec_PM = number of decimals given for posterior mode and CIs of fixed and random effects
  #responses = specify response variables can take multiple values for multi response
  #pvalues = exclusion of pMCMC values for fixed effects - "exclude" = exclude all, "include" = include all or index pvalues to be excluded e.g. "c(1,3)" removes 1st and 3rd, c(1:7) removes 1 to 7. Note pMCMC will still be calculated for fixed effect comparisons.
  #cor_diff = calculates differences between correlations. Should be specified in the same way as fixed_diffs e.g. c("cor1 vs cor2","cor1 vs cor3"...) 
  
  #Naming aid ----
  #variances=c("trait1:trait1.animal","trait2:trait2.animal","trait3:trait3.animal","trait4:trait4.animal","trait5:trait5.animal","trait6:trait6.animal","trait7:trait7.animal","trait8:trait8.animal","trait9:trait9.animal","trait1:trait1.units","trait2:trait2.units","trait3:trait3.units","trait4:trait4.units","trait5:trait5.units","trait6:trait6.units","trait7:trait7.units","trait8:trait8.units","trait9:trait9.units"),
  #covariances =c("trait2:trait1.animal","trait3:trait1.animal","trait4:trait1.animal","trait5:trait1.animal","trait6:trait1.animal","trait7:trait1.animal","trait8:trait1.animal","trait9:trait1.animal","trait2:trait1.units","trait3:trait1.units","trait4:trait1.units","trait5:trait1.units","trait6:trait1.units","trait7:trait1.units","trait8:trait1.units","trait9:trait1.units"),
  #randomvar_names=c("Phylogeny Trait1","Phylogeny Trait2","Phylogeny Trait3","Phylogeny Trait4","Phylogeny Trait5","Phylogeny Trait6","Phylogeny Trait7", "Phylogeny Trait8","Phylogeny Trait9","Residual Trait1","Residual Trait2","Residual Trait3","Residual Trait4","Residual Trait5","Residual Trait6","Residual Trait7","Residual Trait8","Residual Trait9"),
  #randomcovar_names =c("Phylogeny Trait2 : Phylogeny Trait1","Phylogeny Trait3 : Phylogeny Trait1","Phylogeny Trait4 : Phylogeny Trait1","Phylogeny Trait5 : Phylogeny Trait1","Phylogeny Trait6 : Phylogeny Trait1","Phylogeny Trait7 : Phylogeny Trait1","Phylogeny Trait8 : Phylogeny Trait1","Phylogeny Trait9 : Phylogeny Trait1","Residual Trait2 : Residual Trait1","Residual Trait3 : Residual Trait1","Residual Trait4 : Residual Trait1","Residual Trait5 : Residual Trait1","Residual Trait6 : Residual Trait1","Residual Trait7 : Residual Trait1","Residual Trait8 : Residual Trait1","Residual Trait9 : Residual Trait1"),
  
  #Load packages and naming ----
  pacman::p_load(MCMCglmm,coda,openxlsx,stringdist,kableExtra)

  #****************************************************
  #Check terms are specified correctly and relabel terms ----
  #****************************************************
  if (is.null(fixed_names) & any(fixed_diffinc != "none")) {
    stop("fixed names needs to be specified for fixed_diffinc to be calculated")
  } else {
  }

  #Make sure distribution variances are specified
  if (is.null(dist_var)) {
    stop("please specify distribution variances e.g. distr=c(0) for a model with a single gaussian response, distr=c(0,0) for a model with two gaussian responses ect..")
  } else {
  }
  
  #If response(s) not specified
  if (is.null(responses)) {
    
    if(length(dist_var) ==1){
      responses = model$Fixed$formula[2]
      responses = sub("()","",responses)
      }
    else {
      stop("Need to specify responses for multi-response models")
    }
    
  } else {
    responses = responses
  }
  
  #Remove unwanted effects
  #Remove mev from random effect
  model$VCV<-model$VCV[, colnames(model$VCV) !="sqrt(mev):sqrt(mev).meta"]
  
  #Rename model fixed effects
  if(is.null(fixed_names)) {
    fixed_names<- colnames(model$Sol)
  } else  {colnames(model$Sol)<-fixed_names
  }
  
  #****************************************************
  #Fixed effects and Differences between levels ----
  #****************************************************
  #Main effects
  nF=dim(model$Sol)[2]
  fe1=paste(round(posterior.mode(model$Sol),dec_PM)," (",round(HPDinterval(model$Sol)[,1],dec_PM), ", ",round(HPDinterval(model$Sol)[,2],dec_PM),")",sep="")
  
  #P values using summary.MCMCglmm code
  #Pvalues = option to exclude 
  if(any(pvalues == "include")) {
    fe1_p=pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:nF, drop = FALSE] > 0)/dim(model$Sol)[1], 1 - colSums(model$Sol[, 1:nF, drop = FALSE] > 0)/dim(model$Sol)[1]))*2
    fe1_p=round(as.numeric(fe1_p),3)
  } else  {
    if(any(pvalues == "exclude")) {
      fe1_p=rep("-",length(fixed_names))
    } else  {
      fe1_p=pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:nF, drop = FALSE] > 0)/dim(model$Sol)[1], 1 - colSums(model$Sol[, 1:nF, drop = FALSE] > 0)/dim(model$Sol)[1]))*2
      fe1_p=round(as.numeric(fe1_p),3)
      fe1_p[pvalues]<-"-"
    }
  }
  
  fe1=data.frame(Fixed_Effects=colnames(model$Sol),Estimates=fe1, pMCMC=fe1_p)
  
  #function for calculating differences between all columns of a matrix
  pairwise.diffs <- function(x, nF=1)
  {if(is.matrix(x) & nF>1) {
    #Differences
    # create column combination pairs
    col.diffs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
    #col.diffs <- col.diffs[col.diffs[, 1] < col.diffs[, 2], , drop = FALSE]
    
    #pairwise differences 
    result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
    
    # set colnames
    if(is.null(colnames(x)))
      colnames(x) <- 1:ncol(x)
    colnames(result) <- paste(colnames(x)[col.diffs[, 1]], " vs ", 
                              colnames(x)[col.diffs[, 2]], sep = "")
    ndiffs<-dim(result)[2]
    result<-as.mcmc(result)
    fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
    fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
    fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),3), check.names=FALSE)
    return(fe2)
  } 
  }
  
  ##fixed_diffinc ----
  if(any(fixed_diffinc == "none")) {
    fixed=fe1
  } else  {
    #Differences between fixed effects
    fe2<-pairwise.diffs(model$Sol,nF=nF)
    fixed=rbind(fe1,fe2)
    fixed = fixed %>% dplyr::filter(Fixed_Effects %in% c(fixed_names,fixed_diffinc) == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  
  ##fixed_diff_diffs ----
  if(is.null(fixed_diff_diffs)) {
    fixed=fixed
  } else  {
    #create matrix of diffs
    #function for calculating differences between all columns of a matrix
    pairwise.diffs.mat <- function(x)
    {if(is.matrix(x)) {
      # create column combination pairs
      prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
      col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
      
      #pairwise differences 
      result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
      # set colnames
      if(is.null(colnames(x)))
        colnames(x) <- 1:ncol(x)
      colnames(result) <- paste(colnames(x)[col.diffs[, 1]], " vs ", 
                                colnames(x)[col.diffs[, 2]], sep = "")
      ndiffs<-dim(result)[2]
      result<-as.mcmc(result)
      return(result)
    } 
    }
    
    #Create a matrix of differences
    diffs.mat<-pairwise.diffs.mat(model$Sol)
    
    #Calculate differences between differences
    pairwise.diffs2 <- function(x)
    {if(is.matrix(x)) {
      #Differences
      # create column combination pairs
      prs <- cbind(rep(1:ncol(x), each = ncol(x)), 1:ncol(x))
      col.diffs <- prs[prs[, 1] < prs[, 2], , drop = FALSE]
      
      #pairwise differences 
      result <- x[, col.diffs[, 1]] - x[, col.diffs[, 2], drop = FALSE]
      # set colnames
      if(is.null(colnames(x)))
        colnames(x) <- 1:ncol(x)
      colnames(result) <- paste(colnames(x)[col.diffs[, 1]], " - ", 
                                colnames(x)[col.diffs[, 2]], sep = "")
      ndiffs<-dim(result)[2]
      result<-as.mcmc(result)
      fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
      fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
      fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),3), check.names=FALSE)
      return(fe2)
    } 
    }
    diffs.diffs<-pairwise.diffs2(diffs.mat)
    #Select differences that are specified
    diffs.diffs = diffs.diffs %>% dplyr::filter(Fixed_Effects %in% fixed_diff_diffs == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    #add results to fixed effects
    fixed = rbind(fixed,diffs.diffs)
  }
  
  #Fixed effects
  fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
  fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
  #Fixed differences
  fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
  fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
  
  #****************************************************
  #Excel output: fixed effects ----
  #****************************************************
  #Create excel workbook if not specified
  if(is.null(workbook)) {
    workbook<- createWorkbook()
  } else  {workbook
  }
  #Create excel sheet if not specified
  if(create_sheet == "yes") {
    addWorksheet(workbook, sheet)
  } else  {sheet=sheet
  }
  #Calculate start row
  if(is.null(start_row)) {
    start_row=1
  } else  {start_row=start_row
  }
  
  #Table headers
  hs1 <- createStyle(fgFill = "white", halign = "LEFT", textDecoration = "bold",border = "TopBottom")
  hs2 <- createStyle(halign = "LEFT",border = "TopBottom",textDecoration = "bold")
  #table title
  header=data.frame(col1=c(""),col2=c(""),col3=c(""))
  colnames(header)<-c(title,"","")
  writeData(workbook, sheet, header, startCol = 1, startRow = start_row,headerStyle = hs1)
  
  #Fixed effects
  writeData(workbook, sheet, fixedeff, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
  row_nums = start_row+dim(header)[1] + dim(fixedeff)[1]+1
  
  if(any(fixed_diffinc == "none")) { #Do not write fixeddiff dataframe if fixed_diffinc == "none"
  } else  {
    writeData(workbook, sheet, fixeddiff, startCol = 1, startRow = start_row+dim(header)[1] +dim(fixedeff)[1]+1,headerStyle = hs2)
    row_nums = row_nums + dim(fixeddiff)[1]+1
  }
  
  #Bold pMCMC values less than 0.05
  bolding<-createStyle(textDecoration="bold")
  conditionalFormatting(workbook, sheet, cols=3, rows=1:10000, rule="<0.05", style = bolding)

  #****************************************************
  #Random effects ----
  #****************************************************
  #Random effects: variances
  #Should they be estimated or not?
  if(Include_random == "no") {
  #No random effects or none to be estimated
  } else  {
  
  #If only residual variance  
  if(dim(as.matrix(model$VCV))[2]==1) {
    
    rand1=paste(round(posterior.mode(model$VCV),dec_PM)," (",round(HPDinterval(model$VCV)[,1],dec_PM), ", ",round(HPDinterval(model$VCV)[,2],dec_PM),")",sep="")
    randomVar<-data.frame("Random Effects"="Residual","Posterior Mode (CI)"=c(rand1),"I2 % (CI)"=100, check.names=FALSE)
    
    } else  {  
  
  #Rename phylogeny or pedigree term (ginv=) to animal
  colnames(model$VCV) <- sub(ginv, "animal", colnames(model$VCV))
  variances <- sub(ginv, "animal", variances)
  
  #if variances are not specified make them the same as those in the model
  if(is.null(variances)) {
    
    #Multi-response models
    if(length(responses)>1){
      var_ids<-seq(from=1,to=length(responses)^2,by=length(responses)+1)#indices of variances
      for(i in 1:length(model$Random$nfl)+1){ #The indices corresponding to each response variable
        var_ids<-c(var_ids,var_ids+max(var_ids))
      }
      variances<-variances[var_ids]#Pick out variances
    } else  {
      #if single response model then all random effects assumed to be variances. Might not always be the case (e.g. for random regression ) in which case variances and covariance need to be specified
      variances<-colnames(model$VCV)
    } 
  } else  {
  }
  
  ##Separate out variance and covariance terms ----
  var_terms<-model$VCV[,variances]
  
  #rename random effects if specified
  if(is.null(randomvar_names)) {
    colnames(var_terms)<- variances
    randomvar_names<-variances
  } else  {colnames(var_terms)<-randomvar_names
  }
  
  #****************************************************
  ##Variance output ----
  #Summary of variance components
  rand1=paste(round(posterior.mode(var_terms),dec_PM)," (",round(HPDinterval(var_terms)[,1],dec_PM), ", ",round(HPDinterval(var_terms)[,2],dec_PM),")",sep="")
  
  #Calculate % variation explained by each random effect
  
  #Setup sampling variances for multi-response models
  if(length(responses)>1 & sum(S2var) ==0){
    S2var<-rep(0,length(responses))
  } else  {S2var<-S2var
  }
  
  #Separate variances from each trait for multi-response models and create a var_sum with trait specific values
  icc_all<-mcmc()
  
  if(length(responses)>1){
    
    for(i in 1:length(responses)){
      tvar<-var_terms
      colnames(tvar)<-variances
      tvar<-tvar[,grepl(paste0(responses[i],"."),colnames(tvar),fixed=T)]
      tsumvar<-rowSums(tvar) #Calculate sum of variances
      
      tsumvar<-tsumvar + dist_var[i] #Add distribution variance
      tot_tsumvar<-tsumvar+S2var[i] #Add sampling variance
      
      if(any(grepl("animal",colnames(tvar)))) {
        t_icc<-(tvar/tot_tsumvar)*100
        #If a phylogeny is included then divide all terms apart from phylogeny by total variance. ICC of phylogeny should exclude sampling variance as phylogenetic relatedness is a "fixed random effect": see Nakagawa & Santos 2012
      } else  {
        t_icc<-tvar
        t_icc[,grepl("animal",colnames(t_icc))]<-(t_icc[,grepl("animal",colnames(t_icc))]/(tsumvar))*100
        t_icc[,!grepl("animal",colnames(t_icc))]<-(t_icc[,!grepl("animal",colnames(t_icc))]/(tot_tsumvar))*100
      }
      
      #Need to rename colnames to match randomvar_names & need to account for models with more than 1 residual variance notation
      if(length(model$Residual$nfl>0)) {
        number_random=sum(model$Random$nfl)+sum(model$Residual$nfl) #add up random variances and residual variances
        colnames(t_icc)<-randomvar_names[seq(from=i,to=number_random,by=length(responses))]
        icc_all<-cbind(icc_all,t_icc)
      } else  {
        colnames(t_icc)<-randomvar_names[seq(from=i,to=i+length(responses)*length(model$Random$nfl),by=length(responses))]
        icc_all<-cbind(icc_all,t_icc)
      }             
    } 
    
    icc_all<-as.mcmc(icc_all[,-1])
    
  } else {
    #Single response models   
    tvar<-var_terms
    colnames(tvar)<-variances
    tsumvar<-rowSums(tvar) #calculate sum of variances
    tsumvar<-tsumvar + dist_var[1] #Add distribution variance
    tot_tsumvar<-tsumvar+S2var #Add sampling variance
    
    if(any(grepl("animal",colnames(tvar)))) {
      t_icc<-(tvar/tot_tsumvar)*100
      #If a phylogeny is included then divide all terms apart from phylogeny by total variance. ICC of phylogeny should exclude sampling variance as phylogenetic relatedness is a "fixed random effect": see Nakagawa & Santos 2012
    } else  {
      t_icc<-tvar
      t_icc[,grepl("animal",colnames(t_icc))]<-(t_icc[,grepl("animal",colnames(t_icc))]/(tsumvar))*100
      t_icc[,!grepl("animal",colnames(t_icc))]<-(t_icc[,!grepl("animal",colnames(t_icc))]/(tot_tsumvar))*100
    }
    icc_all<-t_icc
    colnames(icc_all)<-randomvar_names#Need to rename colnames to match randomvar_names
  }             
  
  icc_all<-icc_all[,colnames(var_terms)]#Reorder to match variances
  
  #If there is sampling variance then add ICC calculations
  if(sum(S2var) ==0) {
    icc_S2var<-NULL
  } else  {
    icc_S2var<-numeric()
    icc_S2var<-c(icc_S2var,round(((S2var[i]/mean(tot_tsumvar))*100),dec_PM))
    names(icc_S2var)[i]<-paste("Sampling variance",responses[i])
  }  
  
  #Summaries
  icc1=paste(round(colMeans(icc_all),dec_PM)," (",round(HPDinterval(icc_all)[,1],dec_PM), ", ",round(HPDinterval(icc_all)[,2],dec_PM),")",sep="")

  #Format Random effects output
  if(sum(S2var) ==0){
    randomVar<-data.frame("Random Effects"=c(colnames(var_terms)),"Posterior Mode (CI)"=c(rand1),"I2 % (CI)"=c(icc1), check.names=FALSE)
  } else  {
    randomVar<-data.frame("Random Effects"=c(colnames(var_terms),names(icc_S2var)),"Posterior Mode (CI)"=c(rand1,round(S2var,dec_PM)),"I2 % (CI)"=c(icc1,icc_S2var), check.names=FALSE)
  }
  
}
  

  #****************************************************
  #Excel output: random effects ----
  #****************************************************
  writeData(workbook, sheet, randomVar, startCol = 1, startRow = row_nums,headerStyle = hs2)
  row_nums = row_nums + dim(randomVar)[1]+2
  
  }
  
  #****************************************************
  #Correlations ----
  #****************************************************
  if(Include_random == "no") {
    #If no random effects are estimated then skip correlation estimation
  } else  {
    
  #If no covariances skip this part
  if(is.null(randomcovar_names)) {
  } else  {
    #Covariance identification if not specified
    if(is.null(covariances)) {
      covariances<-colnames(model$VCV)[colnames(model$VCV) != variances]
    } else  {
      covariances<-covariances
    }
    #rename phylogeny term to animal
    covariances<-sub(ginv, "animal", covariances)
    covar_terms<-as.mcmc(as.matrix((model$VCV[,covariances])))
    colnames(covar_terms)<-randomcovar_names
    
    #rename covariances effects
    if(is.null(randomcovar_names)) {
      colnames(covar_terms)<- colnames(covar_terms)
    } else  {
      colnames(covar_terms)<-randomcovar_names
    }
    
    #Need to pull out variances relating to covariances and cycle through covariance terms to calculate correlations
    corrs=matrix(0, nrow = dim(covar_terms)[1], ncol = 1)
    for(i in 1:dim(covar_terms)[2]) {
      #covar
      cov<-covar_terms[,i]
      #find variances 
      cov_vars=unlist(strsplit(colnames(covar_terms)[i], ":"))
      #remove punctuation and find matches
      var_names<-amatch(gsub("[[:punct:][:blank:]]+", "", cov_vars),gsub("[[:punct:][:blank:]]+", "", colnames(var_terms)),maxDist=2) 
      #Check variances are extracted correctly
      if(sum(is.na(var_names))>0) stop("covariance term needs to match variances separated by :")
      vars<-var_terms[,var_names] 
      #Calculate corrs
      tmpcor<-cov/sqrt(vars[,1]*vars[,2])
      corrs<-cbind(corrs,tmpcor)
    }
    corrs<-as.mcmc(corrs[,-1])
    #Cor Summaries
    cor1=paste(round(posterior.mode(corrs),dec_PM)," (",round(HPDinterval(corrs)[,1],dec_PM), ", ",round(HPDinterval(corrs)[,2],dec_PM),")",sep="")
    ncors<-ifelse(is.null(dim(corrs)), 1,dim(corrs)[2])
    nits<-ifelse(is.null(dim(corrs)),length(corrs), dim(corrs)[1])
    if(ncors >1){
      pCor=pmax(0.5/nits, pmin(colSums(corrs[,1:ncors, drop = FALSE] > 0)/nits, 1 - colSums(corrs[, 1:ncors, drop = FALSE] > 0)/nits))*2
    } else  {
      pCor=pmax(0.5/nits, pmin(sum(corrs[,drop = FALSE] > 0)/nits, 1 - sum(corrs[, drop = FALSE] > 0)/nits))*2
    }
    randomCorr<-data.frame("Correlations"=colnames(covar_terms),"Posterior Mode (CI)"=cor1,"pMCMC"=round(as.numeric(pCor),3), check.names=FALSE)
    
    ##Write data to excel sheet ----
    writeData(workbook, sheet, randomCorr, startCol = 1, startRow = row_nums,headerStyle = hs2)
    conditionalFormatting(workbook, sheet, cols=3, rows=row_nums+4:10000, rule="<0.05", style = bolding)
    row_nums = row_nums + dim(randomCorr)[1]+2
  }
  
  ##Cor_diffs ----
  #If no cor_diffs then skip this part
  if(is.null(cor_diffs)) {
    return(workbook)
  } else  {
    #Calculate differences between correlations
    colnames(corrs)=colnames(covar_terms)
    corr_comp<-pairwise.diffs(corrs,nF=length(colnames(corrs)))
    
    #Select specified comparisons
    corr_select = corr_comp %>% dplyr::filter(Fixed_Effects %in% cor_diffs == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC) %>% dplyr::rename("Correlation comparions"="Fixed_Effects")
    
    ##Write data to excel sheet
    writeData(workbook, sheet, corr_select, startCol = 1, startRow = row_nums,headerStyle = hs2)
    conditionalFormatting(workbook, sheet, cols=3, rows=row_nums+4:10000, rule="<0.05", style = bolding)
    row_nums = row_nums + dim(corr_select)[1]+2
  }
}  
  #Some final formatting ----
  for(i in 1:length(sheets(workbook))){
    setColWidths(workbook, sheet = i, cols = 1:20, widths = "auto")
  }
  return(workbook)
}

#******************************************************
#function for extracting df from xl workbook
xl_2_df = function(xltab,sheet=NULL){
  df<-readWorkbook(xltab,sheet=sheet,startRow = 2)
  colnames(df)<-gsub("[.]"," ",colnames(df))
  rownames(df)<-NULL
  return(df)
}

xl_2_df2 = function(xltab,sheet=NULL){
  df<-readWorkbook(xltab,sheet=sheet)
  colnames(df)<-gsub("[.]"," ",colnames(df))
  rownames(df)<-NULL
  return(df)
}

#******************************************************
#function for df to Rmd table
md_table = function(df){
  pacman::p_load(kableExtra)
  if(length(grep("Random",df$`Fixed Effects`))>0 & length(grep("Correlations",df$`Fixed Effects`))>0) {
    Fixed<-df[1:(as.numeric(row.names(df[grepl("Random",df$`Fixed Effects`),]))-1),]
    Random<-df[as.numeric(row.names(df[grepl("Random",df$`Fixed Effects`),])):(as.numeric(row.names(df[grepl("Correlations",df$`Fixed Effects`),]))-1),]
    Corrs<-df[as.numeric(row.names(df[grepl("Correlations",df$`Fixed Effects`),])):dim(df)[1],]
    rows_bold<-c(ifelse(Fixed$pMCMC<0.05,T,F),ifelse(Random$pMCMC<0.05,F,F),ifelse(Corrs$pMCMC<0.05,T,F))
  } else  {if(length(grep("Random",df$`Fixed Effects`))>0) {
    Fixed<-df[1:(as.numeric(row.names(df[grepl("Random",df$`Fixed Effects`),]))-1),]
    Random<-df[as.numeric(row.names(df[grepl("Random",df$`Fixed Effects`),])):dim(df)[1],]
    rows_bold<-c(ifelse(Fixed$pMCMC<0.05,T,F),ifelse(Random$pMCMC<0.05,F,F))
  } else  {
    Fixed<-df
    rows_bold<-ifelse(Fixed$pMCMC<0.05,T,F)}
  }
  kbl(df, align = "l", digits = 3) %>%
    kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) %>%
    row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Fixed Effect Comparisons",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Random",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Correlations",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Correlation comparions",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") %>%
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px") %>%
    column_spec(column=3, bold =rows_bold)}

#******************************************************
md_table2 = function(df){
  pacman::p_load(kableExtra)
  kbl(df, align = "l", digits = 3) %>%
    kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) %>%
    row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Random",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Correlations",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Correlation comparions",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px")}

