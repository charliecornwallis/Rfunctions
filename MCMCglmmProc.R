########################################
#Function for processing MCMCglmm models
########################################
#MCMCglmmProc<-function(model=NULL,S2var=0,start_row=1,workbook=NULL, create_sheet="yes",sheet="Results",title=NULL,fixed_names=NULL,fixed_del="none",fixed_grp=NULL,fixed_diffdel="none",fixed_diffinc=NULL,variances=NULL,covariances=NULL,randomvar_names=NULL,randomcovar_names=NULL,Include_random = "yes",padding=4)

#explanation:
#model = MCMCglmm model
#Does the model have a binary response variable? Changes estimation of random effects and residual variance = fixed
#S2var = sampling variance if known - useful for meta-analyses
#start_row=starting row of workbook to add data to if NULL put data in first empty row 
#workbook = adds data if specified, otherwise will make a new e.g. Results
#create_sheet = should a new sheet be created e.g."yes" vs "no"
#sheet= name of sheet "Analysis 1"
#title = Title of table in e.g. "Table 1"
#fixed_names = what you want fixed effects to be called e.g c("Intercept","Season length")
#fixed_del = any fixed effects that should deleted from output - useful if only assessing higher order interactions from a model. Also note that if this is specified then column headers will be suppressed in output table
#fixed_grp = vector that specifies which differences between fixed effects should be calculated. Needs to be the same length as fixed_names. If not included then all differences will be calculated. e.g. c(1,1,1,2,2,2)
#fixed_diffdel = specify comparisons between fixed effects that should be removed from output e.g. c("effect1 vs effect2","effect3 vs effect4"). Must exactly match names of fixed effects and be separate by " vs "
#fixed_diffinc = same as fixed_diffdel but for specific terms to be included in output
#fixed_diff_diffs = calculates differences between differences e.g. c("effect1 vs effect2 - effect3 vs effect4"). Must exactly match names of fixed effects and be separate by " - "
#variances = names of variance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed
#covariances = names of covariance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed
#randomvar_names = what you want random effect variances to be called. Note reserved term animal is renamed to phylogeny
#randomcovar_names = what you want random effect covariances to be called. Note for correlations to be calculated this has to be same as variance1 : variance2
#Include_random = should random effects be included in output
#Padding = space between tables when outputting multiple models to same sheet
#dec_PM = number of decimals given for posterior mode and CIs of fixed and random effects
#link = allows logit & probit, the default is gaussian, and is used to calculate ICCs correctly. Gaussian can also be used for poisson (log) models to produce ICCs on expected scale but NOT for data scale estimates. See de Villemereuil 2016 Genetics & QGglmm package for details.
#For Multi-response models can provide a list of link functions (e.g. c("gaussian","logit")) corresponding to each response trait
#responses - specify response variables can take multiple values for multi response
MCMCglmmProc<-function(model=NULL,responses=NULL,link=c("gaussian"),S2var=0,start_row=1,workbook=NULL, create_sheet="yes",sheet="Results",title=NULL,fixed_names=NULL,fixed_del="none",fixed_grp=NULL,fixed_diffdel="none",fixed_diffinc="all",fixed_diff_diffs =NULL,variances=NULL,covariances=NULL,randomvar_names=NULL,randomcovar_names=NULL,Include_random = "yes",padding=4,dec_PM=2)
{ #Load packages
  pacman::p_load(MCMCglmm,coda)
  
  #Remove unwanted effects
  #Remove mev from random effect
  model$VCV<-model$VCV[, colnames(model$VCV) !="sqrt(mev):sqrt(mev).meta"]
  
  #rename model fixed effects
  if(is.null(fixed_names)) {
    colnames(model$Sol)<- colnames(model$Sol)
  } else  {colnames(model$Sol)<-fixed_names
  }
  
  #if variances and covariances are not specified make them the same as random names
  if(is.null(variances)) {
    variances<-randomvar_names
  } else  {variances<-variances
  }
  if(is.null(covariances)) {
    covariances<-randomcovar_names
  } else  {covariances<-covariances
  }
  #****************************************************
  #Fixed effects and Differences between levels
  #****************************************************
  #Main effects
  nF=model$Fixed$nfl
  fe1=paste(round(posterior.mode(model$Sol),dec_PM)," (",round(HPDinterval(model$Sol)[,1],dec_PM), " , ",round(HPDinterval(model$Sol)[,2],dec_PM),")",sep="")
  
  #P values using summary.MCMCglmm code
  fe1_p=pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:nF, drop = FALSE] > 0)/dim(model$Sol)[1], 1 - colSums(model$Sol[, 1:nF, drop = FALSE] > 0)/dim(model$Sol)[1]))*2
  fe1=data.frame(Fixed_Effects=colnames(model$Sol),Estimates=fe1, pMCMC=fe1_p)
  
  #Do any fixed effects need to be deleted
  if(any(fixed_del == "none")) {
    fe1=fe1
  } else  {
    fe1 = fe1 %>% dplyr::filter(Fixed_Effects %in% fixed_del == F) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  #function for calculating differences between all columns of a matrix
  pairwise.diffs <- function(x, nF=1)
  {if(is.matrix(x) & nF>1) {
    #Differences
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
    fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), " , ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
    fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
    fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=fe2_p, check.names=FALSE)
    return(fe2)
  } 
  }
  
  #Differences between fixed effects
  fe2<-pairwise.diffs(model$Sol,nF=model$Fixed$nfl)
  
  #Estimate differences between specified groups of fixed effects
  if(is.null(fixed_grp)) {
    fixed=rbind(fe1,fe2)
  } else  {
    tmp_mat<-matrix(1,nrow=dim(model$Sol)[1],ncol=dim(model$Sol)[2])
    colnames(tmp_mat)<-fixed_grp
    combinations<-pairwise.diffs(x=tmp_mat,nF=2)
    combinations<-combinations %>% dplyr::select(Fixed_Effects) %>% separate(Fixed_Effects,c("comb1","comb2"),sep =" vs ")
    #dplyr::select desired combinations
    fe2<-cbind(fe2,combinations)
    fe2<-fe2 %>% dplyr::filter(as.numeric(comb1) == as.numeric(comb2)) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    #Combine main effects and differences
    fixed=rbind(fe1,fe2)
  }

  #Do any comparisons need to be deleted
  if(any(fixed_diffdel == "none")) {
    fixed=fixed
  } else  {
    fixed = fixed %>% dplyr::filter(Fixed_Effects %in% fixed_diffdel == F) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  #Should only specific effects be included
  if(any(fixed_diffinc == "all")) {
    fixed=fixed
  } else  {
    fixed = fixed %>% dplyr::filter(Fixed_Effects %in% c(fixed_names,fixed_diffinc) == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
  
  #Should all comparison be deleted
  if(any(fixed_diffdel == "all")) {
    fixed=fe1
  } else  {
    fixed = fixed
  }
  
  #Should any differences of differences be calculated
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
      fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), " , ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
      fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
      fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=fe2_p, check.names=FALSE)
      return(fe2)
    } 
    }
    diffs.diffs<-pairwise.diffs2(diffs.mat)
    #Select differences that are specified
    diffs.diffs = diffs.diffs %>% dplyr::filter(Fixed_Effects %in% fixed_diff_diffs == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    #add results to fixed effects
    fixed = rbind(fixed,diffs.diffs)
  }
  
  #****************************************************
  #Random effects
  #****************************************************
  #Separate out variance and covariance terms
  var_terms<-model$VCV[,variances]
  
  #rename random effects if specified
  if(is.null(randomvar_names)) {
    colnames(var_terms)<- colnames(var_terms)
    randomvar_names<-variances
  } else  {colnames(var_terms)<-randomvar_names
  }
  
  #****************************************************
  #Variance output
  #Summary of variance components
  rand1=paste(round(posterior.mode(var_terms),dec_PM)," (",round(HPDinterval(var_terms)[,1],dec_PM), " , ",round(HPDinterval(var_terms)[,2],dec_PM),")",sep="")
  
  #Calculate % variation explained by each random effect
  
  #Distribution variances
  link_var<-link
  link_var[link_var=="gaussian"]<-0
  link_var[link_var=="poisson"]<-0
  link_var[link_var=="logit"]<-pi^2/3
  link_var[link_var=="probit"]<-1
  
  #Setup sampling variances for multi-response models
  if(length(responses)>1 & sum(S2var) ==0){
    S2var<-rep(0,length(responses))
  } else  {S2var<-S2var
  }
  
  #Separate variances from each trait for multi-response models and create a var_sum with trait specific values
  icc_all<-mcmc()
  icc_S2var<-numeric()
  
  if(length(responses)>1){
    
    for(i in 1:length(responses)){
      tvar<-var_terms
      colnames(tvar)<-variances
      tvar<-tvar[,grepl(responses[i],colnames(tvar))]
      tsumvar<-rowSums(tvar) #Calculate sum of variances
      tsumvar<-tsumvar + as.numeric(link_var[i]) #Add distribution variance
      tot_tsumvar<-tsumvar+S2var #Add sampling variance
      
      if(any(grepl("animal",colnames(tvar)))) {
        t_icc<-(tvar/tot_tsumvar)*100
        #If a phylogeny is included then divide all terms apart from phylogeny by total variance. ICC of phylogeny should exclude sampling variance as phylogenetic relatedness is a "fixed random effect": see Nakagawa & Santos 2012
      } else  {
        t_icc<-tvar
        t_icc[,grepl("animal",colnames(t_icc))]<-(t_icc[,grepl("animal",colnames(t_icc))]/(tsumvar))*100
        t_icc[,!grepl("animal",colnames(t_icc))]<-(t_icc[,!grepl("animal",colnames(t_icc))]/(tot_tsumvar))*100
      }
      colnames(t_icc)<-randomvar_names[c(i,(i+(length(randomvar_names)/(1+length(model$Random$nrt)))))]#Need to rename colnames to match randomvar_names
      icc_all<-cbind(icc_all,t_icc)
      icc_S2var<-c(icc_S2var,round(((S2var[i]/mean(tot_tsumvar))*100),dec_PM))
      names(icc_S2var)[i]<-paste("Sampling variance",responses[i])
    }             
    icc_all<-as.mcmc(icc_all[,-1])
    
  } else {
    
    tvar<-var_terms
    colnames(tvar)<-variances
    tsumvar<-rowSums(tvar) #calculate sum of variances
    tsumvar<-tsumvar + as.numeric(link_var[1]) #Add distribution variance
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
    icc_S2var<-c(icc_S2var,round(((S2var[1]/mean(tot_tsumvar))*100),dec_PM))
    names(icc_S2var)[1]<-paste("Sampling variance",responses[1])
  }             
  
  icc_all<-icc_all[,colnames(var_terms)]#Reorder to match variances
  
  #Summaries
  icc1=paste(round(colMeans(icc_all),dec_PM)," (",round(HPDinterval(icc_all)[,1],dec_PM), " , ",round(HPDinterval(icc_all)[,2],dec_PM),")",sep="")
  
  #Output to excel
  fixed<-data.frame("Fixed Effects"=fixed$Fixed_Effects,"Posterior Mode (CI)"=fixed$Estimates,"pMCMC"=round(as.numeric(fixed$pMCMC),3),check.names=FALSE)
  randomVar<-data.frame("Random Effects"=c(colnames(var_terms),names(icc_S2var)),"Posterior Mode (CI)"=c(rand1,round(S2var,dec_PM)),"I2 % (CI)"=c(icc1,icc_S2var), check.names=FALSE)
  
  #Remove sampling variances if they weren't specified
  randomVar<-randomVar[!grepl("Sampling variance",randomVar$`Random Effects`) | randomVar$`Posterior Mode (CI)`!=0,]
  
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
    start_row=dim(readWorkbook(workbook,sheet =sheet,skipEmptyRows = F))[1]+padding
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
  #Remove column headings if deleting fixed effects as it will be assessing higher order interactions where column names are not needed
  if(any(fixed_del == "none")) {
    writeData(workbook, sheet, fixed, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
  } else  {
    writeData(workbook, sheet, fixed, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2,colNames =FALSE)
  }
  
  #Bold pMCMC values less than 0.05
  bolding<-createStyle(textDecoration="bold")
  conditionalFormatting(workbook, sheet, cols=3, rows=1:10000, rule="<0.05", style = bolding)
  
  #Random effects: variances
  #Should they be outputted or not
  if(Include_random == "yes") {
    writeData(workbook, sheet, randomVar, startCol = 1, startRow = start_row+dim(header)[1]+dim(fixed)[1]+1,headerStyle = hs2)
  } else  {
    workbook=workbook
  }
  
  #Random effects: correlations
  #****************************************************
  #Covariance output - correlations
  
  #If no covariances skip this part
  if(is.null(covariances)) {
    return(workbook)
  }
  else  {
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
    cor1=paste(round(posterior.mode(corrs),dec_PM)," (",round(HPDinterval(corrs)[,1],dec_PM), " , ",round(HPDinterval(corrs)[,2],dec_PM),")",sep="")
    ncors<-ifelse(is.null(dim(corrs)), 1,dim(corrs)[2])
    nits<-ifelse(is.null(dim(corrs)),length(corrs), dim(corrs)[1])
    if(ncors >1){
      pCor=pmax(0.5/nits, pmin(colSums(corrs[,1:ncors, drop = FALSE] > 0)/nits, 1 - colSums(corrs[, 1:ncors, drop = FALSE] > 0)/nits))*2
    } else  {
      pCor=pmax(0.5/nits, pmin(sum(corrs[,drop = FALSE] > 0)/nits, 1 - sum(corrs[, drop = FALSE] > 0)/nits))*2
    }
    randomCorr<-data.frame("Correlations"=colnames(covar_terms),"Posterior Mode (CI)"=cor1,"pMCMC"=round(as.numeric(pCor),3), check.names=FALSE)
    
    #Write data to excel sheet
    writeData(workbook, sheet, randomCorr, startCol = 1, startRow = start_row+dim(header)[1]+dim(fixed)[1]+dim(randomVar)[1]+padding,headerStyle = hs2)
    conditionalFormatting(workbook, sheet, cols=3, rows=start_row+dim(header)[1]+dim(fixed)[1]+dim(randomVar)[1]+4:10000, rule="<0.05", style = bolding)
  }
  
  #Some formatting
  for(i in 1:length(sheets(workbook))){
    setColWidths(workbook, sheet = i, cols = 1:20, widths = "auto")
  }
  #return(workbook)
  return(workbook)
}
