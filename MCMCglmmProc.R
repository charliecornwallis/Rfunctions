#===========================================================
#Function for processing MCMCglmm models#
#===========================================================

#Trouble shooting tools----
# model=m
# responses=NULL
# dist_var=c(0)
# fixed_names=NULL
# Include_random = "yes"
# variances=
# covariances =
# randomvar_names= 
# cor_vcvs =
# cor_names =
# cor_select = NULL
# cor_diffs = 
# partial_vcvs = NULL 
# partial_names = NULL
# partial_select = NULL
# partial_diffs = NULL
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

MCMCglmmProc<-function(start_row=NULL,workbook=NULL, create_sheet="yes",sheet="sheet1",title="", padding=4,dec_PM=2, #specify output details
                       model=NULL,responses=NULL,dist_var=NULL,ginv="animal",S2var=0,#specify model details
                       fixed_names=NULL,fixed_diffinc="none",fixed_diff_diffs =NULL, pvalues="include", #fixed effect specifications
                       Include_random = "no", variances=NULL,randomvar_names=NULL, #random effect specifications
                       cor_vcvs=NULL,cor_names=NULL,cor_select=NULL, cor_diffs=NULL,#covariance specifications
                       partial_vcvs = NULL, partial_names = NULL, partial_select=NULL, partial_diffs = NULL) #partial correlation specifications
{ 
  #Explanation of terms ----
  #model = MCMCglmm model
  #responses = specify response variables can take multiple values for multi response
  #dist_var = distribution variance associated with link function: e.g."gaussian" = 0, "log" = log(1 + log(exp(intercept + 0.5*sumRE)), "logit" = pi^2/3, "probit" = 1. This is a complicated issue and should be given careful thought: see Nakagawa et al 2017. #For Multi-response models can provide a list of dist_var and S2var corresponding to each response trait
  #S2var = sampling variance if known - useful for meta-analyses
  #start_row=starting row of workbook to add data to if NULL put data in first empty row 
  #workbook = adds data if specified, otherwise will make a new e.g. Results
  #create_sheet = should a new sheet be created e.g."yes" vs "no"
  #sheet= name of sheet "Analysis 1"
  #title = Title of table in e.g. "Table 1"
  #Padding = space between tables when outputting multiple models to same sheet
  #dec_PM = number of decimals given for posterior mode and CIs of fixed and random effects
  #pvalues = exclusion of pMCMC values for fixed effects - "exclude" = exclude all, "include" = include all or index pvalues to be excluded e.g. "c(1,3)" removes 1st and 3rd, c(1:7) removes 1 to 7. Note pMCMC will still be calculated for fixed effect comparisons.
  #fixed_names = what you want fixed effects to be called e.g c("Intercept","Season length")
  #fixed_diffinc = differences between fixed effects to be included in output. Terms that come first in fixed_names have to come first in the comparison.
  #fixed_diff_diffs = calculates differences between differences e.g. c("effect1 vs effect2 - effect3 vs effect4"). Must exactly match names of fixed effects and be separate by " - "
  #variances = names of variance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed. If left as NULL then taken from the model object.
  #randomvar_names = what you want random effect variances to be called. Note reserved term animal is renamed to phylogeny
  #Include_random = should random effects be included in output
  #cor_vcvs = specify the variance covariance (vcv) matrices that you want to calculate correlations for as a list [e.g. list(c(1:9),c(10:18)) ]   
  #cor_names = name the traits for each set of correlations as a list [e.g. list(c("trait1 phylogeny","trait2 phylogeny","trait3 phylogeny"),c("trait1 residual","trait2 residual","trait3 residual")) ]
  #cor_select = names of the correlations that you want presented in tables
  #cor_diff = calculates differences between correlations. Should be specified in the same way as fixed_diffs e.g. c("cor1 vs cor2","cor1 vs cor3"...) 
  #partial_vcvs = specify the variance covariance (vcv) matrices that you want to calculate partial correlations for as a list [e.g. list(c(1:9),c(10:18)) ]   
  #partial_names = name the traits for each set of partial correlations as a list [e.g. list(c("trait1 phylogeny","trait2 phylogeny","trait3 phylogeny"),c("trait1 residual","trait2 residual","trait3 residual")) ]
  #partial_select = names of the partial correlations that you want presented in tables
  #partial_diffs = specify the partial correlations that you want to compare. Names must match the partial names separated by vs [e.g. c("trait2 phylogeny_trait1 phylogeny vs trait2 residual_trait1 residual") ]

  #Naming aid ----
  #variances=c("trait1:trait1.animal","trait2:trait2.animal","trait3:trait3.animal","trait4:trait4.animal","trait5:trait5.animal","trait6:trait6.animal","trait7:trait7.animal","trait8:trait8.animal","trait9:trait9.animal","trait1:trait1.units","trait2:trait2.units","trait3:trait3.units","trait4:trait4.units","trait5:trait5.units","trait6:trait6.units","trait7:trait7.units","trait8:trait8.units","trait9:trait9.units"),
  #randomvar_names=c("Phylogeny Trait1","Phylogeny Trait2","Phylogeny Trait3","Phylogeny Trait4","Phylogeny Trait5","Phylogeny Trait6","Phylogeny Trait7", "Phylogeny Trait8","Phylogeny Trait9","Residual Trait1","Residual Trait2","Residual Trait3","Residual Trait4","Residual Trait5","Residual Trait6","Residual Trait7","Residual Trait8","Residual Trait9"),
  #cor_vcvs = list(c(1:9),c(10:18))    
  #cor_names = list(c("trait1 phylogeny","trait2 phylogeny","trait3 phylogeny"),c("trait1 residual","trait2 residual","trait3 residual"))
  #cor_diffs = c("trait2 phylogeny : trait1 phylogeny vs trait2 residual : trait1 residual")
  #partial_vcvs = list(c(1:9),c(10:18))    
  #partial_names = list(c("trait1 phylogeny","trait2 phylogeny","trait3 phylogeny"),c("trait1 residual","trait2 residual","trait3 residual"))
  #partial_diffs = c("trait2 phylogeny : trait1 phylogeny vs trait2 residual : trait1 residual")

  #Load packages and naming ----
  pacman::p_load(MCMCglmm,coda,openxlsx,stringdist,kableExtra,corpcor,tidyverse)

  #===========================================================
  #Check terms are specified correctly and relabel terms ----
  #===========================================================
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
  if (any(colnames(model$VCV) =="sqrt(mev):sqrt(mev).meta")) {
    model$Random$nfl <- model$Random$nfl[1:length(model$Random$nfl)-1]
    model$Random$nrl <- model$Random$nrl[1:length(model$Random$nrl)-1]
    model$Random$nat <- model$Random$nat[1:length(model$Random$nat)-1]
    model$Random$nrt <- model$Random$nrt[1:length(model$Random$nrt)-1]

    model$VCV<-model$VCV[, colnames(model$VCV) !="sqrt(mev):sqrt(mev).meta"]

  } else {
  }

  #Rename model fixed effects
  if(is.null(fixed_names)) {
    fixed_names<- colnames(model$Sol)
  } else  {colnames(model$Sol)<-fixed_names
  }
  
  #===========================================================
  #Fixed effects and Differences between levels ----
  #===========================================================
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
    fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),5), check.names=FALSE)
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
  
    #Include all else specified differences  
  if(any(fixed_diffinc == "all")) {
    fixed = fixed %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  } else  {
    fixed = fixed %>% dplyr::filter(Fixed_Effects %in% c(fixed_names,fixed_diffinc) == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
  }
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
      fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),4), check.names=FALSE)
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
  fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=as.numeric(fixeddiff$pMCMC),check.names=FALSE)
  
  #===========================================================
  #Excel output: fixed effects ----
  #===========================================================
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

  #===========================================================
  #Random effects ----
  #===========================================================
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
 
  #warn if variances are not specified
  if (is.null(variances)) {
    stop("variances must be specified either as indices of VCV or column names of VCV")
  } else {
  }

  #if numeric string submitted then pull out names, otherwise leave as is
  if(is.numeric(variances)) {
    variances = colnames(model$VCV)[variances]
    } else  {
  }  
  
  #Rename phylogeny or pedigree term (ginv=) to animal
  colnames(model$VCV) <- sub(ginv, "animal", colnames(model$VCV))
  variances <- sub(ginv, "animal", variances)
  
  ##Separate out variance and covariance terms ----
  var_terms<-model$VCV[,variances]
  
  #rename random effects if specified
  if(is.null(randomvar_names)) {
    colnames(var_terms)<- variances
    randomvar_names<-variances
  } else  {colnames(var_terms)<-randomvar_names
  }
  
  #===========================================================
  ##Variance output ----
  #===========================================================
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
      if(length(model$Residual$nfl)>1) {
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
  
  #===========================================================
  #Excel output: random effects ----
  #===========================================================
  writeData(workbook, sheet, randomVar, startCol = 1, startRow = row_nums,headerStyle = hs2)
  row_nums = row_nums + dim(randomVar)[1]+2
  
  }
  
  #===========================================================
  #Correlations ----
  #===========================================================
  if(is.null(cor_vcvs)) {
  } else  {
    
    #Function for calculating partial correlations for each iteration from variance-covariance matrix
    MCMCcor = function (VCV, ntraits,trait_names) 
                  {result = data.frame(iteration = as.numeric(), trait_comb = as.character(), cor = as.numeric())
                      for (i in 1:dim(VCV)[1]) {
                          tmp = cov2cor(matrix(VCV[i, ], ntraits, ntraits)) #Calculate partial correlations
                          colnames(tmp) = trait_names #Name rows and columns
                          rownames(tmp) = trait_names
                          tmp <- data.frame(Row = rownames(tmp)[row(tmp)[lower.tri(tmp, diag = FALSE)]],
                                            Col = colnames(tmp)[col(tmp)[lower.tri(tmp, diag = FALSE)]],
                                            Value = tmp[lower.tri(tmp, diag = FALSE)]) #Place results in data frame
                        tmp <- data.frame(iteration = i,trait_comb = paste(tmp$Row, tmp$Col, sep = " : "), cor = tmp$Value) #Create trait combination names
                        result = rbind(result, tmp)
                      }
      
                      # transform to wide mcmc object
                      result = result %>% pivot_wider(id_cols = iteration, names_from = trait_comb, 
                                                        values_from = cor)
                      result = as.mcmc(result[, -1]) 
      
                      # calculate pMCMC values
                      ncors<-ifelse(is.null(dim(result)), 1,dim(result)[2])
                      nits<-ifelse(is.null(dim(result)),length(result), dim(result)[1])
      
                      if(ncors >1){
                        cor=pmax(0.5/nits, pmin(colSums(result[,1:ncors, drop = FALSE] > 0)/nits, 1 - colSums(result[, 1:ncors, drop = FALSE] > 0)/nits))*2
                      } else  {
                        cor=pmax(0.5/nits, pmin(sum(result[,drop = FALSE] > 0)/nits, 1 - sum(result[, drop = FALSE] > 0)/nits))*2
                      }
                      
                      # summarise mcmc samples
                      sum_res = data.frame(trait_comb = colnames(result), 
                                           post_mode_CI = paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep=""),
                                           pMCMC =round(cor,4))
                      
                    return(list(result,sum_res))
                    }
    
    #Create objects to write correlations to 
    corrs = mcmc()
    corrs_sum = data.frame()
    
    # extract vcvs, calculate correlations using function above and combine results from vcvs
    for(i in 1:length(cor_vcvs)) {
      # extract covariance matrix
      vcv = model$VCV[,cor_vcvs[[i]]]

      #calculate corrs for each iteration
      tmp = MCMCcor(VCV = vcv, ntraits = length(cor_names[[i]]),trait_names = cor_names[[i]])
      corrs =cbind(corrs,tmp[[1]])
      corrs_sum = rbind(corrs_sum,tmp[[2]])
      }
    
    #sort formatting 
    corrs = mcmc(corrs[,-1]) 
    corrs_sum = data.frame("Correlations"=corrs_sum$trait_comb,"Posterior Mode (CI)"=corrs_sum$post_mode_CI,"pMCMC"=corrs_sum$pMCMC, check.names=FALSE)
    
    #Select specified comparisons
    if(is.null(cor_select)) {
      corrs_sum
      } else  {
    corrs_sum = corrs_sum %>% dplyr::filter(Correlations %in% cor_select == T) 
    }

    ##Write data to excel sheet ----
    writeData(workbook, sheet, corrs_sum, startCol = 1, startRow = row_nums,headerStyle = hs2)
    conditionalFormatting(workbook, sheet, cols=3, rows=row_nums+4:10000, rule="<0.05", style = bolding)
    row_nums = row_nums + dim(corrs_sum)[1]+2
  }
     
  #===========================================================
  #Cor_diffs ----
  #===========================================================
  #If no cor_diffs then skip this part
  if(is.null(cor_diffs)) {
  } else  {

    #Calculate differences between correlations
    corr_comp<-pairwise.diffs(corrs,nF=length(colnames(corrs)))

    #Select specified comparisons
    corr_select = corr_comp %>% dplyr::filter(Fixed_Effects %in% cor_diffs == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC) %>% dplyr::rename("Correlation Comparisons"="Fixed_Effects", "Posterior Mode (CI)"="Estimates")
    
    ##Write data to excel sheet
    writeData(workbook, sheet, corr_select, startCol = 1, startRow = row_nums,headerStyle = hs2)
    conditionalFormatting(workbook, sheet, cols=3, rows=row_nums+4:10000, rule="<0.05", style = bolding)
    row_nums = row_nums + dim(corr_select)[1]+2
  }
 
  #===========================================================
  #Partial Correlations ----
  #===========================================================
  if(is.null(partial_vcvs)) {
  } else  {
    
    #Function for calculating partial correlations for each iteration from variance-covariance matrix
    MCMCpcor = function (VCV, ntraits,trait_names) 
                  {result = data.frame(iteration = as.numeric(), trait_comb = as.character(), partial_cor = as.numeric())
                      for (i in 1:dim(VCV)[1]) {
                          tmp = cor2pcor(matrix(VCV[i, ], ntraits, ntraits)) #Calculate partial correlations
                          colnames(tmp) = trait_names #Name rows and columns
                          rownames(tmp) = trait_names
                          tmp <- data.frame(Row = rownames(tmp)[row(tmp)[lower.tri(tmp, diag = FALSE)]],
                                            Col = colnames(tmp)[col(tmp)[lower.tri(tmp, diag = FALSE)]],
                                            Value = tmp[lower.tri(tmp, diag = FALSE)]) #Place results in data frame
                        tmp <- data.frame(iteration = i,trait_comb = paste(tmp$Row, tmp$Col, sep = " : "), partial_cor = tmp$Value) #Create trait combination names
                        result = rbind(result, tmp)
                      }
      
                      # transform to wide mcmc object
                      result = result %>% pivot_wider(id_cols = iteration, names_from = trait_comb, 
                                                        values_from = partial_cor)
                      result = as.mcmc(result[, -1]) 
      
                      # calculate pMCMC values
                      ncors<-ifelse(is.null(dim(result)), 1,dim(result)[2])
                      nits<-ifelse(is.null(dim(result)),length(result), dim(result)[1])
      
                      if(ncors >1){
                        pCor=pmax(0.5/nits, pmin(colSums(result[,1:ncors, drop = FALSE] > 0)/nits, 1 - colSums(result[, 1:ncors, drop = FALSE] > 0)/nits))*2
                      } else  {
                        pCor=pmax(0.5/nits, pmin(sum(result[,drop = FALSE] > 0)/nits, 1 - sum(result[, drop = FALSE] > 0)/nits))*2
                      }
                      
                      # summarise mcmc samples
                      sum_res = data.frame(trait_comb = colnames(result), 
                                           post_mode_CI = paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep=""),
                                           pMCMC =round(pCor,4))
                      
                    return(list(result,sum_res))
                    }
    
    #Create objects to write partial correlations to 
    pcorrs = mcmc()
    pcorrs_sum = data.frame()
    
    # extract vcvs, calculate partial correlations using function above and combine results from vcvs
    for(i in 1:length(partial_vcvs)) {
      # extract covariance matrix
      vcv = model$VCV[,partial_vcvs[[i]]]

      #calculate partial corrs for each iteration
      tmp = MCMCpcor(VCV = vcv, ntraits = length(partial_names[[i]]),trait_names = partial_names[[i]])
      pcorrs =cbind(pcorrs,tmp[[1]])
      pcorrs_sum = rbind(pcorrs_sum,tmp[[2]])
      }
    
    #sort formatting 
    pcorrs = mcmc(pcorrs[,-1]) 
    pcorrs_sum = data.frame("Partial Correlations"=pcorrs_sum$trait_comb,"Posterior Mode (CI)"=pcorrs_sum$post_mode_CI,"pMCMC"=pcorrs_sum$pMCMC, check.names=FALSE)
    
   #Select specified comparisons
    if(is.null(partial_select)) {
      pcorrs_sum
      } else  {
      pcorrs_sum = pcorrs_sum %>% dplyr::filter(`Partial Correlations` %in% partial_select == T) 
    }

    ##Write data to excel sheet ----
    writeData(workbook, sheet, pcorrs_sum, startCol = 1, startRow = row_nums,headerStyle = hs2)
    conditionalFormatting(workbook, sheet, cols=3, rows=row_nums+4:10000, rule="<0.05", style = bolding)
    row_nums = row_nums + dim(pcorrs_sum)[1]+2
  }
  
#===========================================================
#Partial  Cor_diffs ----
#===========================================================
#If no cor_diffs then skip this part
if(is.null(partial_diffs)) {
} else  {
  
  #Calculate differences between partial correlations
  pcorr_comp<-pairwise.diffs(pcorrs,nF=length(colnames(pcorrs)))
  
  #Select specified comparisons
  pcorr_select = pcorr_comp %>% dplyr::filter(Fixed_Effects %in% partial_diffs == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC) %>% dplyr::rename("Partial Correlation Comparisons"="Fixed_Effects","Posterior Mode (CI)"="Estimates")
  
  ##Write data to excel sheet
  writeData(workbook, sheet, pcorr_select, startCol = 1, startRow = row_nums,headerStyle = hs2)
  conditionalFormatting(workbook, sheet, cols=3, rows=row_nums+4:10000, rule="<0.05", style = bolding)
  row_nums = row_nums + dim(pcorr_select)[1]+2
}  

#===========================================================
#Some final formatting ----
#============================================================
for(i in 1:length(sheets(workbook))){
  setColWidths(workbook, sheet = i, cols = 1:20, widths = "auto")
}
return(workbook)
}

###############################################
#### Helper Functions for Table Formatting ####
###############################################

#===========================================================
#function for extracting df from xl workbook
#===========================================================
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

#===========================================================
#function for renaming sheets in xl workbook
#===========================================================
rename_xlsheets = function(wb,name,start_sheet=1) {
  # Get the names of all sheets
  sheet_names = names(wb)
  sheet_names = sheet_names[start_sheet:length(sheet_names)]
  
  # Loop through each sheet and rename it
  for (i in seq_along(sheet_names)) {
    new_name <- paste0(name, i)
    sheet_no = (start_sheet-1)+i
    names(wb)[[sheet_no]]=new_name
  }
}

#===========================================================
#function for df to md table that can handle html and other formats e.g. word ----
#===========================================================
md <- function(data,stats=FALSE) {
  pacman::p_load(flextable,officer)

##Output if presenting mixed model stats: bolding significant values and highlighting headings
  if(stats == TRUE){
        #bold rows if less than 0.05 excluding cells with brackets which will be random effect
        rows_bold = data |> mutate(signif = !grepl("\\(", pMCMC) & pMCMC < 0.05,
                               signif = ifelse(is.na(signif),"FALSE",signif)) |> pull(signif)
  } else {
    #Other tables 
  } 
  
  #Output if html format
  if (knitr::is_html_output()) {
  
  if(stats == TRUE){
    pacman::p_load(kableExtra)
      kbl(data, align = "l", digits = 3) |>
        kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) |>
        row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid")|>
        row_spec(grep("^Fixed Effect Comparisons",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")|>
        row_spec(grep("^Random",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")|>
        row_spec(grep("^Correlations",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")|>
        row_spec(grep("^Correlation Comparisons",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |>
        row_spec(grep("^Partial Correlations",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")|>
        row_spec(grep("^Partial Correlation Comparisons",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |>
        row_spec(nrow(data), extra_css = "border-bottom: 1px solid;margin-bottom:1000px") |>
        column_spec(column=which(names(data) == "pMCMC"), bold =rows_bold) |> 
        row_spec(1:nrow(data), extra_css = "height: 1em; white-space: nowrap;") |> 
        column_spec(1:ncol(data), width = "auto")
    }
    #Output if not mixed model output
    else  {  
      kbl(data,align = "l")  |> 
        kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11,fixed_thead = T) |>
        row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid") |>
        row_spec(nrow(data), extra_css = "border-bottom: 1px solid;margin-bottom:1000px") |>
        row_spec(1:nrow(data), extra_css = "height: 1em; white-space: nowrap;") |> 
        column_spec(1:ncol(data), width = "auto") |> 
        scroll_box(width = "1000px", height = "1000px")
    }
 }
  else  {
  #Output for other formats
  if(stats == TRUE){ 

  #Create table and highlight different sections
  ft <-  flextable(data)
   b <- fp_border(color = "black", width = 1)
      # Fixed Effect Comparisons
      fixedcomp_row <- grep("^Fixed Effect Comparisons", data[,1])
      ft  <- ft |> bold(i = fixedcomp_row, part = "body", bold = TRUE) |> 
                   bg(i = fixedcomp_row, part = "body", bg   = "#E7E5E5") |> 
                   hline(i = fixedcomp_row-1, part = "body", border = b) |> 
                   hline(i = fixedcomp_row, part = "body", border = b)

    # Random
    random_row <- grep("^Random", data[,1])
    ft  <- ft |> bold(i = random_row, part = "body", bold = TRUE) |> 
                            bg(i = random_row, part = "body", bg   = "#E7E5E5") |> 
                            hline(i = random_row-1, part = "body", border = b) |> 
                            hline(i = random_row, part = "body", border = b)
    # Correlations
    corr_row <- grep("^Correlations", data[,1])
    ft  <- ft |> bold(i = corr_row, part = "body", bold = TRUE) |> 
                                  bg(i = corr_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = corr_row-1, part = "body", border = b) |> 
                                  hline(i = corr_row, part = "body", border = b)
    # Correlation comparisons
    corrcomp_row <- grep("^Correlation Comparisons", data[,1])
    ft  <- ft |> bold(i = corrcomp_row, part = "body", bold = TRUE) |> 
                                  bg(i = corrcomp_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = corrcomp_row-1, part = "body", border = b) |> 
                                  hline(i = corrcomp_row, part = "body", border = b)
     # Partial Correlations
    pcorr_row <- grep("^Partial Correlations", data[,1])
    ft  <- ft |> bold(i = pcorr_row, part = "body", bold = TRUE) |> 
                                  bg(i = pcorr_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = pcorr_row-1, part = "body", border = b) |> 
                                  hline(i = pcorr_row, part = "body", border = b)
    # Correlation comparisons
    pcorrcomp_row <- grep("^Partial Correlation Comparisons", data[,1])
    ft  <- ft |> bold(i = pcorrcomp_row, part = "body", bold = TRUE) |> 
                                  bg(i = pcorrcomp_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = pcorrcomp_row-1, part = "body", border = b) |> 
                                  hline(i = pcorrcomp_row, part = "body", border = b)
    #Overall formatting
    ft <- ft |>
          theme_vanilla() |>
          flextable::fontsize(size = 10, part = "header") |>
          flextable::fontsize(size = 8, part = "body") |>
          bold(part = "header") |>
          #Bold pMCMC less than 0.05
          bold(i = which(rows_bold), j = 3, bold = TRUE, part = "body") |> 
          bg(part = "header", bg = "#E7E5E5") |>
          flextable::border(border.top = fp_border(color = "black", width = 1), part = "header") |>
          flextable::border(border.bottom = fp_border(color = "black", width = 1), part = "header") |>
          flextable::border(border.bottom = fp_border(color = "black", width = 1), i = nrow(data)) |>
          set_table_properties(width=1,layout = "autofit",opts_html = list(scroll = list(height = "1000px", freeze_first_column = TRUE)),opts_word = list(keep_with_next = TRUE))
          autofit(ft)
  }
    else {
    ft <- flextable(data) |>
          theme_vanilla() |>
          flextable::fontsize(size = 10, part = "header") |>
          flextable::fontsize(size = 8, part = "body") |>
          bold(part = "header") |>
          bg(part = "header", bg = "#E7E5E5") |>
          flextable::border(border.top = fp_border(color = "black", width = 1), part = "header") |>
          flextable::border(border.bottom = fp_border(color = "black", width = 1), part = "header") |>
          flextable::border(border.bottom = fp_border(color = "black", width = 1), i = nrow(data)) |>
          set_table_properties(width=1,layout = "autofit",opts_html = list(scroll = list(height = "1000px", freeze_first_column = TRUE)),opts_word = list(keep_with_next = TRUE))

      autofit(ft)
    }
    }
 }


#===========================================================
#function for df to Rmd table ----
#===========================================================
md_table = function(df){
  pacman::p_load(kableExtra)
  rows_bold = data |> mutate(signif = !grepl("\\(", pMCMC) & pMCMC < 0.05) |> pull(signif)
    kbl(df, align = "l", digits = 3) %>%    
    kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) %>%    
    row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Fixed Effect Comparisons",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Random",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Correlations",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Correlation Comparisons",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") %>% 
    row_spec(grep("^Partial Correlations",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Partial Correlation Comparisons",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") %>%   
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px") %>%    
    column_spec(column=3, bold =rows_bold)}

#===========================================================
## 2nd version
#===========================================================
md_table2 = function(df){pacman::p_load(kableExtra)  
    kbl(df, align = "l", digits = 3) %>%    
    kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) %>%    
    row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Random",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Correlations",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Correlation Comparisons",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%  
    row_spec(grep("^Partial Correlations",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(grep("^Partial Correlation Comparisons",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%    
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px")}
