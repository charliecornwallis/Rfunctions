#***************************************
#Function for processing Hmsc models#
#***************************************

HmscProc<-function(model=NULL,start_row=NULL,workbook=NULL, create_sheet="yes",sheet="sheet1",title="",fixed_names=NULL,fixed_diffinc="none",fixed_diff_diffs =NULL,fixed_diffinc_species="none",pvalues = "include",traits="exclude",trait_names=NULL,pvalues_traits= "exclude",VP_ave = "include",VPnames=NULL,randomvar_names=NULL,Include_random = "yes",Include_species ="exclude",pvalues_species="exclude", VP_species = "include",random_names_species=NULL,Include_random_species = "yes",community_comparisons = NULL,randomisation = "poisson",response_logged = FALSE,padding=4,dec_PM=2)
{ 
  #Explanation ----
  #1. Takes an Hmsc model and combines estimates from multiple chains and output 2 excel sheets: 1) averages across species; 2) Per species values. If there are multiple species these are averaged per mcmc sample using rowMeans to produce posterior distribution of average effects.
  #2. Estimates posterior modes and HPDintervals for effects
  #3. Calculates specified differences for fixed effects
  #4. If community_comparisons is specified (list of composition_metric,composition_variables, composition_comparisons (factors) and ngrid (# cut points for continuous variables default is 5) then it will calculate community composition differences and test if they are different from randomised data. Random communities are created in two possible ways: randomisation = "poisson" creates communities by sampling each taxa from a Poisson distribution with mean = mean taxa presence / abundance.response_logged=TRUE is used if log counts are modelled. Randomisation != poisson reruns model using randomised occurrence matrix (warning this can be slow).
  #5. For models with trait data, combines runs and estimates trait effects
  #6. Calculates average variance explained by random effects
  #7. Calculates % variation of explained by fixed and random effects 
  #8. For models with phylogenetic effects it outputs Rho
  #9. Calculates fit statistics
  #10. Output results to excel file
  
  #Example Terms
  # model=hm1 # model object
  # start_row=NULL #start row of excel sheet for species averages
  # workbook=NULL #name of workbook to add results to if already exists
  # create_sheet="yes" #create a new sheet for species average values
  # sheet="Table 1" #what to name sheet
  # title="Life is good" #title of tables
  # fixed_names=c("A","B","C","D") #names of fixed effects
  # fixed_diffinc=c("A vs B") #the differences between fixed effects for species averages to include e.g. fixed_diffinc=c("A vs B")
  # fixed_diff_diffs =NULL #compare differences of differences between fixed effects of species averages e.g. fixed_diff_diffs=c("A vs B - C vs D")
  #pvalues = exclusion of pMCMC values for fixed effects - "exclude" = exclude all, "include" = include all or index pvalues to be excluded e.g. "c(1,3)" removes 1st and 3rd, c(1:7) removes 1 to 7. Note pMCMC will still be calculated for fixed effect comparisons.
  # traits="include" #include trait effects: this examines the influence of species traits on community composition (e.g. richness if presence / absence of each species)
  # pvalues_traits = as for pvalues, but in relation to traits. Default is exclude.
  # Include_random = "yes" #include random effect estimates or not
  # randomvar_names=c("R1","R2","R3") #names of random effects - will take from model object if not specified
  # VP_ave = "include" #include variance partitioning for species averages
  # VPnames = renaming of terms in Variance partitioning table
  # Include_species ="include" #should a separate sheet with species estimates be included?
  # fixed_diffinc_species="none" #differences between fixed effects to include for specific species e.g. fixed_diffinc_species=c(c("A: species1 vs B: species 1"))
  # VP_species = "include" #include variance partitioning for all species
  # pvalues_species = "exclude". Same as above but for species level estimates. Default is exclude.
  # Include_random_species = "yes" #include random effect variances for each species
  # random_names_species = NULL #names for random effects for each species
  # padding=4 #spacing in excel file
  # dec_PM=2) #decimal places of estimates
  # community_comparisons = NULL # can take a list of lists e.g. "community_comparisons = list(Factor1 = list(composition_metric = "jaccard",comp_names = c("A","B","C","D"),composition_comp = c("A vs B", "B vs C")), Continuous1 = list(composition_metric = "jaccard",comp_names = "A",ngrid = 5))". Any "composition_metric" in vegdist function of vegan package (e.g."bray" or "jaccard") is allowed.
  # composition_var #Factor1, Continuous1 above = name of variable to compare composition metrics as it appears in the model formula 
  # comp_names #the names of the variable as you want to appear in the table. For factors in must include all levels even if a level is not included in the comparisons (e.g. "D" above).

  #****************************************************
  #Section 1: averages per species ----
  #****************************************************
  
  #Load packages and naming
  pacman::p_load(Hmsc,coda,stringdist,matrixStats,openxlsx,vegan,reshape2)
  
  #Convert model object
  post_model = convertToCodaObject(model)
  
  #****************************************************
  #Check terms are specified correctly ----
  #****************************************************
  if (is.null(fixed_names) &  any(fixed_diffinc != "none")) {
    stop("fixed names needs to be specified for fixed_diffinc to be calculated")
  } else {
  }
  
  #****************************************************
  #Fixed effects ----
  #****************************************************
  #combine chains
  tmp1 <- do.call(rbind, post_model$Beta)

  #If multiple species average across
  fixed_mod<-matrix(nrow = nrow(tmp1),ncol = (ncol(tmp1)/model$ns))
  
  #extract estimates and average estimates
  if(model$ns >1) {
  for(i in 1:model$nc){
    tmp2<-tmp1[,c(seq(from = i, to = model$ns*model$nc, by = model$nc))]
    fixed_mod[,i]<-rowMeans(tmp2)
  } 
  }   else  {
    fixed_mod = tmp1
  }
  
  fixed_mod = as.mcmc(fixed_mod)
  
  #Rename model fixed effects
  if(is.null(fixed_names)) {
    colnames(fixed_mod)<-model$covNames
  } else  {colnames(fixed_mod)<-fixed_names
  }
  
  #P values using summary.MCMCglmm code
  nF=dim(fixed_mod)[2]
  #Pvalues = option to exclude 
  if(any(pvalues == "include")) {
    fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nF, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nF, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
    fe1_p=round(as.numeric(fe1_p),3)
  } else  {
    if(any(pvalues == "exclude")) {
      fe1_p=rep("-",length(fixed_names))
    } else  {
      fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nF, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nF, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
      fe1_p=round(as.numeric(fe1_p),3)
      fe1_p[pvalues]<-"-"
    }
  }
  
  fe1=paste(round(posterior.mode(fixed_mod),dec_PM)," (",round(HPDinterval(fixed_mod)[,1],dec_PM), ", ",round(HPDinterval(fixed_mod)[,2],dec_PM),")",sep="")
  fe1=data.frame(Fixed_Effects=colnames(fixed_mod),Estimates=fe1, pMCMC=fe1_p)
  
  #****************************************************
  ##Differences between fixed effects ----
  #****************************************************
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
    fe2=paste(round(posterior.mode(result),dec_PM)," (",round(HPDinterval(result)[,1],dec_PM), ", ",round(HPDinterval(result)[,2],dec_PM),")",sep="")
    fe2_p=pmax(0.5/dim(result)[1], pmin(colSums(result[,1:ndiffs, drop = FALSE] > 0)/dim(result)[1], 1 - colSums(result[, 1:ndiffs, drop = FALSE] > 0)/dim(result)[1]))*2
    fe2=data.frame(Fixed_Effects=colnames(result),Estimates=fe2, pMCMC=round(as.numeric(fe2_p),3), check.names=FALSE)
    return(fe2)
  } 
  }
  
  #All differences between fixed effects
  fe2<-pairwise.diffs(fixed_mod,nF=nF)
  
  ##fixed_diffinc ----
  if(any(fixed_diffinc == "none")) {
    fixed=fe1
  } else  {
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
    diffs.mat<-pairwise.diffs.mat(fixed_mod)
    
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
  
  #****************************************************
  #Trait effects ----
  #****************************************************
  #combine chains
  gamma = do.call(rbind, post_model$Gamma)

  #Rename model trait names
  if(is.null(trait_names)) {
  } else  {model$trNames<-trait_names
  }
  
  #name columns more clearly
  gamma_names = character()
  for(i in 1:length(model$trNames)){
    tmp = paste(model$trNames[i],fixed_names,sep="_")
    gamma_names = c(gamma_names,tmp)
  } 
  
  colnames(gamma) = gamma_names
  gamma = as.mcmc(gamma)
  
  #P values using summary.MCMCglmm code
  nG=dim(gamma)[2]
  #Pvalues = option to exclude 
  if(any(pvalues_traits == "include")) {
    ge1_p=pmax(0.5/dim(gamma)[1], pmin(colSums(gamma[,1:nG, drop = FALSE] > 0)/dim(gamma)[1], 1 - colSums(gamma[, 1:nG, drop = FALSE] > 0)/dim(gamma)[1]))*2
    ge1_p=round(as.numeric(ge1_p),3)
  } else  {
    if(any(pvalues_traits == "exclude")) {
      ge1_p=rep("-",length(gamma_names))
    } else  {
      ge1_p=pmax(0.5/dim(gamma)[1], pmin(colSums(gamma[,1:nG, drop = FALSE] > 0)/dim(gamma)[1], 1 - colSums(gamma[, 1:nG, drop = FALSE] > 0)/dim(gamma)[1]))*2
      ge1_p=round(as.numeric(ge1_p),3)
      ge1_p[pvalues_traits]<-"-"
    }
  }
  
  ge1=paste(round(posterior.mode(gamma),dec_PM)," (",round(HPDinterval(gamma)[,1],dec_PM), ", ",round(HPDinterval(gamma)[,2],dec_PM),")",sep="")
  ge1=data.frame("Trait Effects"=colnames(gamma),Estimates=ge1, pMCMC=ge1_p,check.names=FALSE)
  
  #****************************************************
  #Random effects ----
  #****************************************************
  if(Include_random == "yes") {
  #matrix for placing estimates into
  random_mod = matrix(nrow = model$samples*length(post_model$Omega[[1]]), ncol =length(post_model$Omega))
  
  for(i in 1:length(post_model$Omega)){
    #combine chains 
    tmp1 = do.call(rbind, post_model$Omega[[i]])
    tmp2= matrix(nrow = dim(tmp1)[1], ncol =model$ns)
    
    #extract variances for each species for each random effect
    for(j in 1:dim(tmp1)[1]){
    tmp3 = matrix(tmp1[j,],nrow = model$ns, ncol =model$ns)
    tmp3 = diag(tmp3)
    tmp2[j,] = tmp3 
    }
    
    #Average across species
    tmp2 = rowMeans(tmp2)
    random_mod[,i] = tmp2
  }
  
    random_mod<-as.mcmc(random_mod)
    
    #Name columns
    #if variances are not specified make them the same as those in the model
    if(is.null(randomvar_names)) {
      colnames(random_mod)<-names(model$ranLevels)
    } else  {colnames(random_mod)<-randomvar_names
    }
    
    #Summary of variance components
    randomVar = paste(round(posterior.mode(random_mod),dec_PM)," (",round(HPDinterval(random_mod)[,1],dec_PM), ", ",round(HPDinterval(random_mod)[,2],dec_PM),")",sep="")
    randomVar = data.frame("Random Effects"=colnames(random_mod),"Posterior Mode (CI)"=randomVar,"-"="",check.names=FALSE)
    
  } else  {
  }
  
  #****************************************************
  #% variation of explained by fixed and random effects averaged across species  ----
  #****************************************************
  if(VP_ave == "include") {
    VP = suppressWarnings(computeVariancePartitioning(model)$vals)
    VP = data.frame(rowMeans(VP))
    VP = VP %>% mutate(across(everything(), ~round(., 4)*100))
    
    if(is.null(VPnames)) {
    } else  {
      rownames(VP) = VPnames
    }
    
    VP = data.frame("Variance Partitioning"=rownames(VP),"%"=VP[,1],check.names=FALSE)
    rownames(VP)<-NULL
  } else  {
  }
  
  #****************************************************
  #Fit statistics ----
  #****************************************************
  # To assess model fit in terms of $R^2$, we apply the `evaluateModelFit` function to the posterior predictive distribution computed by the function `computePredictedValues`.
  model_preds = computePredictedValues(model)
  model_fit = as.data.frame(evaluateModelFit(model, predY=model_preds))
  
  #Round fit stats
  model_fit = model_fit %>% dplyr::summarise(across(everything(), ~mean(., na.rm=T))) %>% mutate(across(everything(), ~round(., 2))) 
  
  #****************************************************
  #Phylo effects ----
  #****************************************************
  #Rho = does not ask whether the species traits are correlated with the phylogeny, a question that is often in the focus of phylogenetic comparative analyses. Instead, the phylogenetic signal parameter ρ measures whether the species niches (i.e., their responses to the environmental covariates, as measured by the β parameters) show phylogenetic correlations.
  if(is.null(post_model$Rho)) {
  } else  {
    rho = as.mcmc(unlist(post_model$Rho))
    rho = paste(round(posterior.mode(rho),dec_PM)," (",round(HPDinterval(rho)[,1],dec_PM), ", ",round(HPDinterval(rho)[,2],dec_PM),")",sep="")
    rho = data.frame("Phylogenetic Effects"="Rho", "Posterior Mode (CI)"=rho,check.names=FALSE)
  }
  
  #****************************************************
  ##Excel output: fixed and random effects ----
  #****************************************************
  #Fixed effects
  fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
  fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
  #Fixed differences
  fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
  fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
  
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
  header=data.frame(col1=c(""),col2=c(""),col3=c(""),col4=c(""),col5=c(""),col6=c(""),col7=c(""),col8=c(""),col9=c(""))
  colnames(header)<-c(title,"","","","","","","","")
  writeData(workbook, sheet, header, startCol = 1, startRow = start_row,headerStyle = hs1)
  
  #Fixed effects
  writeData(workbook, sheet, fixedeff, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
  row_nums = start_row+dim(header)[1] + dim(fixedeff)[1]+1
  #Remove column headings if deleting fixed effects as it will be assessing higher order interactions where column names are not needed
  if(any(fixed_diffinc == "none")) { #Do not write fixeddiffinc dataframe if = "none"
    } else  {
    writeData(workbook, sheet, fixeddiff, startCol = 1, startRow = start_row+dim(header)[1] +dim(fixedeff)[1]+1,headerStyle = hs2)
    row_nums = row_nums + dim(fixeddiff)[1]+1
    }
  
  #Trait effects
  if(traits=="include") { #Write trait effects if = "include"
    writeData(workbook, sheet, ge1, startCol = 1, startRow = row_nums,headerStyle = hs2)
    row_nums = row_nums + dim(ge1)[1]+1
    } else  {
 }
  
  #Bold pMCMC values less than 0.05
  bolding<-createStyle(textDecoration="bold")
  conditionalFormatting(workbook, sheet, cols=3, rows=1:10000, rule="<0.05", style = bolding)
  
  #Random effects: variances
  #Should they be outputted or not
  if(Include_random == "yes") {
    writeData(workbook, sheet, randomVar, startCol = 1, startRow = row_nums,headerStyle = hs2)
    row_nums = row_nums + dim(randomVar)[1]+2
  } else  {
    randomVar = data.frame()
  }
  
  if(VP_ave == "include") {
    writeData(workbook, sheet, VP, startCol = 1, startRow = row_nums,headerStyle = hs2)
    row_nums = row_nums + dim(VP)[1]+3
  } else  {
  }
  
  if(is.null(post_model$Rho)) {
    workbook=workbook  
  } else  {
    writeData(workbook, sheet, rho, startCol = 1, startRow = row_nums,headerStyle = hs2)
    row_nums = row_nums + dim(rho)[1]+4
    }
  
  writeData(workbook, sheet, "Fit statistics", startCol = 1, startRow = start_row+dim(header)[1]+dim(fixedeff)[1]+1+dim(ge1)[1]+1+ifelse(dim(fixeddiff)[1] > 0,dim(fixeddiff)[1]+1,0)+dim(randomVar)[1]+2+dim(VP)[1]+5,headerStyle = hs2)
  writeData(workbook, sheet, model_fit, startCol = 1, startRow = start_row+dim(header)[1]+dim(fixedeff)[1]+1+dim(ge1)[1]+1+ifelse(dim(fixeddiff)[1] > 0,dim(fixeddiff)[1]+1,0)+dim(randomVar)[1]+2+dim(VP)[1]+7,headerStyle = hs2)
  
#****************************************************
#Section 2: Community comparisons ----
#****************************************************

  if(is.null(community_comparisons)) {
  #nothing to do done
  } else  {
  #Calculate community comparisons  
  community = data.frame()
  
  for(i in 1:length(community_comparisons)){
    composition_var = names(community_comparisons)[i]
    composition_metric = community_comparisons[[i]]$composition_metric
    comp_names = community_comparisons[[i]]$comp_names
    
    #setup predictions: need to setup differently for continuous and categorical variables
    if(is.null(community_comparisons[[i]]$composition_comp)) {
      #continuous effects
      ngrid = ifelse(is.null(community_comparisons[[i]]$ngrid),5,community_comparisons[[i]]$ngrid) #number of points along the gradient 
      
      #Construct gradient for predictions
      comp_1 = constructGradient(model, focalVariable = composition_var,      
                                 non.focalVariables = 1,
                                 ngrid = ngrid)
      
      #Predict values
      comp_1_predictions = predict(model, 
                                   XData=comp_1$XDataNew, 
                                   studyDesign=comp_1$studyDesignNew, 
                                   ranLevels=comp_1$rLNew)
      
      #rename
      colnames(comp_1$XDataNew)[1] = comp_names[1]
      composition_var = comp_names[1]
      
    } else  {
      #categorical effects
      composition_comp = community_comparisons[[i]]$composition_comp
      #Construct predictions
      comp_1 = constructGradient(model, focalVariable = composition_var,      
                                 non.focalVariables = 1)
      
      #Predict values
      comp_1_predictions = predict(model, 
                                   XData=comp_1$XDataNew, 
                                   studyDesign=comp_1$studyDesignNew, 
                                   ranLevels=comp_1$rLNew)
      #rename
      comp_1$XDataNew[,1] = comp_names
      comp_1$XDataNew[,1] = as.factor(comp_1$XDataNew[,1])
    }
    
    #Create a dataframe to write dissimilarity metric to
    comp_res = data.frame()
    comp_diff = data.frame()
    
    #Comparing communities to randomised community
    if(randomisation == "poisson") {
    
    #Compared to communities sample from Poisson distribution with mean = to prevalence of each taxa
    #for each iteration calculate dissimilarity and write to results to dataframe
    for(j in 1:length(comp_1_predictions)){
      dissim = as.data.frame(comp_1_predictions[j]) #model predictions
      
      #Code below creates random expectation by sampling from a Poisson distribution with mean = to prevalence of each taxa. If jaccard then set max to 1 (e.g. present)
      tax_means = as.data.frame(model$Y)
      tax_means = colMeans(tax_means) #mean prevalence/abundance of each taxa
      
      if (composition_metric == "jaccard") {
        #dissimR = matrix(rpois(length(dissim), mean(as.matrix(dissim))), nrow = nrow(dissim), ncol = ncol(dissim))
        dissimR <- matrix(nrow = nrow(dissim), ncol = ncol(dissim))
        for (k in 1:ncol(dissim)) {
          dissimR[, k] <- rpois(nrow(dissim), tax_means[k])
        }
        dissimR = pmin(dissimR, 1)
      } else {
        if (response_logged == TRUE) {
          #Backtransform if response is logged
          dissim = exp(as.matrix(dissim))
          dissimR <- matrix(nrow = nrow(dissim), ncol = ncol(dissim))
          for (k in 1:ncol(dissim)) {
            dissimR[, k] <- rpois(nrow(dissim), tax_means[k])
          }
          dissimR = log(dissimR+1)
        } else {
          dissimR <- matrix(nrow = nrow(dissim), ncol = ncol(dissim))
          for (k in 1:ncol(dissim)) {
            dissimR[, k] <- rpois(nrow(dissim), tax_means[k])
          }
        }
      }
      
      #rename
      dissim_names = data.frame(code =rownames(comp_1$XDataNew),names=comp_1$XDataNew[,composition_var])
      dissim_names = dissim_names %>% mutate_if(is.numeric, round, 1)
      
      #Calculate similarity
      real = melt(as.matrix(vegdist(dissim, method = composition_metric)), varnames = c("code1", "code2"), value.name = composition_metric) #dissimilarity measure on real data
      ran = melt(as.matrix(vegdist(dissimR, method = composition_metric)), varnames = c("code1", "code2"), value.name = paste(composition_metric,"R",sep="")) #dissimilarity measure on randomised data
      
      tmp = left_join(real,ran, by = c("code1" = "code1", "code2" = "code2")) 
      tmp = tmp %>% mutate(diff=tmp[,3]-tmp[,4], #calculate difference in dissimlarity metric between real and randomised data
                           name_1=dissim_names$names[match(code1,dissim_names$code)], #change actual names
                           name_2=dissim_names$names[match(code2,dissim_names$code)], #change actual names
                           name_12=paste0(name_1," vs ",name_2)) %>% 
        filter(name_1 != name_2) %>%
        dplyr::select(-c(code1,code2,name_1,name_2)) 

      tmp1 = tmp %>% dplyr::select(name_12, composition_metric) %>% pivot_wider(names_from = name_12, values_from = composition_metric) #matrix of real values
      tmp2 = tmp %>% dplyr::select(name_12, diff) %>% pivot_wider(names_from = name_12, values_from = diff) #matrix of differences between real and randomised values

      
      comp_res = rbind(comp_res,tmp1) #combine estimates from different iterations
      comp_diff = rbind(comp_diff,tmp2) #combine difference estimates from different iterations
      }
      
    }  else  {
      
    #*******************************************  
    #Randomise response data: takes along time
    rand_list <- vector("list", length(comp_1_predictions))
    for (j in 1:length(comp_1_predictions)) {
      rand_list[[j]] <- model$Y[sample(nrow(model$Y)), sample(ncol(model$Y))]
    }
    
    #for each iteration calculate dissimilarity and write to results to dataframe
    for(k in 1:length(comp_1_predictions)){
      #Predict real data for each iteration
      dissim = as.data.frame(comp_1_predictions[k]) #model predictions 
      
      #run model on randomised data for each model
      modelR = Hmsc(Y = rand_list[[k]], XData=model$XData,
                 XFormula = as.formula(model$XFormula),
                 studyDesign = model$studyDesign,
                 ranLevels=model$ranLevels, #this models residual correlations between species
                 phyloTree = model$phyloTree,
                 distr = model$distr)
      #Run
      modelR = sampleMcmc(modelR, thin = 1, samples = 2, transient = 0,
                       nChains = 1, verbose = 0, nParallel = 1,updater=list(GammaEta=FALSE))
      #Construct predictions
      comp_1R = constructGradient(modelR, focalVariable = composition_var,non.focalVariables = 1)
      
      #Predict values
      comp_1_predictionsR = predict(modelR, 
                                    XData=comp_1R$XDataNew, 
                                    studyDesign=comp_1R$studyDesignNew, 
                                    ranLevels=comp_1R$rLNew)
      dissimR = as.data.frame(comp_1_predictionsR[2]) 
     
      #rename
      dissim_names = data.frame(code =rownames(comp_1$XDataNew),names=comp_1$XDataNew[,composition_var])
      dissim_names = dissim_names %>% mutate_if(is.numeric, round, 1)
      
      #Calculate similarity
      real = melt(as.matrix(vegdist(dissim, method = composition_metric)), varnames = c("code1", "code2"), value.name = composition_metric) #dissimilarity measure on real data
      ran = melt(as.matrix(vegdist(dissimR, method = composition_metric)), varnames = c("code1", "code2"), value.name = paste(composition_metric,"R",sep="")) #dissimilarity measure on randomised data
      
      tmp = left_join(real,ran, by = c("code1" = "code1", "code2" = "code2")) 
      tmp = tmp %>% mutate(diff=tmp[,3]-tmp[,4], #calculate difference in dissimlarity metric between real and randomised data
                           name_1=dissim_names$names[match(code1,dissim_names$code)], #change actual names
                           name_2=dissim_names$names[match(code2,dissim_names$code)], #change actual names
                           name_12=paste0(name_1," vs ",name_2)) %>% 
        filter(name_1 != name_2) %>%
        dplyr::select(-c(code1,code2,name_1,name_2)) 
      
        tmp1 = tmp %>% dplyr::select(name_12, 1) %>% pivot_wider(names_from = name_12, values_from = composition_metric) #matrix of real values
        tmp2 = tmp %>% dplyr::select(name_12, diff) %>% pivot_wider(names_from = name_12, values_from = diff) #matrix of differences between real and randomised values

      comp_res = rbind(comp_res,tmp1) #combine estimates from different iterations
      comp_diff = rbind(comp_diff,tmp2) #combine difference estimates from different iterations
      
    }  
    }  
    
    #remove rows where differences between communities weren't possible to calculate. This can happen where categories to be compared both have 0 presence.
    comp_res = comp_res[complete.cases(comp_res),] 
    comp_diff = comp_diff[complete.cases(comp_diff),] 
    
    #select specified comparisons if specified if not keep all
    if(is.null(community_comparisons[[i]]$composition_comp)) {
      comp_res = as.mcmc(comp_res)
      comp_diff = as.mcmc(comp_diff) 
    } else  {
      comp_res = as.mcmc(comp_res[,composition_comp])
      comp_diff = as.mcmc(comp_diff[,composition_comp]) 
    } 
    
    #Pvalue of differences
    nC = dim(comp_diff)[2]
    cc1_p = pmax(0.5/dim(comp_diff)[1], pmin(colSums(comp_diff[,1:nC, drop = FALSE] > 0)/dim(comp_diff)[1], 1 - colSums(comp_diff[, 1:nC, drop = FALSE] > 0)/dim(comp_diff)[1]))*2
    
    #Format estimates
    cc1_res=paste(round(posterior.mode(comp_res),dec_PM)," (",round(HPDinterval(comp_res)[,1],dec_PM), ", ",round(HPDinterval(comp_res)[,2],dec_PM),")",sep="")
    cc1_diff=paste(round(posterior.mode(comp_diff),dec_PM)," (",round(HPDinterval(comp_diff)[,1],dec_PM), ", ",round(HPDinterval(comp_diff)[,2],dec_PM),")",sep="")
    
    cc1=data.frame("Variable"=composition_var,"Community Comparison"=colnames(comp_res),"Posterior Mode (CI)"=cc1_res, "Difference from Random Posterior Mode (CI)"=cc1_diff, "pMCMC"=round(as.numeric(cc1_p),3),check.names=FALSE)
    
    #Combine results
    community = rbind(community,cc1)
  }
  
  #****************************************************
  ##Excel output: community comparisons ----
  #****************************************************
  #Create new sheet 
  sheet2 = paste(sheet,"Community_comparisons",sep="_")
  start_row = 1
  addWorksheet(workbook, sheet2)
  
  header=data.frame(col1=c(""),col2=c(""),col3=c(""),col4=c(""),col5=c(""))
  colnames(header)<-c(paste(title,": Community comparisons",sep=" "),"","","","")
  
  #table title
  writeData(workbook, sheet2, header, startCol = 1, startRow = start_row,headerStyle = hs1)
  
  #Community comparisons
  writeData(workbook, sheet2, community, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
  
  #Bold pMCMC values less than 0.05
  bolding<-createStyle(textDecoration="bold")
  conditionalFormatting(workbook, sheet2, cols=5, rows=1:10000, rule="<0.05", style = bolding)
  
  }

  
  #****************************************************
  #Section 3: per species values ----
  #****************************************************
  #should species estimates be included then proceed, otherwise return workbook
  if(Include_species == "include") {
    
    #****************************************************
    #Fixed effects ----
    #****************************************************
    fixed_mod = do.call(rbind, post_model$Beta)
    fixed_mod = as.mcmc(fixed_mod)
    
    #Rename model fixed effects
    if(is.null(fixed_names)) {
      
      fixed_names = character()
      for(i in 1:length(model$spNames)){
        tmp = paste(model$covNames,model$spNames[i],sep=": ")
        fixed_names = c(fixed_names,tmp)
      }
      colnames(fixed_mod)<-fixed_names
    } else  {
      tmp = character()
      for(i in 1:length(model$spNames)){
        tmp2 = paste(fixed_names,model$spNames[i],sep=": ")
        tmp = c(tmp,tmp2)
      }
      colnames(fixed_mod)<-tmp
    }
    
    #P values using summary.MCMCglmm code
    nG=dim(fixed_mod)[2]
    #Pvalues = option to exclude 
    if(any(pvalues_species == "include")) {
      fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
      fe1_p=round(as.numeric(fe1_p),3)
    } else  {
      if(any(pvalues_species == "exclude")) {
        fe1_p=rep("-",length(fixed_names))
      } else  {
        fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
        fe1_p=round(as.numeric(fe1_p),3)
        fe1_p[pvalues_species]<-"-"
      }
    }
    
    fe1=paste(round(posterior.mode(fixed_mod),dec_PM)," (",round(HPDinterval(fixed_mod)[,1],dec_PM), ", ",round(HPDinterval(fixed_mod)[,2],dec_PM),")",sep="")
    fe1=data.frame(Fixed_Effects=colnames(fixed_mod),Estimates=fe1, pMCMC=fe1_p)
    
    #****************************************************
    ##Differences between fixed effects ----
    #****************************************************
    
    ##fixed_diffinc ----
    if(any(fixed_diffinc_species == "none")) {
      fixed=fe1
    } else  {
      #All differences between fixed effects: note this creates a matrix of all fixedeffects X all fixedeffects so will max out if there are many species.
      fe2<-pairwise.diffs(fixed_mod,nF=nF)
      fixed=rbind(fe1,fe2)
      fixed = fixed %>% dplyr::filter(Fixed_Effects %in% c(fixed_names,fixed_diffinc_species) == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    }
    
    #****************************************************
    #Random effects ----
    #****************************************************
    if(Include_random_species == "yes") {
    #matrix for placing estimates into
    random_mod = matrix(nrow = model$samples*length(post_model$Omega[[1]]))
    
    for(i in 1:length(post_model$Omega)){
      #combine chains 
      tmp1 = do.call(rbind, post_model$Omega[[i]])
      tmp2= matrix(nrow = dim(tmp1)[1], ncol =model$ns)
      
      #extract variances for each species for each random effect
      for(j in 1:dim(tmp1)[1]){
        tmp3 = matrix(tmp1[j,],nrow = model$ns, ncol =model$ns)
        tmp3 = diag(tmp3)
        tmp2[j,] = tmp3 
      }
      
      random_mod = cbind(random_mod,tmp2)
    }
    random_mod = random_mod[,-1]
    random_mod<-as.mcmc(random_mod)
    
    #Rename model random effects
    if(is.null(random_names_species)) {
      
      random_names_species = character()
      for(i in 1:length(model$rLNames)){
        tmp = paste(model$rLNames[i],model$spNames,sep=": ")
        random_names_species = c(random_names_species,tmp)
      }
      colnames(random_mod)<-random_names_species
    } else  {
      colnames(random_mod)<-random_names_species
    }
    
      #Summary of variance components
      randomVarspp = paste(round(posterior.mode(random_mod),dec_PM)," (",round(HPDinterval(random_mod)[,1],dec_PM), ", ",round(HPDinterval(random_mod)[,2],dec_PM),")",sep="")
      randomVarspp = data.frame("Random Effects"=colnames(random_mod),"Posterior Mode (CI)"=randomVarspp,"-"="",check.names=FALSE)
      
    } else  {
    }
    
    #****************************************************
    #% variation of explained by fixed and random effects for each species  ----
    #****************************************************
    if(VP_species == "include") {
      VP = as.data.frame(suppressWarnings(computeVariancePartitioning(model)$vals))
      VP = VP %>% mutate(across(everything(), ~round(., 4)*100))
      VP = data.frame("Variance Partitioning"=rownames(VP),"%"=VP,check.names=FALSE)
      colnames(VP) = sub("%.","% ",colnames(VP))
      rownames(VP)<-NULL
    } else  {
    }
    
    #****************************************************
    #Fit statistics ----
    #****************************************************
    # To assess model fit in terms of $R^2$, we apply the `evaluateModelFit` function to the posterior predictive distribution computed by the function `computePredictedValues`.
    model_preds = computePredictedValues(model)
    model_fit = as.data.frame(evaluateModelFit(model, predY=model_preds))
    
    #Round fit stats
    model_fit = model_fit %>% mutate(across(everything(), ~round(., 2)))
    model_fit = data.frame(Species=model$spNames,model_fit)
    
    #****************************************************
    ##Excel output: fixed and random effects ----
    #****************************************************
    #Fixed effects
    fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
    fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
    #Fixed differences
    fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
    fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
    
    #Create new sheet 
    sheet3 = paste(sheet,"species",sep="_")
    start_row = 1
    addWorksheet(workbook, sheet3)
    
    #
    header=data.frame(col1=c(""),col2=c(""),col3=c(""),col4=c(""),col5=c(""),col6=c(""),col7=c(""),col8=c(""),col9=c(""))
    colnames(header)<-c(paste(title,": species estimates",sep=" "),"","","","","","","","")
    
    #table title
    writeData(workbook, sheet3, header, startCol = 1, startRow = start_row,headerStyle = hs1)

    #Fixed effects
    writeData(workbook, sheet3, fixedeff, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
    row_nums = start_row+dim(header)[1]+dim(fixedeff)[1]+2
    
    #Remove column headings if deleting fixed effects as it will be assessing higher order interactions where column names are not needed
    if(any(fixed_diffinc_species == "none")) { #Do not write fixed_diffinc_species if "none"
      } else  {
      writeData(workbook, sheet3, fixeddiff, startCol = 1, startRow = start_row+dim(header)[1] +dim(fixedeff)[1]+1,headerStyle = hs2)
      row_nums = row_nums + dim(fixeddiff)[1] + 2
        }
    
    #Bold pMCMC values less than 0.05
    bolding<-createStyle(textDecoration="bold")
    conditionalFormatting(workbook, sheet3, cols=3, rows=1:10000, rule="<0.05", style = bolding)
    
    #Random effects: variances
    #Should they be outputted or not
    if(Include_random_species == "yes") {
      writeData(workbook, sheet3, randomVarspp, startCol = 1, startRow = row_nums,headerStyle = hs2)
      row_nums = row_nums + dim(randomVarspp)[1]+2
    } else  {
    }
    
    if(VP_species == "include") {
      writeData(workbook, sheet3, VP, startCol = 1, startRow = row_nums,headerStyle = hs2)
      row_nums = row_nums + dim(VP)[1]+2
      } else  {
    }
    
    writeData(workbook, sheet3, "Fit statistics", startCol = 1, startRow = row_nums,headerStyle = hs2)
    writeData(workbook, sheet3, model_fit, startCol = 1, startRow = row_nums+2,headerStyle = hs2)
    return(workbook)
    
    #If species estimate != "include" then just returns a workbook with estimates averaged across species
  } else  {
    return(workbook)
  }
}


#function for extracting df from xl workbook
xl_2_df = function(xltab,sheet=NULL){
  df<-readWorkbook(xltab,sheet=sheet,startRow = 2)
  colnames(df)<-gsub("[.]"," ",colnames(df))
  rownames(df)<-NULL
  return(df)
}

#function for df to Rmd table
hmsc_md = function(df){
  pacman::p_load(kableExtra)
  df = df %>% replace(is.na(.), "")
  kbl(df, align = "l", digits = 3) %>%
    kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) %>%
    row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Fixed Effect",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Trait Effects",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Random Effects",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Variance",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Phylogenetic Effects",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Fit statistics",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px")}
