#***************************************
#Function for processing Hmsc models#
#***************************************

HmscProc<-function(model=NULL,start_row=NULL,workbook=NULL, create_sheet="yes",sheet="sheet1",title="",fixed_names=NULL,fixed_diffinc="none",fixed_diff_diffs =NULL,fixed_diffinc_species="none",traits="exclude",pvalues = "include",VP_ave = "include",VPnames=NULL,randomvar_names=NULL,Include_random = "yes",Include_species ="exclude", VP_species = "include",random_names_species=NULL,Include_random_species = "yes",padding=4,dec_PM=2)
{ 
  #Explanation ----
  #1. Takes an Hmsc model and combines estimates from multiple chains and output 2 excel sheets: 1) averages across species; 2) Per species values. If there are multiple species these are averaged per mcmc sample using rowMeans to produce posterior distribution of average effects.
  #2. Estimates posterior modes and HPDintervals for effects
  #3. Calculates specified differences for fixed effects
  #4. For models with trait data, combines runs and estimates trait effects
  #5. Calculates average variance explained by random effects
  #6. Calculates % variation of explained by fixed and random effects 
  #7. For models with phylogenetic effects it outputs Rho
  #8. Calculates fit statistics
  #9. Output results to excel file
  
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
  # traits="include" #include trait effects: this examines the influence of species traits on community composition (e.g. richness if presence / absence of each species)
  # pvalues = "include" #include pvalues for comparisons or not
  # Include_random = "yes" #include random effect estimates or not
  # randomvar_names=c("R1","R2","R3") #names of random effects - will take from model object if not specified
  # VP_ave = "include" #include variance partitioning for species averages
  #VPnames = renaming of terms in Variance partitioning table
  # Include_species ="include" #should a separate sheet with species estimates be included?
  # fixed_diffinc_species="none" #differences between fixed effects to include for specific species e.g. fixed_diffinc_species=c(c("A: species1 vs B: species 1"))
  # VP_species = "include" #include variance partitioning for all species
  # Include_random_species = "yes" #include random effect variances for each species
  # random_names_species = NULL #names for random effects for each species
  # padding=4 #spacing in excel file
  # dec_PM=2) #decimal places of estimates
  
  #****************************************************
  #Section 1: averages per species ----
  #****************************************************
  
  #Load packages and naming
  pacman::p_load(Hmsc,coda,stringdist,runjags,matrixStats,openxlsx)
  
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
  tmp1 = combine.mcmc(post_model$Beta)
  
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
  if(pvalues == "include") {
    fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nF, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nF, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
    fe1_p=round(as.numeric(fe1_p),3)
  } else  {
    if(pvalues == "exclude") {
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
  gamma = combine.mcmc(post_model$Gamma)
  
  #name columns more clearly
  gamma_names = character()
  for(i in 1:length(model$trNames)){
    tmp = paste(model$trNames[i],model$covNames,sep="_")
    gamma_names = c(gamma_names,tmp)
  } 
  
  colnames(gamma) = gamma_names
  gamma = as.mcmc(gamma)
  
  #P values using summary.MCMCglmm code
  nG=dim(gamma)[2]
  #Pvalues = option to exclude 
  if(pvalues == "include") {
    ge1_p=pmax(0.5/dim(gamma)[1], pmin(colSums(gamma[,1:nG, drop = FALSE] > 0)/dim(gamma)[1], 1 - colSums(gamma[, 1:nG, drop = FALSE] > 0)/dim(gamma)[1]))*2
    ge1_p=round(as.numeric(ge1_p),3)
  } else  {
    if(pvalues == "exclude") {
      ge1_p=rep("-",length(gamma_names))
    } else  {
      ge1_p=pmax(0.5/dim(gamma)[1], pmin(colSums(gamma[,1:nG, drop = FALSE] > 0)/dim(gamma)[1], 1 - colSums(gamma[, 1:nG, drop = FALSE] > 0)/dim(gamma)[1]))*2
      ge1_p=round(as.numeric(ge1_p),3)
      ge1_p[pvalues]<-"-"
    }
  }
  
  ge1=paste(round(posterior.mode(gamma),dec_PM)," (",round(HPDinterval(gamma)[,1],dec_PM), ", ",round(HPDinterval(gamma)[,2],dec_PM),")",sep="")
  ge1=data.frame("Trait Effects"=colnames(gamma),Estimates=ge1, pMCMC=ge1_p,check.names=FALSE)
  
  #****************************************************
  #Random effects ----
  #****************************************************
  #matrix for placing estimates into
  random_mod = matrix(nrow = model$samples*length(model$repList), ncol =length(post_model$Omega))
  
  for(i in 1:length(post_model$Omega)){
    #combine chains 
    tmp1 = combine.mcmc(post_model$Omega[i])
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
    
    if(Include_random == "yes") {
    #Summary of variance components
    rand1 = paste(round(posterior.mode(random_mod),dec_PM)," (",round(HPDinterval(random_mod)[,1],dec_PM), ", ",round(HPDinterval(random_mod)[,2],dec_PM),")",sep="")
    rand1 = data.frame("Random Effects"=colnames(random_mod),"Posterior Mode (CI)"=rand1,"-"="",check.names=FALSE)
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
  #Excel output: fixed and random effects ----
  #****************************************************
  #Fixed effects
  fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
  fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
  #Fixed differences
  fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
  fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
  #Random
  randomVar<-rand1
  
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
  #Section 2: per species values ----
  #****************************************************
  #should species estimates be included then proceed, otherwise return workbook
  if(Include_species == "include") {
    
    #****************************************************
    #Fixed effects ----
    #****************************************************
    fixed_mod = combine.mcmc(post_model$Beta)
    
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
    if(pvalues == "include") {
      fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
      fe1_p=round(as.numeric(fe1_p),3)
    } else  {
      if(pvalues == "exclude") {
        fe1_p=rep("-",length(fixed_names))
      } else  {
        fe1_p=pmax(0.5/dim(fixed_mod)[1], pmin(colSums(fixed_mod[,1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1], 1 - colSums(fixed_mod[, 1:nG, drop = FALSE] > 0)/dim(fixed_mod)[1]))*2
        fe1_p=round(as.numeric(fe1_p),3)
        fe1_p[pvalues]<-"-"
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
    #matrix for placing estimates into
    random_mod = matrix(nrow = model$samples*length(model$repList))
    
    for(i in 1:length(post_model$Omega)){
      #combine chains 
      tmp1 = combine.mcmc(post_model$Omega[1])
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
    
    
    if(Include_random_species == "yes") {
      #Summary of variance components
      rand1 = paste(round(posterior.mode(random_mod),dec_PM)," (",round(HPDinterval(random_mod)[,1],dec_PM), ", ",round(HPDinterval(random_mod)[,2],dec_PM),")",sep="")
      rand1 = data.frame("Random Effects"=colnames(random_mod),"Posterior Mode (CI)"=rand1,"-"="",check.names=FALSE)
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
    #Excel output: fixed and random effects ----
    #****************************************************
    #Fixed effects
    fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
    fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
    #Fixed differences
    fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
    fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
    
    #Random
    randomVar<-rand1
    
    #Create new sheet 
    sheet2 = paste(sheet,"species",sep="_")
    start_row = 1
    addWorksheet(workbook, sheet2)
    
    #
    header=data.frame(col1=c(""),col2=c(""),col3=c(""),col4=c(""),col5=c(""),col6=c(""),col7=c(""),col8=c(""),col9=c(""))
    colnames(header)<-c(paste(title,"species estimates",sep=" "),"","","","","","","","")
    
    #table title
    writeData(workbook, sheet2, header, startCol = 1, startRow = start_row,headerStyle = hs1)

    #Fixed effects
    writeData(workbook, sheet2, fixedeff, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
    row_nums = start_row+dim(header)[1]+dim(fixedeff)[1]+2
    
    #Remove column headings if deleting fixed effects as it will be assessing higher order interactions where column names are not needed
    if(any(fixed_diffinc_species == "none")) { #Do not write fixed_diffinc_species if "none"
      } else  {
      writeData(workbook, sheet2, fixeddiff, startCol = 1, startRow = start_row+dim(header)[1] +dim(fixedeff)[1]+1,headerStyle = hs2)
      row_nums = row_nums + dim(fixeddiff)[1] + 2
        }
    
    #Bold pMCMC values less than 0.05
    bolding<-createStyle(textDecoration="bold")
    conditionalFormatting(workbook, sheet2, cols=3, rows=1:10000, rule="<0.05", style = bolding)
    
    #Random effects: variances
    #Should they be outputted or not
    if(Include_random_species == "yes") {
      writeData(workbook, sheet2, randomVar, startCol = 1, startRow = row_nums,headerStyle = hs2)
      row_nums = row_nums + dim(randomVar)[1]+2
    } else  {
    }
    
    if(VP_species == "include") {
      writeData(workbook, sheet2, VP, startCol = 1, startRow = row_nums,headerStyle = hs2)
      row_nums = row_nums + dim(VP)[1]+2
      } else  {
    }
    
    writeData(workbook, sheet2, "Fit statistics", startCol = 1, startRow = row_nums,headerStyle = hs2)
    writeData(workbook, sheet2, model_fit, startCol = 1, startRow = row_nums+2,headerStyle = hs2)
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
