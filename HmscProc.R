#***************************************
#Function for processing Hmsc models#
#***************************************

#Function ----

HmscProc<-function(model=NULL,start_row=NULL,workbook=NULL, create_sheet="yes",sheet="sheet1",title="",fixed_names=NULL,fixed_diffinc="all",fixed_diff_diffs =NULL,pvalues = "include",VP_species = "include",randomvar_names=NULL,Include_random = "yes",padding=4,dec_PM=2)
{ 
  #Explanation ----
  #1. Takes an Hmsc model and combines estimates from multiple chains
  #2. Estimates posterior modes and HPDintervals for effects
  #3. Calculates specified differences for fixed effects
  #4. Calculates % variation explained by random effects [not yet implemented]
  #5. Calculates % variation explained per species
  #6. Calculates fit statistics
  #7. Output results to data.frame
  
  #Troubleshooting aid
  # model=model
  # workbook=NULL
  # start_row=1
  # create_sheet="yes"
  # sheet="Hm_1"
  # title="Table for X"
  # fixed_names=c()
  # fixed_diffinc=c()
  # pvalues = "include"
  # VP_species = "include"
  # randomvar_names=c()
  # fixed_diff_diffs =NULL
  # Include_random = "yes"
  # padding=4
  # dec_PM=2
  
  #Load packages and naming ----
  pacman::p_load(Hmsc,coda,stringdist,runjags,openxlsx)
  
  #Convert model object
  model=model
  post_model = convertToCodaObject(model)
  
  #****************************************************
  #Fixed effects ----
  #****************************************************
  fixed_mod = combine.mcmc(post_model$Beta)
  
  #Rename model fixed effects
  if(is.null(fixed_names)) {
    colnames(fixed_mod)<- colnames(fixed_mod)
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
  fixed=rbind(fe1,fe2)
  
  ##fixed_diffinc ----
  if(any(fixed_diffinc == "all")) {
  } else  {
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
  #Random effects ----
  #****************************************************
  if(length(post_model$Omega)>1) {
    
    #Create matrix to put estimates in
    random_mod = matrix(NA,length(post_model$Omega[[1]][[1]])*length(model$repList),length(post_model$Omega))
    
    #extract estimates
    for(i in 1:length(post_model$Omega)){
      random_mod[,i]<-combine.mcmc(post_model$Omega[i])
      random_mod<-as.mcmc(random_mod)
    }
    
    #Name columns
    #if variances are not specified make them the same as those in the model
    if(is.null(randomvar_names)) {
      colnames(random_mod)<-names(model$ranLevels)
    } else  {colnames(random_mod)<-randomvar_names
    }
    
    #Summary of variance components
    rand1 = paste(round(posterior.mode(random_mod),dec_PM)," (",round(HPDinterval(random_mod)[,1],dec_PM), ", ",round(HPDinterval(random_mod)[,2],dec_PM),")",sep="")
    rand1 = data.frame("Random Effects"=colnames(random_mod),"Posterior Mode (CI)"=rand1,"-"="",check.names=FALSE)
  } else  {
  }
  
  #% variation per species ----
  if(VP_species == "include") {
    VP = data.frame(computeVariancePartitioning(model)$vals)
    VP = VP %>% mutate(across(everything(), ~round(., 2)))
    VP = data.frame("Variance Partitioning"=rownames(VP),VP,check.names=FALSE)
    rownames(VP)<-NULL
  } else  {
  }
  
  #Fit statistics ----
  # To assess model fit in terms of $R^2$, we apply the `evaluateModelFit` function to the posterior predictive distribution computed by the function `computePredictedValues`.
  model_preds = computePredictedValues(model)
  model_fit = as.data.frame(evaluateModelFit(model, predY=model_preds))
  
  #Round fit stats
  model_fit = model_fit %>% mutate(across(everything(), ~round(., 2)))
  
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
  writeData(workbook, sheet, fixedeff, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
  
  #Remove column headings if deleting fixed effects as it will be assessing higher order interactions where column names are not needed
  if(any(fixed_diffinc == "none")) { #Do not write fixeddiff dataframe if fixed_diffinc == "none"
  } else  {
    writeData(workbook, sheet, fixeddiff, startCol = 1, startRow = start_row+dim(header)[1] +dim(fixedeff)[1]+1,headerStyle = hs2)
  }
  
  
  #Bold pMCMC values less than 0.05
  bolding<-createStyle(textDecoration="bold")
  conditionalFormatting(workbook, sheet, cols=3, rows=1:10000, rule="<0.05", style = bolding)
  
  #Random effects: variances
  #Should they be outputted or not
  if(Include_random == "yes") {
    writeData(workbook, sheet, randomVar, startCol = 1, startRow = start_row+dim(header)[1]+dim(fixedeff)[1]+1+ifelse(dim(fixeddiff)[1] > 0,dim(fixeddiff)[1]+1,0),headerStyle = hs2)
  } else  {
    workbook=workbook
  }
  
  if(VP_species == "include") {
    writeData(workbook, sheet, VP, startCol = 1, startRow = start_row+dim(header)[1]+dim(fixedeff)[1]+1+ifelse(dim(fixeddiff)[1] > 0,dim(fixeddiff)[1]+1,0)+dim(randomVar)[1]+2,headerStyle = hs2)
  } else  {
    workbook=workbook
  }
  
  writeData(workbook, sheet, "Fit statistics", startCol = 1, startRow = start_row+dim(header)[1]+dim(fixedeff)[1]+1+ifelse(dim(fixeddiff)[1] > 0,dim(fixeddiff)[1]+1,0)+dim(randomVar)[1]+2+dim(VP)[1]+2,headerStyle = hs2)
  writeData(workbook, sheet, model_fit, startCol = 1, startRow = start_row+dim(header)[1]+dim(fixedeff)[1]+1+ifelse(dim(fixeddiff)[1] > 0,dim(fixeddiff)[1]+1,0)+dim(randomVar)[1]+2+dim(VP)[1]+3,headerStyle = hs2)
  
  return(workbook)
}


#function for extracting df from xl workbook
xl_2_df = function(xltab,sheet=NULL){
  df<-readWorkbook(xltab,sheet=sheet)
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
    row_spec(grep("^Random Effects",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Variance",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(grep("^Fit statistics",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px")}
