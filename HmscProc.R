#***************************************
#Function for processing Hmsc models#
#***************************************

HmscProc<-function(model=NULL,start_row=NULL,workbook=NULL, create_sheet="yes",sheet="sheet1",title="",fixed_names=NULL,fixed_diffinc="none",fixed_diff_diffs =NULL,fixed_diffinc_species="none",pvalues = "include",traits="exclude",trait_names=NULL,pvalues_traits= "exclude",VP_ave = "include",VPnames=NULL,randomvar_names=NULL,Include_random = "yes",Include_species ="exclude",pvalues_species="exclude", VP_species = "include",random_names_species=NULL,Include_random_species = "yes",community_comparisons = NULL,community_null = "constant",response_logged = FALSE,padding=4,dec_PM=2)
{ 
  #Explanation ----
  #1. Takes an Hmsc model and combines estimates from multiple chains and outputs up to 4 excel sheets:
  #   averages across species; community comparisons; per-species estimates; per-species variance partitioning.
  #   Community averages use rowMeans across species per MCMC sample to produce a posterior of average effects.
  #2. Estimates posterior modes and HPD intervals for all effects.
  #3. Calculates specified pairwise differences for fixed effects.
  #4. If community_comparisons is specified, calculates pairwise community dissimilarity along focal
  #   covariate gradients and tests whether dissimilarity exceeds a model-based null (see community_null).
  #5. For models with trait data, combines chains and estimates trait effects (Gamma).
  #6. Calculates average variance explained by random effects (Omega diagonals).
  #7. Calculates % variation explained by fixed and random effects via computeVariancePartitioning.
  #8. For models with phylogenetic effects, outputs the phylogenetic signal parameter Rho.
  #9. Calculates model fit statistics (R^2) via evaluateModelFit on posterior predictive distribution.
  #10. Writes all results to an Excel workbook.
  # community_null = "constant" # null model used to test whether community dissimilarity exceeds chance.
  #   "constant": the focal covariate is fixed at its sample mean (continuous variables) or first
  #               factor level (categorical variables) for all prediction points. The null community
  #               represents expected composition when the focal covariate has no gradient — all
  #               sites experience the average value. Remaining dissimilarity under the null arises
  #               from stochastic sampling variation only. Use this when covariates are on any scale.
  #   "zero":     the focal covariate is fixed at zero for all prediction points. For mean-centred
  #               or standardised covariates (mean ≈ 0), this is equivalent to "constant". For raw
  #               uncentred covariates, zero may lie outside the observed data range — a warning is
  #               issued. Directly tests H0: the focal covariate's contribution to the linear
  #               predictor (beta × X = 0) generates no community turnover. Use this when the
  #               model covariates are mean-centred or standardised.
  # response_logged = FALSE # set TRUE if you log-transformed Y before fitting (e.g. log(Y+1) as
  #               a continuous response); back-transforms predictions before dissimilarity is computed.
  #               Not required for Hmsc Poisson family — predictions from expected=FALSE are already
  #               on the count scale.
  # community_comparisons = NULL # can take a list of lists e.g. "community_comparisons = list(Factor1 = list(composition_metric = "jaccard",comp_names = c("A","B","C","D"),composition_comp = c("A vs B", "B vs C")), Continuous1 = list(composition_metric = "jaccard",comp_names = "A",ngrid = 5))". Any "composition_metric" in vegdist function of vegan package (e.g."bray" or "jaccard") is allowed.
  # composition_var #Factor1, Continuous1 above = name of variable to compare composition metrics as it appears in the model formula 
  # comp_names #the names of the variable as you want to appear in the table. For factors in must include all levels even if a level is not included in the comparisons (e.g. "D" above).

  #****************************************************
  #Section 1: averages per species ----
  #****************************************************
  
  #Load packages
  pacman::p_load(Hmsc,coda,stringdist,matrixStats,openxlsx,vegan,reshape2,
                 dplyr,tidyr,purrr)
  
  #Convert model object
  post_model = convertToCodaObject(model)
  
  #****************************************************
  #Check terms are specified correctly ----
  #****************************************************
  if(is.null(fixed_names) & any(fixed_diffinc != "none")) {
    stop("fixed names needs to be specified for fixed_diffinc to be calculated")
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
  
  # Helper: compute pMCMC for columns of an MCMC sample matrix.
  # pvalues = "include"  → return values for all columns
  # pvalues = "exclude"  → return "-" for all columns
  # pvalues = integer index → return "-" for indexed columns, values for the rest
  calc_pmcmc <- function(samples_mat, pvalues, n_cols) {
    n_its <- nrow(samples_mat)
    p <- pmax(0.5 / n_its,
              pmin(colSums(samples_mat[, 1:n_cols, drop = FALSE] > 0) / n_its,
                   1 - colSums(samples_mat[, 1:n_cols, drop = FALSE] > 0) / n_its)) * 2
    p <- round(as.numeric(p), 3)
    if(any(pvalues == "exclude")) return(rep("-", n_cols))
    if(!any(pvalues == "include")) p[pvalues] <- "-"
    p
  }

  nF    <- dim(fixed_mod)[2]
  fe1_p <- calc_pmcmc(fixed_mod, pvalues, nF)
  
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
  if(!is.null(fixed_diff_diffs)) {
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
  if(!is.null(trait_names)) model$trNames <- trait_names

  # Build trait × covariate combination names for gamma columns.
  # Uses colnames(fixed_mod) rather than fixed_names to handle the case
  # where fixed_names was NULL (model$covNames was used as the column names).
  gamma_names <- character()
  for(i in seq_along(model$trNames)){
    gamma_names <- c(gamma_names, paste(model$trNames[i], colnames(fixed_mod), sep = "_"))
  }
  
  colnames(gamma) = gamma_names
  gamma = as.mcmc(gamma)
  
  nG    <- dim(gamma)[2]
  ge1_p <- calc_pmcmc(gamma, pvalues_traits, nG)
  
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
    
    # Extract the diagonal of each MCMC sample's species covariance matrix.
    # Each row of tmp1 is the flattened ns×ns Omega matrix for one sample;
    # the diagonal gives the variance (not covariance) for each species.
    tmp2 <- t(apply(tmp1, 1, function(row) diag(matrix(row, nrow = model$ns, ncol = model$ns))))

    #Average variances across species for each sample
    random_mod[, i] <- rowMeans(tmp2)
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
    
  }

  # Compute predictions and variance partitioning once here; both are reused in
  # Section 3 (per-species) to avoid running the same expensive call twice.
  model_preds <- computePredictedValues(model)
  VP_all      <- suppressWarnings(computeVariancePartitioning(model)$vals)

  #****************************************************
  #% variation explained by fixed and random effects averaged across species  ----
  #****************************************************
  if(VP_ave == "include") {
    VP <- data.frame(rowMeans(VP_all))
    VP <- VP |> mutate(across(everything(), ~round(., 4) * 100))
    if(!is.null(VPnames)) rownames(VP) <- VPnames
    VP <- data.frame("Variance Partitioning" = rownames(VP), "%" = VP[, 1], check.names = FALSE)
    rownames(VP) <- NULL
  }

  #****************************************************
  #Phylo effects ----
  #****************************************************
  # Rho measures whether species niches (β parameters) show phylogenetic signal —
  # i.e. whether closely related species respond similarly to environmental covariates.
  if(!is.null(post_model$Rho)) {
    rho <- as.mcmc(unlist(post_model$Rho))
    rho <- paste(round(posterior.mode(rho), dec_PM), " (",
                 round(HPDinterval(rho)[, 1], dec_PM), ", ",
                 round(HPDinterval(rho)[, 2], dec_PM), ")", sep = "")
    rho <- data.frame("Phylogenetic Effects" = "Rho", "Posterior Mode (CI)" = rho, check.names = FALSE)
  }

  #****************************************************
  #Fit statistics ----
  #****************************************************
  model_fit <- as.data.frame(evaluateModelFit(model, predY = model_preds))
  
  #Round fit stats
  model_fit = model_fit %>% dplyr::summarise(across(everything(), ~mean(., na.rm=T))) %>% mutate(across(everything(), ~round(., 2))) 
  
  #****************************************************
  ##Excel output: fixed and random effects ----
  #****************************************************
  #Fixed effects
  fixedeff <- fixed[!grepl(" vs ",fixed$Fixed_Effects),]
  fixedeff<-data.frame("Fixed Effects"=fixedeff$Fixed_Effects,"Posterior Mode (CI)"=fixedeff$Estimates,"pMCMC"=fixedeff$pMCMC,check.names=FALSE)
  #Fixed differences
  fixeddiff <- fixed[grepl(" vs ",fixed$Fixed_Effects),]
  fixeddiff<-data.frame("Fixed Effect Comparisons"=fixeddiff$Fixed_Effects,"Posterior Mode (CI)"=fixeddiff$Estimates,"pMCMC"=round(as.numeric(fixeddiff$pMCMC),3),check.names=FALSE)
  
  if(is.null(workbook))    workbook  <- createWorkbook()
  if(create_sheet == "yes") addWorksheet(workbook, sheet)
  if(is.null(start_row))   start_row <- 1
  
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
    row_nums = row_nums + dim(randomVar)[1]+1
  } else  {
    randomVar = data.frame()
  }
  
  if(VP_ave == "include") {
    writeData(workbook, sheet, VP, startCol = 1, startRow = row_nums, headerStyle = hs2)
    row_nums <- row_nums + dim(VP)[1] + 1
  }

  if(!is.null(post_model$Rho)) {
    writeData(workbook, sheet, rho, startCol = 1, startRow = row_nums, headerStyle = hs2)
    row_nums <- row_nums + dim(rho)[1] + 2
  }
  
  writeData(workbook, sheet, "Fit Statistics",startCol = 1, startRow = row_nums,headerStyle = hs2)
  addStyle(workbook, sheet, style = hs2, rows = row_nums, cols = 1:3, gridExpand = TRUE)

  writeData(workbook, sheet, model_fit, startCol = 1, startRow = row_nums+1,headerStyle = hs2)
  
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
    
    # Preserve the original model formula variable name before any display relabelling.
    # The original name is needed to look up values in model$XData and to pass the
    # correct column name to predict(), which uses the model formula internally.
    composition_var_orig <- composition_var

    # -------------------------------------------------------
    # Build the prediction gradient (constructGradient fixes all
    # non-focal variables at their means; the focal variable spans
    # its observed range for continuous variables or all levels for
    # categorical variables). Predict calls are deferred until after
    # the null design is built, so both use identical settings.
    # -------------------------------------------------------
    if(is.null(community_comparisons[[i]]$composition_comp)) {
      # Continuous focal variable: evenly spaced grid across observed range
      ngrid <- ifelse(is.null(community_comparisons[[i]]$ngrid), 5,
                      community_comparisons[[i]]$ngrid)
      comp_1 <- constructGradient(model,
                                  focalVariable      = composition_var_orig,
                                  non.focalVariables = 1,
                                  ngrid              = ngrid)
    } else {
      # Categorical focal variable: one row per factor level
      composition_comp <- community_comparisons[[i]]$composition_comp
      comp_1 <- constructGradient(model,
                                  focalVariable      = composition_var_orig,
                                  non.focalVariables = 1)
    }

    # -------------------------------------------------------
    # Build the NULL prediction design
    # -------------------------------------------------------
    # The null represents the community when the focal covariate is
    # held at a single constant value across all gradient points.
    # Because the only thing varying across gradient points is the
    # focal covariate (all others are fixed at their means by
    # constructGradient), fixing the focal covariate at a constant
    # value makes all prediction rows identical. The null community
    # therefore differs between gradient points only due to stochastic
    # sampling variation (from expected = FALSE, see below).
    #
    # community_null = "constant":
    #   Continuous: focal covariate fixed at its sample mean across
    #   all training sites. Null asks: "How dissimilar would
    #   communities be if every site experienced the average value
    #   of the focal covariate?"
    #   Categorical: all rows set to the first (reference) factor level.
    #
    # community_null = "zero":
    #   Continuous: focal covariate fixed at 0. For mean-centred /
    #   standardised covariates this equals the mean (identical to
    #   "constant"). For raw uncentred covariates, 0 may be outside
    #   the observed range — a warning is issued.
    #   Categorical: same as "constant" (reference level = X = 0
    #   in dummy coding).
    # -------------------------------------------------------
    XDataNull <- comp_1$XDataNew   # copy; we modify only the focal column

    if(community_null == "constant") {

      if(is.factor(comp_1$XDataNew[, composition_var_orig])) {
        # Categorical: fix all rows at the first factor level
        null_level <- levels(comp_1$XDataNew[, composition_var_orig])[1]
        XDataNull[, composition_var_orig] <- factor(
          null_level, levels = levels(comp_1$XDataNew[, composition_var_orig]))
      } else {
        # Continuous: fix at sample mean of the training data
        XDataNull[, composition_var_orig] <- mean(model$XData[, composition_var_orig])
      }

    } else if(community_null == "zero") {

      if(is.factor(comp_1$XDataNew[, composition_var_orig])) {
        # Categorical: reference level corresponds to all-zero dummy coding
        null_level <- levels(comp_1$XDataNew[, composition_var_orig])[1]
        XDataNull[, composition_var_orig] <- factor(
          null_level, levels = levels(comp_1$XDataNew[, composition_var_orig]))
      } else {
        # Continuous: fix at 0; warn if covariate is not mean-centred,
        # as X = 0 may then lie outside the observed data range
        covar_mean <- mean(model$XData[, composition_var_orig])
        covar_sd   <- sd(  model$XData[, composition_var_orig])
        if(abs(covar_mean) > 0.1 * covar_sd) {
          warning(paste0(
            "community_null = 'zero': covariate '", composition_var_orig,
            "' does not appear to be mean-centred ",
            "(mean = ", round(covar_mean, 3), ", SD = ", round(covar_sd, 3), "). ",
            "Fixing X = 0 may extrapolate outside the observed range. ",
            "Consider mean-centring the covariate or using community_null = 'constant'."
          ))
        }
        XDataNull[, composition_var_orig] <- 0
      }

    } else {
      stop("community_null must be 'constant' or 'zero'.")
    }

    # -------------------------------------------------------
    # Predict real and null communities from the posterior
    # -------------------------------------------------------
    # expected = FALSE samples actual observations from the posterior
    # predictive distribution — 0/1 for presence-absence (probit),
    # non-negative integers for counts (Poisson), continuous values
    # for normal responses. Using expected = FALSE for BOTH predictions
    # ensures they are on the same observed-data scale, making their
    # dissimilarity values directly comparable.
    #
    # When expected = TRUE (the default), all null prediction rows are
    # identical (same covariate value) so pairwise dissimilarity is
    # identically zero — the comparison would be trivial. expected = FALSE
    # introduces stochastic sampling variation, giving a meaningful
    # non-zero null dissimilarity that reflects chance fluctuations.
    # -------------------------------------------------------

    # Real community predictions: vary because the focal covariate varies
    comp_1_predictions <- predict(model,
                                  XData       = comp_1$XDataNew,
                                  studyDesign = comp_1$studyDesignNew,
                                  ranLevels   = comp_1$rLNew,
                                  expected    = FALSE)

    # Null community predictions: focal covariate held constant
    null_predictions <- predict(model,
                                XData       = XDataNull,
                                studyDesign = comp_1$studyDesignNew,
                                ranLevels   = comp_1$rLNew,
                                expected    = FALSE)

    # -------------------------------------------------------
    # Display relabelling (AFTER predict calls, for output only)
    # Renaming is purely cosmetic — it does not affect the predictions.
    # -------------------------------------------------------
    if(is.null(community_comparisons[[i]]$composition_comp)) {
      # Continuous: give the gradient column a user-friendly axis label
      colnames(comp_1$XDataNew)[1] <- comp_names[1]
      composition_var <- comp_names[1]
    } else {
      # Categorical: replace internal factor codes with user-supplied names
      comp_1$XDataNew[, 1] <- factor(comp_names, levels = comp_names)
    }

    # -------------------------------------------------------
    # For each posterior sample, compute pairwise dissimilarity
    # between real communities and between null communities, then
    # take the difference (real - null).
    # -------------------------------------------------------
    n_iter         <- length(comp_1_predictions)
    comp_res_list  <- vector("list", n_iter)   # raw dissimilarity, real communities
    comp_diff_list <- vector("list", n_iter)   # real dissimilarity minus null dissimilarity

    for(j in seq_len(n_iter)){

      # Extract predicted community matrix for this posterior sample.
      # Rows = gradient levels (sites); columns = species.
      dissim  <- as.data.frame(comp_1_predictions[[j]])   # real community
      dissimR <- as.data.frame(null_predictions[[j]])     # null community

      # Back-transform if the response was modelled on a log scale.
      # Only needed when the user manually log-transformed Y before
      # fitting as a continuous response; not needed for Hmsc Poisson.
      if(response_logged == TRUE) {
        dissim  <- exp(dissim)  - 1
        dissimR <- exp(dissimR) - 1
      }

      # Mapping from internal gradient row codes to display-level names
      dissim_names <- data.frame(
        code  = rownames(comp_1$XDataNew),
        names = comp_1$XDataNew[, composition_var]
      ) |>
        mutate(across(where(is.numeric), ~round(., 1)))

      # Pairwise dissimilarity between real predicted communities
      real <- melt(
        as.matrix(vegdist(dissim,  method = composition_metric)),
        varnames   = c("code1", "code2"),
        value.name = composition_metric
      )

      # Pairwise dissimilarity between null communities.
      # Under the null, all rows of dissimR are independently sampled
      # from the same community (focal covariate fixed), so this captures
      # stochastic variation only — the baseline we compare against.
      ran <- melt(
        as.matrix(vegdist(dissimR, method = composition_metric)),
        varnames   = c("code1", "code2"),
        value.name = paste0(composition_metric, "_null")
      )

      # Difference = real dissimilarity - null dissimilarity.
      # Positive values mean communities are MORE different than the
      # null predicts: i.e., the focal covariate drives extra turnover.
      tmp <- left_join(real, ran, by = c("code1", "code2")) |>
        mutate(
          diff    = .data[[composition_metric]] - .data[[paste0(composition_metric, "_null")]],
          name_1  = dissim_names$names[match(code1, dissim_names$code)],
          name_2  = dissim_names$names[match(code2, dissim_names$code)],
          name_12 = paste0(name_1, " vs ", name_2)
        ) |>
        filter(name_1 != name_2) |>
        dplyr::select(-c(code1, code2, name_1, name_2))

      # Reshape to wide format: one column per community pair
      comp_res_list[[j]]  <- tmp |>
        dplyr::select(name_12, all_of(composition_metric)) |>
        pivot_wider(names_from = name_12, values_from = all_of(composition_metric))

      comp_diff_list[[j]] <- tmp |>
        dplyr::select(name_12, diff) |>
        pivot_wider(names_from = name_12, values_from = diff)
    }

    # Combine across posterior samples (rows = samples, columns = community pairs)
    comp_res  <- do.call(rbind, comp_res_list)
    comp_diff <- do.call(rbind, comp_diff_list)

    # Remove samples where dissimilarity could not be computed
    # (occurs when both communities have zero presence for all species)
    comp_res  <- comp_res[ complete.cases(comp_res),  ]
    comp_diff <- comp_diff[complete.cases(comp_diff), ]

    # Select only the specified pairwise comparisons; keep all if not specified
    if(is.null(community_comparisons[[i]]$composition_comp)) {
      comp_res  <- as.mcmc(comp_res)
      comp_diff <- as.mcmc(comp_diff)
    } else {
      comp_res  <- as.mcmc(comp_res[,  composition_comp, drop = FALSE])
      comp_diff <- as.mcmc(comp_diff[, composition_comp, drop = FALSE])
    }

    # -------------------------------------------------------
    # pMCMC for the difference from null
    # -------------------------------------------------------
    # Tests whether communities are credibly MORE dissimilar than the
    # null model predicts. pMCMC < 0.05 means the focal covariate drives
    # community turnover beyond what stochastic variation alone generates.
    # The test is two-tailed: negative values (communities less dissimilar
    # than null) would also give small pMCMC but indicate the opposite.
    nC    <- ncol(as.matrix(comp_diff))
    n_its <- nrow(as.matrix(comp_diff))
    cc1_p <- pmax(0.5 / n_its,
                  pmin(colSums(as.matrix(comp_diff)[, 1:nC, drop = FALSE] > 0) / n_its,
                       1 - colSums(as.matrix(comp_diff)[, 1:nC, drop = FALSE] > 0) / n_its)) * 2

    # Posterior mode and 95% HPD interval summaries
    cc1_res  <- paste(round(posterior.mode(comp_res),  dec_PM), " (",
                      round(HPDinterval(comp_res)[,  1], dec_PM), ", ",
                      round(HPDinterval(comp_res)[,  2], dec_PM), ")", sep = "")
    cc1_diff <- paste(round(posterior.mode(comp_diff), dec_PM), " (",
                      round(HPDinterval(comp_diff)[, 1], dec_PM), ", ",
                      round(HPDinterval(comp_diff)[, 2], dec_PM), ")", sep = "")

    # Output table columns:
    #   "Posterior Mode (CI)"             real dissimilarity between predicted communities
    #   "Difference from Null ..."        real minus null dissimilarity; positive = focal
    #                                     variable drives more turnover than chance
    #   "pMCMC"                           probability that the difference is zero (two-tailed)
    cc1 <- data.frame(
      "Variable"                                 = composition_var,
      "Community Comparison"                     = colnames(comp_res),
      "Posterior Mode (CI)"                      = cc1_res,
      "Difference from Null Posterior Mode (CI)" = cc1_diff,
      "pMCMC"                                    = round(as.numeric(cc1_p), 3),
      check.names = FALSE
    )

    community <- rbind(community, cc1)
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
    
    nG    <- dim(fixed_mod)[2]
    fe1_p <- calc_pmcmc(fixed_mod, pvalues_species, nG)
    
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
      fe2<-pairwise.diffs(fixed_mod,nF=nG)
      fixed=rbind(fe1,fe2)
      fixed = fixed %>% dplyr::filter(Fixed_Effects %in% c(fixed_names,fixed_diffinc_species) == T) %>% dplyr::select(Fixed_Effects,Estimates,pMCMC)
    }
    
    #****************************************************
    #Random effects ----
    #****************************************************
    if(Include_random_species == "yes") {
      # Pre-allocate: rows = combined MCMC samples, cols = species × random effects
      n_spp_samples <- nrow(do.call(rbind, post_model$Omega[[1]]))
      random_mod    <- matrix(nrow = n_spp_samples,
                              ncol = model$ns * length(post_model$Omega))

      for(i in seq_along(post_model$Omega)){
        tmp1 <- do.call(rbind, post_model$Omega[[i]])
        # Vectorized diagonal extraction: each row of tmp1 is a flattened ns×ns
        # covariance matrix; extract the diagonal (per-species variances) for each sample
        tmp2 <- t(apply(tmp1, 1, function(row) diag(matrix(row, nrow = model$ns, ncol = model$ns))))
        col_idx <- ((i - 1) * model$ns + 1):(i * model$ns)
        random_mod[, col_idx] <- tmp2
      }
      random_mod <- as.mcmc(random_mod)

      # Name columns: random effect label : species name
      if(is.null(random_names_species)) {
        random_names_species <- unlist(lapply(model$rLNames, function(nm) paste(nm, model$spNames, sep = ": ")))
      }
      colnames(random_mod) <- random_names_species
    
      #Summary of variance components
      randomVarspp = paste(round(posterior.mode(random_mod),dec_PM)," (",round(HPDinterval(random_mod)[,1],dec_PM), ", ",round(HPDinterval(random_mod)[,2],dec_PM),")",sep="")
      randomVarspp = data.frame("Random Effects"=colnames(random_mod),"Posterior Mode (CI)"=randomVarspp,"-"="",check.names=FALSE)
      
    } else  {
    }
    
    #****************************************************
    ##Excel output: fixed and random effects for species ----
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
    
    #table title
    header=data.frame(col1=c(""),col2=c(""),col3=c(""))
    colnames(header)<-c(paste(title,": species estimates",sep=" "),"","")

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
    
    #Random effects: variances

    #Should they be outputted or not
    if(Include_random_species == "yes") {
      writeData(workbook, sheet3, randomVarspp, startCol = 1, startRow = row_nums, headerStyle = hs2)
      row_nums <- row_nums + dim(randomVarspp)[1] + 2
    }

    bolding <- createStyle(textDecoration = "bold")
    conditionalFormatting(workbook, sheet3, cols = 3, rows = 1:10000, rule = "<0.05", style = bolding)

    #****************************************************
    #Variance Partitioning (% variation) of explained by fixed and random effects for each species  ----
    #****************************************************
    # Reuses VP_all computed in Section 1 — no second call to computeVariancePartitioning
    VP <- as.data.frame(VP_all)
    VP <- VP |> mutate(across(everything(), ~round(., 4) * 100))
    VP <- data.frame("Variance Partitioning" = rownames(VP), "%" = VP, check.names = FALSE)
    colnames(VP) <- sub("%\\.", "% ", colnames(VP))   # fixed regex: escape the literal dot
    rownames(VP) <- NULL

    #****************************************************
    #Species Fit statistics ----
    #****************************************************
    # Reuses model_preds computed in Section 1 — no second call to computePredictedValues
    model_fit_spp <- as.data.frame(evaluateModelFit(model, predY = model_preds))
    
    #Round fit stats
    model_fit_spp = model_fit_spp %>% mutate(across(everything(), ~round(., 2)))
    model_fit_spp = data.frame(Species=model$spNames,model_fit_spp)

    #****************************************************
    ##Excel output: variance partitioning for species ----
    #****************************************************
    if(VP_species == "include") {
    #Create new sheet 
    sheet4 = paste(sheet,"species_VP",sep="_")
    start_row = 1
    addWorksheet(workbook, sheet4)
    
    #table title
    header=data.frame(col1=c(""),col2=c(""),col3=c(""))
    colnames(header)<-c(paste(title,": species variances",sep=" "),"","")

    writeData(workbook, sheet4, header, startCol = 1, startRow = start_row,headerStyle = hs1)

    #Variance Partitioning
      writeData(workbook, sheet4, VP, startCol = 1, startRow = start_row+dim(header)[1],headerStyle = hs2)
      row_nums = start_row+dim(header)[1]+dim(VP)[1]+2


    writeData(workbook, sheet4, "Fit Statistics", startCol = 1, startRow = row_nums,headerStyle = hs2)
    addStyle(workbook, sheet4, style = hs2, rows = row_nums, cols = 1, gridExpand = TRUE)

    writeData(workbook, sheet4, model_fit_spp, startCol = 1, startRow = row_nums+1,headerStyle = hs2)
    #If species estimate != "include" then just returns a workbook with estimates averaged across species   
    return(workbook)
  }
  }
  return(workbook)
}


#function for extracting df from xl workbook
xl_2_df = function(xltab,sheet=NULL){
  df<-readWorkbook(xltab,sheet=sheet,startRow = 2)
  colnames(df)<-gsub("[.]"," ",colnames(df))
  rownames(df)<-NULL
  return(df)
}

#===========================================================
#function for df to md table that can handle html and other formats e.g. word ----
#===========================================================
hmsc_md <- function(data,stats=FALSE) {
  pacman::p_load(flextable,officer)

##Output if presenting mixed model stats: bolding significant values and highlighting headings
  if(stats == TRUE){
        #bold rows if less than 0.05 excluding cells with brackets which will be random effect
        rows_bold = data |> mutate(signif = !grepl("\\(", pMCMC) & pMCMC < 0.05,
                               signif = ifelse(is.na(signif),"FALSE",signif)) |> pull(signif)
  } else {
    #Other tables 
  } 
  
  #Replace any NAs with blanks
  data[is.na(data)] <- ""

  #Output if html format
  if (knitr::is_html_output()) {
  
  if(stats == TRUE){
    pacman::p_load(kableExtra)
      kbl(data, align = "l", digits = 3) |>
        kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11) |>
        row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Fixed Effect",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Trait Effects",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Random Effects",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Variance Partitioning",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Phylogenetic Effects",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Fit Statistics",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(nrow(data), extra_css = "border-bottom: 1px solid;margin-bottom:1000px") |> 
        column_spec(column=which(names(data) == "pMCMC"), bold =rows_bold) |> 
        row_spec(1:nrow(data), extra_css = "height: 1em; white-space: nowrap;") |> 
        column_spec(1:ncol(data), width = "auto")
    }
    #Output if not mixed model output or not stats testing
    else  {  
      kbl(data,align = "l")  |> 
        kable_styling(bootstrap_options = c("hover", "condensed"),html_font="helvetica",font_size = 11,fixed_thead = T) |>
        row_spec(0, bold=T,background="#E7E5E5", extra_css = "border-top: 1px solid; border-bottom: 1px solid") |>
        row_spec(nrow(data), extra_css = "border-bottom: 1px solid;margin-bottom:1000px") |>
        row_spec(1:nrow(data), extra_css = "height: 1em; white-space: nowrap;") |> 
        row_spec(grep("^Fixed Effect",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Trait Effects",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Random Effects",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Variance Partitioning",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Phylogenetic Effects",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        row_spec(grep("^Fit Statistics",data[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid") |> 
        column_spec(1:ncol(data), width = "auto")
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

    # Trait Effects
    trait_row <- grep("^Trait Effects", data[,1])
    ft  <- ft |> bold(i = trait_row, part = "body", bold = TRUE) |> 
                                  bg(i = trait_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = trait_row-1, part = "body", border = b) |> 
                                  hline(i = trait_row, part = "body", border = b)

    # Random
    random_row <- grep("^Random", data[,1])
    ft  <- ft |> bold(i = random_row, part = "body", bold = TRUE) |> 
                            bg(i = random_row, part = "body", bg   = "#E7E5E5") |> 
                            hline(i = random_row-1, part = "body", border = b) |> 
                            hline(i = random_row, part = "body", border = b)
    # Variance partitioning
    vp_row <- grep("^Variance Partitioning", data[,1])
    ft  <- ft |> bold(i = vp_row, part = "body", bold = TRUE) |> 
                                  bg(i = vp_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = vp_row-1, part = "body", border = b) |> 
                                  hline(i = vp_row, part = "body", border = b)
    # Phylogenetic Effects
    phy_row <- grep("^Phylogenetic Effects", data[,1])
    ft  <- ft |> bold(i = phy_row, part = "body", bold = TRUE) |> 
                                  bg(i = phy_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = phy_row-1, part = "body", border = b) |> 
                                  hline(i = phy_row, part = "body", border = b)
     # Fit statistics
    fit_row <- grep("^Fit Statistics", data[,1])
    ft  <- ft |> bold(i = fit_row, part = "body", bold = TRUE) |> 
                                  bg(i = fit_row, part = "body", bg   = "#E7E5E5") |> 
                                  hline(i = fit_row-1, part = "body", border = b) |> 
                                  hline(i = fit_row, part = "body", border = b)

    #Overall formatting
    ft <- ft |>
          theme_vanilla() |>
          flextable::fontsize(size = 10, part = "header") |>
          flextable::fontsize(size = 8, part = "body") |>
          bold(part = "header") |>
          #Bold pMCMC less than 0.05
          bold(i = ~ pMCMC < 0.05, j = "pMCMC", bold = TRUE, part = "body") |>
          bg(part = "header", bg = "#E7E5E5") |>
          flextable::border(border.top = fp_border(color = "black", width = 1), part = "header") |>
          flextable::border(border.bottom = fp_border(color = "black", width = 1), part = "header") |>
          flextable::border(border.bottom = fp_border(color = "black", width = 1), i = nrow(data)) |>
          set_table_properties(width=1,layout = "autofit",opts_html = list(scroll = list(height = "1000px", freeze_first_column = TRUE)),opts_word = list(keep_with_next = TRUE)) |>
          autofit()
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
          set_table_properties(width=1,layout = "autofit",opts_html = list(scroll = list(height = "1000px", freeze_first_column = TRUE)),opts_word = list(keep_with_next = TRUE)) |>
          autofit()
    }
    }
 }


#function for df to Rmd table
hmsc_md2 = function(df){
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
    row_spec(grep("^Fit Statistics",df[,1]), bold=T,background="#E7E5E5",extra_css = "border-top: 1px solid; border-bottom: 1px solid")%>%
    row_spec(nrow(df), extra_css = "border-bottom: 1px solid;margin-bottom:1000px")}
