This repositories contains functions for phylogenetic analyses:

ASF_functions.R 
- processes output from various models placing node estimates in dataframes facilitating parent descendent values

MCMCglmmProc.R
- Processes output from MCMCglmm models into formatted tables.
- Takes the following arguments:

#model = MCMCglmm model
#response = list of responses (e.g. c(trait1,trait2))
#link = link functions used for each response variable (e.g. c("logit","gaussian")). Changes the calculation of ICCs by adding distribution variances.
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
#variances = names of variance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed. If left as NULL then taken from the model object.
#covariances = names of covariance terms in VCV object either indices or written. Must be as they appear in the model object not the renamed. If left as NULL then they are not outputted.
#randomvar_names = what you want random effect variances to be called. Note reserved term animal is renamed to phylogeny
#randomcovar_names = what you want random effect covariances to be called. Note for correlations to be calculated this has to be same as variance1 : variance2
#Include_random = should random effects be included in output
#Padding = space between tables when outputting multiple models to same sheet
#dec_PM = number of decimals given for posterior mode and CIs of fixed and random effects
#link = allows logit & probit, the default is gaussian, and is used to calculate ICCs correctly. Gaussian can also be used for poisson (log) models to produce ICCs on expected scale but NOT for data scale estimates. See de Villemereuil 2016 Genetics & QGglmm package for details.
#For Multi-response models can provide a list of link functions (e.g. c("gaussian","logit")) corresponding to each response trait
#responses = specify response variables can take multiple values for multi response
#pvalues = exclusion of pMCMC values for fixed effects - "exclude", "include" or "c(?,?...)" giving list of which p values to exclude. Note pMCMC will still be calculated for fixed effect comparisons
#levels = if at.level notation is used how many levels are there. Default is 0.
