MCMCeigen = function(VCV,traits){
  result = data.frame(iteration=as.numeric(),
                      eigen_vector=as.numeric(),
                      eigen_value=as.numeric(),
                      sum_values=as.numeric())
  
  for(i in 1:dim(VCV)[1]){ 
    tmp = eigen(matrix(VCV[i,],traits,traits))
    tmp = data.frame(iteration=i,
                     eigen_vector=seq(1:dim(tmp$vectors)[1]),
                     eigen_value=tmp$values,
                     sum_values=sum(tmp$values))
    result = rbind(result,tmp)
  }
  result$prop_values = result$eigen_value / result$sum_values
  
  sum_res = result %>% pivot_wider(id_cols = iteration,names_from = eigen_vector, values_from = prop_values) 
  sum_res = as.mcmc(sum_res[,-1])
  sum_res = data.frame(eigen_vector=unique(result$eigen_vector),post_mode=posterior.mode(sum_res),lowerCI=HPDinterval(sum_res)[,1],upperCI=HPDinterval(sum_res)[,2])
  return(list(result, sum_res))
}