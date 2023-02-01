

drop_feature_correlated <- function(data){
  
  cor_matrix_rm <- round(cor(data), 2)                 # Modify correlation matrix
  cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
  diag(cor_matrix_rm) <- 0
  
  
  data_new <- data[ , !apply(cor_matrix_rm, 2, function(x) any(x>=0.95))]
  return(data_new)
}
