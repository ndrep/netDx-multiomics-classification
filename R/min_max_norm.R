min_max_norm <- function(x){
  (x - min(x)) / (max(x) - min(x))
}