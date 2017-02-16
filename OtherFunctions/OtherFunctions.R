

## Quick function to sort a HH vector of within household data by relate variable
quick_sort_HH <- function(X=x,P=p,R=r){# X =vector, p =no of variables, r =col index of relate var 
  XX = matrix(X,ncol=P,byrow = TRUE)
  XX = XX[order(XX[,R]),]
  return(matrix(t(XX),nrow=1))
}