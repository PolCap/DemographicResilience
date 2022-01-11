# Function to simulate random matrices 

random_matrices <- function(dimension=NA, upd=FALSE){
  
  # Simulate x transitions from 0 to 1
  matUvector <- runif(dimension^2,0,1) 
  
  # Transform into a matrix
  matU <- matrix(matUvector,nrow=dimension) 
  
  # Remove the first row of the matU
  
  matU[1,2:dimension] <- 0  
  
  # Simulate stage-specific survival vector
  
  surv <- runif(dimension,0,1)
  
  # Remove shrinkage if subd is true
  
  if(isTRUE(upd)){matU[upper.tri(matU)] <- 0}
  
  # Make sure that columns add to 1
  
  matU <- apply(matU, 2, function(x) x/sum(x))
  
  # Penalise by stage-specific survival
  
  for(i in 1:dimension) {matU[,i]<-matU[,i]*surv[i]} 
  
  # Create a 0 matrix
  matF <- matU*0
  
  # Simulate random fecundity
  
  matF[1,2:dimension] <- rpois(dimension-1, lambda=sample.int(100, 1))
  
  # Create the Matrix A
  
  matA <- matU+matF
  
  matrices <- list("matrix_A" = matA,
                   "matrix_U" = matU,
                   "matrix_F" = matF)
  
  return(matrices)
}

