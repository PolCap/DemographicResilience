
# Functions to simulate matrix population models

generate_survival_curve <- function(type, values = 10, exp = 10, to_zero = FALSE){
  
  
  # Check if type is valid
  if(type < 1 || type > 3){
    cat("---Argument type is invalid.")
    return()
  }
  
  x <- seq(0,1, length = values + 1)
  
  survival_three <- rev(x^exp)
  survival_two <- rev(x)
  survival_one <- rev((-1*survival_three)+1)
  
  if(type == 2){
    
    survival_curve <- data.frame(cbind(x, survival_two))
    names(survival_curve)[2] <- c("survival_curve")
    
    return(survival_curve)
    
  }
  
  survival_curve <- c()
  
  if(type >= 1 && type < 2){
    
    weight_curve_two <- type - 1
    
    for(i in 1:values + 1){
      
      survival_curve[i] <- weighted.mean(c(survival_one[i], survival_two[i]), c(1 - weight_curve_two, weight_curve_two))
      
    }
    
    
  } else if(type > 2 && type <= 3){
    
    weight_curve_three <- type - 2
    
    for(i in 1:values + 1){
      
      survival_curve[i] <- weighted.mean(c(survival_two[i], survival_three[i]), c(1 - weight_curve_three, weight_curve_three))
      
    }
    
  } else{
    
    cat("---Something went wrong.")
    
    return()
    
  }
  
  survival_curve <- (10^survival_curve)/10
  
  survival_curve[1] <- 1
  
  # Now that the survival curve (lx) has been generated, now we will calculate the vital rate (px)
  
  px <- rep(0, times = length(survival_curve))
  
  px[1] <- 1
  
  for(i in 2:length(px)){
    
    px[i] <- survival_curve[i]/survival_curve[i-1]
    
  }
  
  survival_values <- data.frame(cbind(x, survival_curve, px)) 
  
  if(to_zero == TRUE){
    
    return(survival_values)
  }
  
  survival_values <- survival_values[-nrow(survival_values), ]
  
  return(survival_values)
  
}

generate_reproduction_curve <- function(cutoff, values, slope, starting_value = 0.1){
  
  x <- seq(0,1, length = values)
  
  reproduction <- rep(0, times = values)
  
  cutoff_index <- ceiling(cutoff * values)
  
  reproduction_output <- reproduction[cutoff_index:values]
  
  reproduction_output[1] <- starting_value
  
  if(length(reproduction_output) > 1){
    
    for(i in 2:length(reproduction_output)){
      
      
      reproduction_output[i] <- reproduction_output[i-1] * slope
      
    }
    
  }
  
  reproduction[cutoff_index:values] <- reproduction_output
  
  reproduction_curve <- data.frame(cbind(x, reproduction)) 
  
  return(reproduction_curve)
  
}

build_leslie <- function(survival, reproduction,
                         average_stasis_growth_ratio = FALSE,
                         random_stasis_growth_ratio = FALSE){
  
  matrix_A <- matrix(rep(0, times = length(survival)^2), 
                     nrow = length(survival), 
                     ncol = length(survival))
  matrix_U <- matrix_A
  matrix_F <- matrix_A
  
  matrix_F[1,] <- reproduction
  
  if(average_stasis_growth_ratio == FALSE &&
     random_stasis_growth_ratio == FALSE){
    
    for(i in 2:length(survival)){
      
      matrix_U[[i,i-1]] <- survival[i]
      
    }
    
    matrix_U[[length(survival), length(survival)]] <- survival[length(survival)]
    
    matrix_A <- matrix_U + matrix_F
    
    matrices <- list("matrix_A" = matrix_A,
                     "matrix_U" = matrix_U,
                     "matrix_F" = matrix_F)
    
    return(matrices)
    
  } else if(is.numeric(average_stasis_growth_ratio) == TRUE &&
            random_stasis_growth_ratio == FALSE){
    
    
    for(i in 2:length(survival)){
      
      matrix_U[[i,i-1]] <- survival[i] * (1-average_stasis_growth_ratio)
      matrix_U[[i-1,i-1]] <- survival[i] * (average_stasis_growth_ratio)
      
    }
    
    matrix_U[[length(survival), length(survival)]] <- survival[length(survival)]
    
    matrix_A <- matrix_U + matrix_F
    
    matrices <- list("matrix_A" = matrix_A,
                     "matrix_U" = matrix_U,
                     "matrix_F" = matrix_F)
    
    return(matrices)
    
  } else if(average_stasis_growth_ratio == FALSE &&
            random_stasis_growth_ratio == TRUE){
    
    
    for(i in 2:length(survival)){
      
      random_ratio <- runif(n = 1, min = 0, max = 1)
      
      matrix_U[[i,i-1]] <- survival[i] * (1-random_ratio)
      matrix_U[[i-1,i-1]] <- survival[i] * (random_ratio)
      
    }
    
    matrix_U[[length(survival), length(survival)]] <- survival[length(survival)]
    
    matrix_A <- matrix_U + matrix_F
    
    matrices <- list("matrix_A" = matrix_A,
                     "matrix_U" = matrix_U,
                     "matrix_F" = matrix_F)
    
    return(matrices)
    
  }
  
}

build_lefkovitch <- function(x=NA,
                             survival = NA,
                             reproduction = NA,
                             maturation = 1) {
  
  stage_name <- x

  # Useful indicies
  n_stage <- length(stage_name)
  stage_index <- seq_along(stage_name)
  
  # Build matrix for survival
  matrix_U <- mpmtools::subdiag(n_stage, (survival * maturation)[-n_stage]) +
    diag(survival * (1 - maturation))
  
  # Build matrix for maternity "season"
  mat_mat <- diag(n_stage)
  mat_mat[1,] <- mat_mat[1,] + reproduction
  
  # Build the final matrix  
  
  matrix_A <- matrix_U %*% mat_mat
  matrix_A <- dropzeros(matrix_A)
  matrix_F <- matrix_A
  matrix_F[-1,] <- 0
  matrix_U <- matrix_A
  matrix_U[1,] <- 0
  matrices <- list("matrix_A" = matrix_A,
                   "matrix_U" = matrix_U,
                   "matrix_F" = matrix_F)
  
  return(matrices)

}

