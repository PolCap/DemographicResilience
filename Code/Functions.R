# ---------------------------------------------------------------------------- #
# - FILE NAME:   Functions.R         
# - DATE:        12/08/2019
# - DESCRIPTION:  Functions to estimate the transient dynamics of populations in 
#                 compadre and comadre. 
# - AUTHORS:     Iain Stott & Pol Capdevila Lanzaco
# ---------------------------------------------------------------------------- #


return.time <- function(matA){
  dr(matA, return.time = T)$t
}

AgeAtMaturity <- function(matU, matR, startLife){
  lifeTimeRepEvents(matU, matR, startLife)$La
}


cdb_collapse <-  function(cdb, columns) {
  
  if ("MatrixComposite" %in% names(cdb) &&
      any(cdb$MatrixComposite == "Seasonal")) {
    warning("cdb contains rows with MatrixComposite == 'Seasonal'. This ",
            "method of collapsing is not suitable for seasonal matrices. ",
            "Consider removing prior to collapsing.", call. = FALSE)
  }
  
  # leave other validation to cdb_id
  
  # create a unique integer identifier for each group to be collapsed
  id_collapse <- cdb_id(cdb, columns)
  id_unique <- unique(id_collapse)
  
  # extract data slot
  dat <- cdb@data
  
  # remove any list-columns except for 'mat'
  col_mat <- names(dat) == "mat"
  col_list <- vapply(dat, class, "", USE.NAMES = FALSE) == "list"
  dat <- dat[,col_mat | !col_list]
  
  # coerce factor columns to character
  is_factor <- vapply(dat, is.factor, logical(1), USE.NAMES = FALSE)
  if (any(is_factor)) {
    j <- which(is_factor)
    for(i in j) dat[[i]] <- as.character(dat[[i]])
  }
  
  # collapse to list of 1-row data frames for each group
  coll_l <- lapply(id_unique,
                   FUN = CollapseWrapper,
                   dat = dat,
                   id_collapse = id_collapse)
  
  # bind data frames
  coll_dat <- do.call(plyr::rbind.fill, coll_l)
  
  new('CompadreDB',
      data = coll_dat,
      version = cdb@version)
}
CollapseWrapper <- function(id, dat, id_collapse) {
  dat <- dat[id_collapse %in% id,]
  CollapseFn(dat)
}


CollapseFn <- function(x) {
  mat <- mpm_mean(x$mat)
  d <- as_tibble(lapply(x[,-1], CollapseCol))
  d <- add_column(d, mat = list(mat), .before = 1)
  if (nrow(x) > 1) {
    d$MatrixComposite <- "Collapsed"
    if ("Lat" %in% names(x)) d$Lat <- mean(x$Lat)
    if ("Lon" %in% names(x)) d$Lon <- mean(x$Lon)
    d$SurvivalIssue <- max(colSums(mat@matU))
  }
  return(d)
}


CollapseCol <- function(col) {
  ifelse(length(unique(col)) == 1,
         unique(col),
         paste(unique(col), collapse = "; "))
}

Fecundity <-  function(matU, matF){
  vitalRates(matU, matF, weights = "SSD")$fec
}

# Function to convert post-reproductive census projection matrix models to pre-reproductive census models
# Modified from Jelbert et al. 2019 Nat Comms (requires RCompadre)

convert2pre <- function(matA, matF, matC, matU, problem){
  mat_corrected_pre <- list()
  if(problem=="pre"){mat_corrected_pre<-matA[[1]]}
  else{
    mat <- matA[[1]]
    matF <- matF[[1]]
    matC <- matC[[1]]
    matS<-matU[[1]]
    matrix_size <- nrow(mat)
    if(problem=="error"){
      newmat<-matS%*%matF+matC+matS
    }
    if(problem=="post"){ # persistent seedbank
    surv_vec<-apply(matS,2,sum)
    surv_mat<-matrix(surv_vec,nrow=matrix_size,ncol=matrix_size,byrow=T)
    newmat1<-matF/surv_mat
    newmat1[is.nan(newmat1)] = 0
    newmat<-matS%*%newmat1+matC+matS
    newmat[is.nan(newmat)] = 0}
  if(mat[1,1]==0){
      mat_corrected_pre<-newmat[2:matrix_size,2:matrix_size]
      }else{mat_corrected_pre<-newmat}
    if(any(is.infinite(mat_corrected_pre))){
      mat_corrected_pre <- NA}
    }
  return(list(mat_corrected_pre))
}
