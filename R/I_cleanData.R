###
# Removing all NAs, unused variables and factor levels from data
#
.cleanData <- function(.dat, .form, .drop){
  
  if( inherits(.dat, "groupedData") ) { .dat <- as.data.frame( .dat ) }
  
  .vars <- .allVars1( interpret.gam(.form)$fake.formula )
  
  # Data is a list, but not a data.frame hence it needs special treament
  # NB: assumption here is that .dat contains 1 and only 1 data.frame and 
  # remaining entries are matrices
  if( is.list(.dat) && !is.data.frame(.dat) ){
    
    # Keep only elements that are data.frames or that appear in model formula
    .dat <-  .dat[ which( sapply(.dat, is.data.frame) | (names(.dat) %in% .vars) ) ]
    
    # Check if there are matrices (e.g. for functional effects)
    .matVar <- which( !sapply(.dat, is.data.frame) )
    
    # No matrices in the list, we do the check only on the data.frame element
    if( !length(.matVar) ){
      return( .cleanData(.dat = .dat[ sapply(.dat, is.data.frame) ][[1]], 
                         .form = .form, .drop = .drop) )
    }
    
    .datI <- subset(.dat[ -.matVar ][[1]], 
                    select = .vars[!(.vars %in% names(.dat))] ) # Standard part of data
    .datM <- .dat[ .matVar ]                                    # List of matrices
    
    # Find rows with NAs in two parts
    .badI <- attr(na.omit(.datI), "na.action")
    .badM <- attr(na.omit(do.call("cbind", .datM)), "na.action")
    
    # Now remove all bad rows from all elements of .dat
    .badAll <- union(.badI, .badM) 
    
    if( !is.null(.badAll) ){
      
      .datI <- .datI[-.badAll, ]
      if( is.null(.drop) || .drop ) { .datI <- droplevels( .datI ) }
      
      .datM <- lapply(.datM, function(.X) .X[-.badAll, ])
      
      .datO <- c(list(.datI), .datM)
      
    }
    
  } else{
    
    .datO <- na.omit( subset(.dat, select = .vars) )
    
    if( is.null(.drop) || .drop ) { .dat <- droplevels( .dat ) }
    
  }
  
  return( .datO )
  
}

####
# Test case
####
# .dat <- list(data.frame("x" = 1:10, "y" = 1:10, "u1" = rnorm(10)),
#              "z1" = matrix(1:30, 10, 3),
#              "z2" = matrix(1:30, 10, 3),
#              "u2" = matrix(1:30, 10, 3))
# .dat[[1]]$y[1] <- NA
# .dat$z1[3, 2] <- NaN
# .dat$z2[10, 2] <- NaN
# # .vars <- c("x", "y", "z1", "z2")
# 
# # The 1st, 3rd and 10th rows should be removed, and the u1 and u2 variables should disappear
# qgam:::.cleanData(.dat = .dat, .form = y ~ s(x) + s(y, z1) + z2, .drop = TRUE)