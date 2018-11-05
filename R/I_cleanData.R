###
# Removing all NAs, unused variables and factor levels from data
#
.cleanData <- function(.dat, .form, .drop){
  
  if( inherits(.dat, "groupedData") ) { .dat <- as.data.frame( .dat ) }
  
  .vars <- .allVars1( interpret.gam(.form)$fake.formula )
  
  .dat <- na.omit( subset(.dat, select = .vars) )
  
  if( is.null(.drop) || .drop ) { .dat <- droplevels( .dat ) }
  
  return( .dat )
  
}