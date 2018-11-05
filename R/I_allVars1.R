## 
# Version of all.vars that doesn't split up terms like x$y into x and y
# Copied from all.vars1 function in mgcv 1.18-24
#
.allVars1 <- function(form){
  vars <- all.vars(form)
  vn <- all.names(form)
  vn <- vn[vn%in%c(vars,"$","[[")] ## actual variable related names
  if ("[["%in%vn) stop("can't handle [[ in formula")
  ii <- which(vn%in%"$") ## index of '$'
  if (length(ii)) { ## assemble variable names
    vn1 <- if (ii[1]>1) vn[1:(ii[1]-1)]
    go <- TRUE
    k <- 1
    while (go) {
      n <- 2; 
      while(k<length(ii) && ii[k]==ii[k+1]-1) { k <- k + 1;n <- n + 1 }
      vn1 <- c(vn1,paste(vn[ii[k]+1:n],collapse="$"))
      if (k==length(ii)) {
        go <- FALSE
        ind <- if (ii[k]+n<length(vn)) (ii[k]+n+1):length(vn) else rep(0,0) 
      } else {
        k <- k +  1
        ind <- if (ii[k-1]+n<ii[k]-1) (ii[k-1]+n+1):(ii[k]-1) else rep(0,0)
      }
      vn1 <- c(vn1,vn[ind])
    }
  } else vn1 <- vn
  return( vn1 )
}
