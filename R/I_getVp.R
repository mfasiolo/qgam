.getVp <- function(.fit, .obj, .lsp, .lpi)
{
  .Vp <- if(inherits(.obj$family, "general.family")){
    .fit$gcv.ubre <- as.numeric(.fit$REML) 
    .fit$outer.info <- NULL
    .fit$sp <- exp(.lsp)
    .fit$scale.estimated <- FALSE
    .fit$scale <- 1
    .fit$method <- "REML"
    .Vp <- mgcv:::gam.fit5.post.proc(.fit,.obj$Sl,.obj$L,.obj$lsp0,.obj$S,.obj$off)$Vb
    .Vp <- .Vp[.lpi[[1]], .lpi[[1]]]
  } else {
    .Vp <- .fit$Vb
  } 
  return( .Vp )
}

