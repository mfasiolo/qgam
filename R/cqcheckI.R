##########################
#' Interactive visual checks for additive quantile fits
#' 
#' @description Given an additive quantile model, fitted using \code{qgam}, \code{cqcheck2DI} provides some interactive
#'              2D plots that allow to check what proportion of responses, \code{y}, falls below the fitted quantile.
#'              This is an interactive version of the \code{cqcheck} function.
#'  
#' @param obj the output of a \code{qgam} call. 
#' @param v if a 1D plot is required, \code{v} should be either a single character or a numeric vector. In the first case
#'          \code{v} should be the names of one of the variables in the dataframe \code{X}. In the second case, the length
#'          of \code{v} should be equal to the number of rows of \code{X}. If a 2D plot is required, \code{v} should be 
#'          either a vector of two characters or a matrix with two columns.  
#' @param X a dataframe containing the data used to obtain the conditional quantiles. By default it is NULL, in which
#'          case predictions are made using the model matrix in \code{obj$model}.
#' @param y vector of responses. Its i-th entry corresponds to the i-th row of X.  By default it is NULL, in which
#'          case it is internally set to \code{obj$y}.
#' @param run if TRUE (default) the function produces an interactive plot, otherwise it returns the corresponding shiny app.    
#' @param width the width of the main plot. Default is "100\%".
#' @param height width the width of the main plot. Default is "680px".
#' @return Simply produces an interactive plot.
#' @details This is an interactive version of the \code{cqcheck}, see \code{?cqcheck} for details. The main interactive
#'          feature is that one can select an area by brushing, and then double-click to zoom in. In the 1D case the vertical 
#'          part of the selected area is not use: we zoom only along the x axis. Double-clicking without brushing zooms out. 
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>. 
#' @examples
#' \dontrun{ 
#' #######
#' # Example 1: Bivariate additive model y~1+x+x^2+z+x*z/2+e, e~N(0, 1)
#' #######
#' library(qgam)
#' set.seed(15560)
#' n <- 1000
#' x <- rnorm(n, 0, 1); z <- rnorm(n)
#' X <- cbind(1, x, x^2, z, x*z)
#' beta <- c(0, 1, 1, 1, 0.5)
#' y <- drop(X %*% beta) + rnorm(n) 
#' dataf <- data.frame(cbind(y, x, z))
#' names(dataf) <- c("y", "x", "z")
#' 
#' #### Fit a constant model for median
#' qu <- 0.5
#' fit <- qgam(y~1, qu = qu, err = 0.05, data = dataf)
#' 
#' # Look at what happens along x: clearly there is non linear pattern here
#' cqcheckI(obj = fit, v = c("x"), X = dataf, y = y) 
#' 
#' #### Add a smooth for x
#' fit <- qgam(y~s(x), qu = qu, err = 0.05, data = dataf)
#' cqcheckI(obj = fit, v = c("x"), X = dataf, y = y) # Better!
#' 
#' # Lets look across across x and z. As we move along z (x2 in the plot) 
#' # the colour changes from green to red
#' cqcheckI(obj = fit, v = c("x", "z"), X = dataf, y = y)
#' 
#' # The effect look pretty linear
#' cqcheckI(obj = fit, v = c("z"), X = dataf, y = y)
#' 
#' #### Lets add a linear effect for z 
#' fit <- qgam(y~s(x)+z, qu = qu, err = 0.05, data = dataf)
#' 
#' # Looks better!
#' cqcheckI(obj = fit, v = c("z"))
#' 
#' # Lets look across x and y again: green prevails on the top-left to bottom-right
#' # diagonal, while the other diagonal is mainly red.
#' cqcheckI(obj = fit, v = c("x", "z"))
#' 
#' ### Maybe adding an interaction would help?
#' fit <- qgam(y~s(x)+z+I(x*z), qu = qu, err = 0.05, data = dataf)
#' 
#' # It does! The real model is: y ~ 1 + x + x^2 + z + x*z/2 + e, e ~ N(0, 1)
#' cqcheckI(obj = fit, v = c("x", "z"))
#' }
#'
cqcheckI <- function(obj, v, X = NULL, y = NULL, run = TRUE, width = "100%", height = "680px")
{
  
  #### Set up
  if( is.null(X) ){ 
    X <- obj$model
    if( is.null(y) ){ y <- obj$y }
  } else {
    if( is.null(y) ){ stop("If you provide X you must provide also the corresponding vector of responses y") }
  }
  
  if( length(y)!=nrow(X) ){ stop("length(y)!=nrow(X)") }
  
  ####### Setting up 1D and 2D cases
  if( is.character(v) ){ # Name(s) of variable(s) in X provided OR ...
    if(length(v) == 1){ ## 1D CASE ##
      if( !(v %in% names(X)) ) stop("(v %in% names(X)) == FALSE")
      x1 <- X[[v]]
      x2 <- NULL
    } else {
      if(length(v) == 2){ ## 2D CASE ##
        if( !(v[1] %in% names(X)) ) stop("(v[1] %in% names(X)) == FALSE")
        if( !(v[2] %in% names(X)) ) stop("(v[2] %in% names(X)) == FALSE")
        x1 <- X[[v[1]]]
        x2 <- X[[v[2]]]
      } else { 
        stop("If is.character(v)==TRUE, then length(v) should be either 1 or 2.") 
      }
    }
  } else { # ... actual numeric value of the variable(s) provided 
    if(is.vector(v)){ ## 1D CASE ##
      x1 <- v
      x2 <- NULL
      if(length(v) != nrow(X)){ stop("length(v) != ncol(X)") }
    } else {
      if(is.matrix(v)){ ## 2D CASE ##
        if(ncol(v)!=2){ stop("In the 2D case, v should be a matrix with 2 columns or a vector of 2 characters") } 
        x1 <- v[ , 1]
        x2 <- v[ , 2]
      } else {
        stop("In the 2D case, v should be a matrix with 2 columns or a vector of 2 characters")
      }
    }
  } 
  
  out <- if( is.null(x2) ){ # One dimensional OR ...
    .cqcheck1DI(.obj = obj, .x1 = x1, .X = X, .y = y, .width = width, .height = height)
  } else { # ... two dimensional case
    .cqcheck2DI(.obj = obj, .x1 = x1, .x2 = x2, .X = X, .y = y, .width = width, .height = height)
  }
  
  if( run ){ return(runApp(out)) } else { return(out) }
}

#### Internal function for 1D plot
.cqcheck1DI <- function(.obj, .x1, .X, .y, .width, .height)
{
  # User interface
  ui <- fluidPage(
    sidebarPanel(
      numericInput('nbin', 'Num. bins', 10, min = 1, max = Inf), 
      numericInput('lev', 'Sign. lev.', 0.05, min = 0, max = 1),
      width = 2
    ),
    mainPanel(
      h4("Brush and double-click to zoom in. Double-click to zoom out."),
      plotOutput("plot1", 
                 dblclick = "plot1_dblclick",
                 hover = "plot1_hover",
                 brush = brushOpts(
                   id = "plot1_brush",
                   resetOnNew = TRUE), 
                 width = .width,
                 height = .height),
      
      verbatimTextOutput("info")
    )
  )
  
  # Server side
  server <- function(input, output) {
    
    # Control panel inputs
    ranges <- reactiveValues(x = NULL, y = NULL)
    nbin <- reactive({ input$nbin })
    lev <- reactive({ input$lev })
    
    myPlot <- function(brush = NULL)
    {
      if (!is.null(brush)) {
        good <- which(.x1 > brush$xmin & .x1 < brush$xmax)
        out <-  cqcheck(obj = .obj, v = .x1[good], X = .X[good, ], y = .y[good], nbin = nbin(), lev = lev())
      } else {
        out <- cqcheck(obj = .obj, v = .x1, X = .X, y = .y, nbin = nbin(), lev = lev())
      }
      return( out )
    }
    
    # Initial plot
    output$plot1 <- renderPlot({ myPlot() })
    
    # Update if double click or double click + brush
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      output$plot1 <- renderPlot({ myPlot(brush = brush) })
    })
    
    # Print some info
    output$info <- renderText({
      x_str <- function(e) {
        if(is.null(e)) return("NULL\n")
        paste0("x=", round(e$x, 1), "\n")
      }
      x_range_str <- function(e) {
        if(is.null(e)) return("NULL\n")
        paste0("xmin=", round(e$xmin, 1), ", xmax=", round(e$xmax, 1))
      }
      
      paste0(
        "hover: ", x_str(input$plot1_hover),
        "brush: ", x_range_str(input$plot1_brush)
      )
    })
  }
  return( shinyApp(ui, server) )
}



#### Internal function for 2D plot
.cqcheck2DI <- function(.obj, .x1, .x2, .X, .y, .width, .height)
{
  # User interface
  ui <- fluidPage(
    sidebarPanel(
      numericInput('nbin1', 'N. bins x1', 10, min = 1, max = Inf), 
      numericInput('nbin2', 'N. bins x2', 10, min = 1, max = Inf), 
      numericInput('lev', 'Sign. lev.', 0.05, min = 0, max = 1), 
      checkboxInput('scatter', 'Add scatter', value = FALSE),
      width = 2
    ),
    mainPanel(
      h4("Brush and double-click to zoom in. Double-click to zoom out."),
      plotOutput("plot1", 
                 dblclick = "plot1_dblclick",
                 hover = "plot1_hover",
                 brush = brushOpts(
                   id = "plot1_brush",
                   resetOnNew = TRUE), 
                 width = .width,
                 height = .height),
      
      verbatimTextOutput("info")
    )
  )
  
  # Server side
  server <- function(input, output) {
    
    # Control panel inputs
    ranges <- reactiveValues(x = NULL, y = NULL)
    nbin1 <- reactive({ input$nbin1 })
    nbin2 <- reactive({ input$nbin2 })
    lev <- reactive({ input$lev })
    scatter <- reactive({ input$scatter })
    
    myPlot <- function(brush = NULL)
    {
      if (!is.null(brush)) {
        good <- which(.x1 > brush$xmin & .x1 < brush$xmax & .x2 > brush$ymin & .x2 < brush$ymax)
        out <-  cqcheck(obj = .obj, v = cbind(.x1[good], .x2[good]), X = .X[good, ], y = .y[good], 
                        nbin = c(nbin1(), nbin2()), scatter = scatter(), lev = lev())
      } else {
        out <- cqcheck(obj = .obj, v = cbind(.x1, .x2), X = .X, y = .y, 
                       nbin = c(nbin1(), nbin2()), scatter = scatter(), lev = lev())
      }
      return( out )
    }
    
    # Initial plot
    output$plot1 <- renderPlot({ myPlot() })
    
    # Update if double click or double click + brush
    observeEvent(input$plot1_dblclick, {
      brush <- input$plot1_brush
      output$plot1 <- renderPlot({ myPlot(brush = brush) })
    })
    
    # Print some info
    output$info <- renderText({
      xy_str <- function(e) {
        if(is.null(e)) return("NULL\n")
        paste0("x1=", round(e$x, 1), ", x2=", round(e$y, 1), "\n")
      }
      xy_range_str <- function(e) {
        if(is.null(e)) return("NULL\n")
        paste0("x1min=", round(e$xmin, 1), ", x1max=", round(e$xmax, 1),
               ", x2min=", round(e$ymin, 1), ", x2max=", round(e$ymax, 1))
      }
      
      paste0(
        "hover: ", xy_str(input$plot1_hover),
        "brush: ", xy_range_str(input$plot1_brush)
      )
    })
  }
  
  return( shinyApp(ui, server) )
}