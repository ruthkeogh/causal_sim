#------------------------------------------------------
#------------------------------------------------------
# functions used: expit, dp3, table.func, addlinetoplot, addlinetoplot.2
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#' expit function
expit=function(x){exp(x)/(1+exp(x))}

#' dp3: Function for printing x to 3 decimal places
dp3=function(x){sprintf("%.3f",x)}

#' table.func Function for printing "estimate (empirical standard deviation)" used to create Table 2
#' @param s Estimate (a vector of 5 numbers)
#' @param  sd Standard deviation of estimate (a vector of 5 numbers)
table.func=function(s,sd){
  c(paste0(dp3(s[1])," (",dp3(sd[1]),")"),
    paste0(dp3(s[2])," (",dp3(sd[2]),")"),
    paste0(dp3(s[3])," (",dp3(sd[3]),")"),
    paste0(dp3(s[4])," (",dp3(sd[4]),")"),
    paste0(dp3(s[5])," (",dp3(sd[5]),")"))
}

#' addlinetoplot: Function for adding a line to a ggplot using geom_line
#' allowing the user to specify the colour of the line
#' @param dataset Data frame corresponding to 'data' argument used in geom_line
#' @param varx x variable corresponding 'x' argument used in aes_string component of geom_line
#' @param vary y variable corresponding 'y' argument used in aes_string component of geom_line
#' @param vcol colour to be used for the line, corresponding to the 'colour' used in aes_string component of geom_line
addlinetoplot <- function(dataset, varx, vary,vcol) { 
  list(geom_line(data=dataset, aes_string(x=varx, y=vary,group=1,colour=vcol),size=1.5))
}

#' addlinetoplot.2: Function for adding a line to a ggplot using geom_line
#' allowing the user to specify the colour and linetype (solid, dashed, etc) of the line
#' @param dataset Data frame corresponding to 'data' argument used in geom_line
#' @param varx x variable corresponding 'x' argument used in aes_string component of geom_line
#' @param vary y variable corresponding 'y' argument used in aes_string component of geom_line
#' @param vcol colour to be used for the line, corresponding to the 'colour' used in aes_string component of geom_line
#' @param vline linetype to be used for the line, corresponding to the 'linetype' used in aes_string component of geom_line
addlinetoplot.2 <- function(dataset, varx, vary,vcol,vline) { 
  list(geom_line(data=dataset, aes_string(x=varx, y=vary,group=1,colour=vcol,linetype=vline),size=1.5))
}
