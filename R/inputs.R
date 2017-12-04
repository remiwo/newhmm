###################################################
### window size
###################################################

is_window <- function(x) {
  assert_that(is.numeric(x), length(x)==1)
  (x==3) | (x==5) | (x==7) | (x==9)
}

on_failure(is_window) <- function(call, env) {
  paste0(deparse(call$x)," is not a proper window size. It must be 3 or 5 or 7 or 9.")
}

###################################################
### protein sequence
###################################################

is_protein <- function(x) {
  assert_that(is.string(x), not_empty(x), is.character(x), nchar(x)<=30)
  (grepl('U',x)==F) & (grepl('O',x)==F)  & (grepl('J',x)==F) & (grepl('V',x)==F) & (grepl('X',x)==F) & (grepl('B',x)==F)
}

on_failure(is_protein) <- function(call,env) {
  paste0(deparse(call$x)," is not a valid peptide sequence. See documentation for details")
}

###################################################
### turn angle
###################################################

is_turn_angle <-function(x) {
  assert_that(is.numeric(x), length(x)==1)
  (x>=0) & (x<=2*pi)
}

on_failure(is_turn_angle) <- function(call,env) {
  paste0(deparse(call$x)," is not a valid turn angle. See documentation for details")
}

###################################################
### hydro scale
###################################################

is_hydro_scale <-function(x) {
  assert_that(is.list(x), length(x)==20)
  not_empty(x$A)==T & not_empty(x$R)==T & not_empty(x$N)==T & not_empty(x$D)==T & not_empty(x$C)==T & not_empty(x$E)==T & not_empty(x$Q)==T & not_empty(x$G)==T & not_empty(x$H)==T & not_empty(x$I)==T & not_empty(x$L)==T & not_empty(x$K)==T & not_empty(x$M)==T & not_empty(x$F)==T & not_empty(x$P)==T & not_empty(x$S)==T & not_empty(x$T)==T & not_empty(x$W)==T & not_empty(x$Y)==T & not_empty(x$V)==T
}

on_failure(is_hydro_scale) <- function(call,env) {
  paste0(deparse(call$x)," is not a proper hydrophobicity scale. See documentation for details")
}
