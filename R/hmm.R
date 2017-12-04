library(assertthat)
library(testthat)

source('R/inputs.R')

#######################################################
### Examples of peptide sequences:
### H1 - influenza strain H1 fusion peptide
### gp41 - HIV fusion peptide

H1<-"GLFGAIAGFIEGGWTGMIDGWYG"
gp41 <- "AVGIGALFLGFLGAAGSTMGAAS"


#######################################################
### Wimley-White octanol scale

current_scale <- list('A'=0.5, 'R'=1.81, 'N'=0.85, 'D'=3.64, 'C'=-0.02, 'E'=3.63, 'Q'=0.77, 'G'=1.15, 'H'=2.33, 'I'=-1.12, 'L'=-1.25, 'K'=2.8, 'M'=-0.67, 'F'=-1.71, 'P'=0.14, 'S'=0.46, 'T'=0.25, 'W'=-2.09, 'Y'=-0.71, 'V'=-0.46)




#######################################################
### H_average

#' Calculates average hydrophobicity
#'
#' Calculates average hydrophobicity of a given peptide sequence based on a hydrophobicity scale
#' @param hydro_scale Hydrophobicity scale in a form of list (Wimley_White octanol scale used as default (current_scale))
#' @param prot_sequence Protein (peptide) sequence given as string (FASTA format with no description line, just the sequence)
#' @return Sequence-based average hydrophobicity
#' @export


H_average <- function(hydro_scale, prot_sequence){

  hydro_aa = NULL

  for (i in (1:nchar(prot_sequence))){
    hydro_aa[i] <- as.numeric(current_scale[substr(prot_sequence,i,i)])
  }
  return(mean(hydro_aa))
}


H_moment_norm <- function(hydro_scale, prot_sequence, turn_angle){

  assert_that(is_protein(prot_sequence))
  assert_that(is_turn_angle(turn_angle))
  assert_that(is_hydro_scale(hydro_scale))

  sin_part = NULL
  cos_part = NULL

  for (i in (1:nchar(prot_sequence))){
    hydr_aa <- as.numeric(current_scale[substr(prot_sequence,i,i)])
    sin_part[i] <- hydr_aa * sin((i)*turn_angle)
    cos_part[i] <- hydr_aa * cos((i)*turn_angle)
  }
  mju_H <- (sum(sin_part)^2 + sum(cos_part)^2)^0.5
  mju_H_norm <- mju_H/nchar(prot_sequence)

  return(mju_H_norm)
}

#######################################################
### H_moment_norm

#' Calculates normalized hydrophobicity moment
#'
#' Calculates normalized hydrophobicity moment for a given peptide sequence for a given hydrophobicity scale and helix trun angle.
#' @param hydro_scale Hydrophobicity scale in a form of list (Wimley_White octanol scale used as default (current_scale))
#' @param prot_sequence Protein (peptide) sequence given as string (FASTA format with no description line, just the sequence)
#' @param turn_angle Helix turn angle in radians
#' @return Single value of hydrophobic moment for a given sequence
#' @export


H_moment_norm <- function(hydro_scale, prot_sequence, turn_angle){

  assert_that(is_protein(prot_sequence))
  assert_that(is_turn_angle(turn_angle))
  assert_that(is_hydro_scale(hydro_scale))

  sin_part = NULL
  cos_part = NULL

  for (i in (1:nchar(prot_sequence))){
    hydr_aa <- as.numeric(current_scale[substr(prot_sequence,i,i)])
    sin_part[i] <- hydr_aa * sin((i)*turn_angle)
    cos_part[i] <- hydr_aa * cos((i)*turn_angle)
  }
  mju_H <- (sum(sin_part)^2 + sum(cos_part)^2)^0.5
  mju_H_norm <- mju_H/nchar(prot_sequence)

  return(mju_H_norm)
}


#######################################################
### H_moment_window

#' Calculates normalized hydrophobicity moment of a peptide for a given window size
#'
#' Calculates normalized hydrophobicity moment for a given peptide sequence for a given hydrophobicity scale, helix trun angle and window size. For the given peptide sequence, hydrophobic moment is calculated for a window (of 3 or 5 or 7 or 9 amino acid) assuming helical structure with a given turn angle.
#' @param hydro_scale Hydrophobicity scale in a form of list (Wimley_White octanol scale used as default (current_scale))
#' @param prot_sequence Protein (peptide) sequence given as string (FASTA format with no description line, just the sequence)
#' @param turn_angle Helix turn angle in radians
#' @param window Size of window (3, 5, 7, 9)
#' @return A vector of hydrophobic moment values for a given sequence (has a length of: (number of amino acid - window size)+1 )
#' @export

H_moment_window <- function(hydro_scale, prot_sequence, window, turn_angle){

  assert_that(is_window(window))

  H_value = NULL

  for (i in (1:(nchar(prot_sequence)-window+1))){
    H_value[i] <- H_moment_norm(hydro_scale,substr(prot_sequence,i,window+i-1),turn_angle)
  }
  return(H_value)
}


#######################################################
### H_moment_map

#' Draws a map of hydrophobic moment of a peptide for a given window size
#'
#' Calculates and draws hydrophobic moment map for a given peptide sequence for a given hydrophobicity scale, and window size. The results are presented graphically as a map for the following coordinates: amino acid number (middle of window) and helical turn angle (in degrees). More about hydrophobic moment maps can be found in the \href{http://onlinelibrary.wiley.com/doi/10.1016/j.febslet.2013.06.054/full}{article}.
#' @param hydro_scale Hydrophobicity scale in a form of list (Wimley_White octanol scale used as default (current_scale))
#' @param prot_sequence Protein (peptide) sequence given as string (FASTA format with no description line, just the sequence). Must be shorter than 30 amino acids.
#' @param window Size of window (3, 5, 7, 9)
#' @return A vector of hydrophobic moment values for a given sequence (has a length of: (number of amino acid - window size)+1 )
#' @export

H_moment_map <- function(hydro_scale, prot_sequence, window){

  H_mom_value <- H_moment_window(current_scale, prot_sequence, window, 0)
  half_window <- window %/% 2

  for (i in seq(10,180,10)){
    H_mom_value_curr <- H_moment_window(current_scale, prot_sequence, window, turn_angle=i*pi/180)
    H_mom_value <- rbind(H_mom_value,H_mom_value_curr)
  }
  row.names(H_mom_value) <- seq(0,180,10)
  colnames(H_mom_value) <- seq(half_window + 1, nchar(prot_sequence) - half_window, 1)
  heatmap(H_mom_value, Colv=NA, Rowv=NA, xlab="aa", ylab="turn angle (deg)")

  return(H_mom_value)
}

#######################################################
### H_moment_cross

#' Draws a cross-section of a hydrophobicity moment map for a given turn angle and window size
#'
#' Calculates and draws a cross-section of a hydrophobicity moment map for a given turn angle and a window size (3,5,7,or 9).
#' @param hydro_scale Hydrophobicity scale in a form of list (Wimley_White octanol scale used as default (current_scale))
#' @param prot_sequence Protein (peptide) sequence given as string (FASTA format with no description line, just the sequence)
#' @param turn_angle Helix turn angle in radians
#' @param window Size of window (3, 5, 7, 9)
#' @return A vector of hydrophobic moment values for a given sequence (has a length of: (number of amino acid - window size)+1 )
#' @export

H_moment_cross <- function(hydro_scale, prot_sequence, window, turn_angle){

  half_window <- window %/% 2

  H_value = NULL
  for (i in (1:(nchar(prot_sequence)-window+1))){
    H_value[i] <- H_moment_norm(hydro_scale,substr(prot_sequence,i,window+i-1),turn_angle)
  }
  plot (seq(half_window + 1, nchar(prot_sequence) - half_window, 1), H_value, type="b", xlab="aa", ylab="mju_H")

  return(H_value)
}
