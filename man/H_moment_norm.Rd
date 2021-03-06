% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hmm.R
\name{H_moment_norm}
\alias{H_moment_norm}
\title{Calculates normalized hydrophobicity moment}
\usage{
H_moment_norm(hydro_scale, prot_sequence, turn_angle)
}
\arguments{
\item{hydro_scale}{Hydrophobicity scale in a form of list (Wimley_White octanol scale used as default (current_scale))}

\item{prot_sequence}{Protein (peptide) sequence given as string (FASTA format with no description line, just the sequence)}

\item{turn_angle}{Helix turn angle in radians}
}
\value{
Single value of hydrophobic moment for a given sequence
}
\description{
Calculates normalized hydrophobicity moment for a given peptide sequence for a given hydrophobicity scale and helix trun angle.
}
