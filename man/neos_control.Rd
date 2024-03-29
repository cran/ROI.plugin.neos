% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plugin.R
\name{neos_control}
\alias{neos_control}
\title{Neos Control Variables}
\usage{
neos_control(
  method = "auto",
  wait = TRUE,
  email = "",
  password = "",
  user = "rneos",
  dry_run = FALSE,
  options = "",
  parameters = "",
  gdx = "",
  restart = "",
  wantgdx = "",
  wantlst = "",
  wantlog = "",
  comments = ""
)
}
\arguments{
\item{method}{a chracter string giving the name of the solver to be selected
on the \code{NEOS} server.}

\item{wait}{a logical indicating whether the R interpreter should wait for the
command to finish, or run it asynchronously. If \code{TRUE}
\pkg{ROI} returns an object of class \code{"neos_job"}.}

\item{email}{a character string giving the email address.}

\item{password}{a character string giving the account password.}

\item{user}{a character string giving the username.}

\item{dry_run}{a logical if \code{TRUE} \pkg{ROI} returns the solver call.}

\item{options}{a character string (default is \code{""}) passed to \code{options}
tag of the \code{GAMS} solver template.}

\item{parameters}{a character string (default is \code{""}) passed to \code{parameters}
tag of the \code{GAMS} solver template.}

\item{gdx}{a character string (default is \code{""}) passed to \code{gdx}
tag of the \code{GAMS} solver template.}

\item{restart}{a character string (default is \code{""}) passed to \code{restart}
tag of the \code{GAMS} solver template.}

\item{wantgdx}{a character string (default is \code{""}) passed to \code{wantgdx}
tag of the \code{GAMS} solver template.}

\item{wantlst}{a character string (default is \code{""}) passed to \code{wantlst}
tag of the \code{GAMS} solver template.}

\item{wantlog}{a character string (default is \code{""}) passed to \code{wantlog}
tag of the \code{GAMS} solver template.}

\item{comments}{a character string (default is \code{""}) passed to \code{comments}
tag of the \code{GAMS} solver template.}
}
\description{
The control variables for \code{ROI.plugin.neos}.
}
