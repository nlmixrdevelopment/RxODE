\name{genShinyApp.template}
\alias{genShinyApp.template}
\alias{write.template.server}
\alias{write.template.ui}
\title{
   Generate an example (template) of a dosing regimen shiny app
}
\description{
   Create a complete shiny application for exploring dosing regimens
   given a (hardcoded) PK/PD model.
}
\usage{
genShinyApp.template(appDir = "shinyExample", verbose = TRUE)
write.template.server(appDir)
write.template.ui(appDir)
}
\arguments{
  \item{appDir}{
  a string with a directory where to store the shiny app, by default
  is \code{"shinyExample"}.
  The directory \code{appDir} will be created if it does not exist.
}
  \item{verbose}{
  logical specifying whether to write messages as the shiny app
  is generated. Defaults to \code{TRUE}.
}
}
\details{
A PK/PD model is defined using \code{\link{RxODE}}, and
a set of parameters and initial values are defined.  Then
the appropriate R scripts for the shiny's user interface \code{ui.R} 
and the server logic \code{server.R} are created in the 
directory \code{appDir}.

The function evaluates the following PK/PD model:
\preformatted{
   C2 = centr/V2;
   C3 = peri/V3;
   d/dt(depot) =-KA*depot;
   d/dt(centr) = KA*depot - CL*C2 - Q*C2 + Q*C3;
   d/dt(peri)  =                    Q*C2 - Q*C3;
   d/dt(eff)  = Kin - Kout*(1-C2/(EC50+C2))*eff; 
}

To launch the shiny app, simply issue the \code{runApp(appDir)}
R command.
}
\value{
   None, these functions are used for their side effects.
}
\note{
   These functions create a simple, but working example of a dosing
   regimen simulation web application. Users may want to modify the code 
   to experiment creating shiny applications for their specific 
   \code{RxODE} models.
}

\seealso{
  \code{\link{RxODE}}, \code{\link{eventTable}}, and the
  package \pkg{shiny} (\url{shiny.rstudio.com}.
}
\examples{
\dontrun{
   # create the shiny app example (template)
   genShinyApp.template(appDir = "myapp")

   # run the shiny app
   runApp("myapp")
}
}
\keyword{simulation}
\keyword{nonlinear}
\concept{PK/PD}
\concept{pharmacometrics}
