##' Return if the current document is a Rstudio notebook
##'
##' @return If we are in a notebook
##' @author  Max Gordon
##' @keywords internal
##' @examples
##'
##' rxIsNotebook() # Check to see if you are running in a notebook.
##'
##' @export
rxIsNotebook <- function() {
  if (!rstudioapi::isAvailable()) {
    return(FALSE)
  }

  ctxt <- rstudioapi::getActiveDocumentContext()
  if (grepl("\\.Rmd$", ctxt$path)) {
    return(TRUE)
  }

  # Look for html_notebook within the header if the file hasn't been saved
  contents <- ctxt$contents
  header <- grep("^---$", contents)
  if (length(header) == 2) {
    return(any(grepl("html_notebook$",
                     contents[min(header) : max(header)])))
  }

  return(FALSE)
}
