################################################################################
### Functions for summarizing reports
##
## Created on: 2016-08-19
## Author: Kazuki Yoshida
################################################################################


##' Sequentially replace strings
##'
##' Wrapper for multiple chained \code{gsub} replacement. See the explanation of \code{pat_rep_vec} for details
##'
##' @param x vector of strings. Transformed to character if of other type.
##' @param pat_rep_vec vector with even number of elements. Each odd number pattern is replaced with the following even number string. The replacement proceed sequentially, so the order of these pattern-replacement pairs matter.
##'
##' @return a vector with replaced strings
##'
##' @author Kazuki Yoshida
##'
##' @export
SeqGsub <- function(x, pat_rep_vec) {
    ## There must be even number of elements
    assertthat::assert_that((length(pat_rep_vec) %% 2) == 0)
    ## Characters are easier to work on.
    x <- as.character(x)
    ##
    ind <- seq_along(pat_rep_vec)
    ## Odd elements are patterns
    patterns     <- pat_rep_vec[(ind %% 2) == 1]
    ## Even elements are replacements
    replacements <- pat_rep_vec[(ind %% 2) == 0]

    ## Sequentially replace
    for (i in seq_along(patterns)) {
        x <- gsub(pattern = patterns[i],
                  replacement = replacements[i],
                  x = x)
    }

    x
}
