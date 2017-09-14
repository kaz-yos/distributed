################################################################################
### Miscellaneous undocumented functions
##
## Created on: 2016-05-17
## Author: Kazuki Yoshida
################################################################################

### Simple math functions
expit <- function(x) {
    e_x <- exp(x)
    e_x / (1 + e_x)
}
##
logit <- function(p) {
    log(p / (1 - p))
}


### Extract model coef and corresponding var
CoefVar <- function(model, variable) {
    assert_that(!is.null(variable))

    if (is.error(model)) {
        ## If the model is broken return numeric c(NA,NA)
        c(coef = as.numeric(NA),
          var  = as.numeric(NA))
    } else {
        c(coef = as.numeric(coef(model)[variable]),
          var  = as.numeric(vcov(model)[variable, variable]))
    }
}


### Detect errors given by try()
## http://adv-r.had.co.nz/Exceptions-Debugging.html#condition-handling
##' Tests for objects of type \code{try-error}
##'
##' Tests for objects of type \code{try-error}
##'
##' @param x return object from a function wrapped by \code{try}.
##'
##' @return logical indicating whether the return object is an error object given by \code{try} or normal return object.
##'
##' @author Kazuki Yoshida
##'
##' @export
is.error <- function(x) {
    inherits(x, "try-error")
}


### var() that returns 0 if variance is NA
safe_var <- function(x) {
    ifelse(length(x) == 0,
           ## No observation, cannot vary, thus, defined 0.
           0,
    ifelse(length(x) == 1 & !is.na(x[1]),
           ## Only one observation, cannot vary, thus, defined 0.
           0,
           var(x)))
}


### Generate two point distribution sample with specified sample mean and variance
## a is the departure from mean
TwoPointSample <- function(n, mean, var) {
    assert_that(length(n) == 1)
    assert_that(length(mean) == 1)
    assert_that(length(var) == 1)
    ## Round for safety
    n <- round(n)

    if (n == 0) {
        ## No element
        return(as.numeric())
    } else if (n == 1) {
        ## Only one element just replicate mean
        return(mean)
    } else if (n %% 2 == 0) {
        ## Even
        ## Departure from origin
        a = sqrt((n-1)/n * var)
        out <- rep(c(a,-a), each = n / 2)
        ## Shift
        return(out + mean)
    } else {
        ## Odd
        ## Departure from origin
        a = sqrt(var)
        ## Center the last element
        out <- c(rep(c(a,-a), each = (n - 1) / 2), 0)
        ## Shift
        return(out + mean)
    }
}



### Generate a skewed two point distribution sample with specified sample mean and variance
## Give sample size n, sample mean m, and sample variance v
## Values are 0 < p1 < p0 = m < p2
## n1 people take on p1
## n2 people take on p2
## n1 >= n2 to make the mean closer to p1.
## n0 people (reminder people) are placed at p0 = m
## distances from the mean:
## d1 = m - p1;
## d2 = p2 - m
## d2 = r * d1

##' Create a sample with skewed two or three point distribution
##'
##' Create a sample with skewed two or three point distribution with a pre-specified non-negative mean and variance. The purpose is to create a two-point or three-point distribution that take on only positive values having a specified sample mean (m) and sample variance (v). By definition of the sample mean, we need to have at least one value below the mean (p1) and at least one value above the mean (p2). As this function is intended to reproduce a distribution for weights, the lower value cannot be negative (p1 > 0). This is achieved by placing more subjects at p1 than p2 such that p1 is closer to the mean (distance d1 = m - p1) than p2 (distance d2 = p2 - m). The algorithm tries to put progressively more people at p1 until it achieves p1 > 0. See also the details.
##'
##' @param n Required number of individuals in the sample. If 0, an empty row data frame is returned.
##' @param m Required positive sample mean. As this is used to reconstruct a weight distribution, it has to be positive. NaN created by 0/0 is allowed if n = 0.
##' @param v Required non-negative sample variance. As the weights can be the same for all individuals, zero is permitted.
##' @param max_r Maximum for the ratio of the individuals below the mean to the individuals above the mean. The algorithm starts at \code{r = 1} and proceed until the desired results are obtained, \code{max_r} is reached, or an error is encountered.
##'
##' @details For a given sample size n, sample mean m, and sample variance v, we need to find values p1 and p2 satisfying 0 < p1 < m < p2 that result in a sample with the desired n, m, and v. The desired squared sum is (n-1)v. The algorithm proceed as follows. Initialize r = 1. Obtain maximum K such that K(r+1) <= n. Let n0 = n - K(r+1). These n0 subjects are reminder subjects to be placed at the sample mean m. Let n1 = Kr. These n1 subjects are to be placed at p1. Let n2 = K. These n2 subjects are to be placed at p2. Notice that n1 / n2 = Kr / K = r. Because of this group size ratio, to maintain the mean m, the distances of p1 and p2 from the mean must have the inverse ratio. That is d1 / d2 = 1 / r and d2 = rd1. Note when r = 1, we have the same number of subjects at both p1 and p2, thus, d1 = d2 and the mean m is at the mid-point between p1 and p2. As r increases, the mean m is progressively closer to p1. Because the mean is at m, the total squared sum contribution of these (n1 + n2) subjects is n1 d1^2 + n2 d2^2 = n1 d1^2 + n2 (rd1)^2. This expression with an unknown d1 is equated with the desired squared sum (n-1)v and solved for d1. This results in d1^2 = (n-1)v / (n1 + r^2 n2). The positive root is the solution we want for d1 (it is a function of r). Importantly, we can ignore the n0 subjects placed at the mean because they all have 0 squared sum contributions. This d1 may result in negative p1 if d1 > m. If this is the case, r is incremented by 1 and algorithm is rerun. This will continue until finding the smallest r that gives a desired all positive two-/three-point sample distribution, r_max is reached, or sample size n is exhausted before finding r (sample size is too small to achieve the very small mean m and high variance v).
##'
##' @return data frame with a \code{W} column containing the values and a \code{count} column containing the number of subjects having value \code{W}.
##'
##' @author Kazuki Yoshida
##'
##' @export
SkewedTwoPointSample <- function(n, m, v, max_r = 10^3) {

    ## Negative sample size is broken
    assert_that(n >= 0)

    ## Handle easy cases first
    if (n <= 0) {
        ## Empty data frame
        return(data.frame(W = as.numeric(),
                          count = as.numeric()))
    } else if (n == 1 | v == 0) {
        ## No choice if only one observation or no variance
        assert_that(m > 0)
        return(data.frame(W = m,
                          count = n))
    }
    ## Otherwise proceed to the iterative algorithm
    ## Test here to allow for n = 0, m = NaN = 0/0 case
    assert_that(m >  0)
    assert_that(v >= 0)

    ## Target squared sum
    ss <- (n - 1) * v

    ## Initial group ratio
    r <- 1

    ## Loop until max_r is reached or not enough n
    while (r <= max_r & r < n) {
        ## Reminder people are placed at the mean
        n0 <- n %% (r + 1)
        ## Number of r:1 units
        K <- n %/% (r + 1)
        ## Number of people placed at p1
        n1 <- r * K
        ## Number of people placed at p2
        n2 <-     K
        ## Distance d1 = m - p1
        d1 <- sqrt(ss / (n1 + n2 * r^2))
        ## Lower point
        p1 <- m - d1
        ## Proceed if p1 > 0
        if (p1 > 0) {
            ## Upper point
            p2 <- m + (r * d1)
            ## Points and numbers
            if (n0 > 0) {
                ## If middle piece exists
                ps <- c(p1, m,  p2)
                ns <- c(n1, n0, n2)
            } else {
                ## Otherwise, just two elements
                ps <- c(p1, p2)
                ns <- c(n1, n2)
            }
            ## Return a data frame
            return(data.frame(W = ps,
                              count = ns))
        } else {
            ## If p1 <= 0, increment r and try again
            r <- r + 1
        }
    }

    ## If getting out of the above while loop, it's an error
    if (r > max_r) {
        stop("max_r reached before achieving the required mean and variance")
    } else if (r >= n) {
        stop("Not enough n to achieve the required mean and variance")
    } else {
        stop("while loop was existed for an unknown reason")
    }
}


### Generate a two point distribution (1 point above mean all others between 0 and mean)
##' Generate a two point distribution sample with only one point above mean.
##'
##' Generate a two point distribution sample with only one point above mean with the desired mean and variance.
##'
##' @param n Required number of individuals in the sample. If 0, an empty row data frame is returned.
##' @param m Required positive sample mean. As this is used to reconstruct a weight distribution, it has to be positive. NaN created by 0/0 is allowed if n = 0.
##' @param v Required non-negative sample variance. As the weights can be the same for all individuals, zero is permitted.
##' @return data frame with a \code{W} column containing the values and a \code{count} column containing the number of subjects having value \code{W}.
##'
##' @author Kazuki Yoshida
ExtremeTwoPointSample <- function(n, m, v) {
    ## Negative sample size is broken
    assert_that(n >= 0)

    ## Handle easy cases first
    if (n <= 0) {
        ## Empty data frame
        return(data.frame(W = as.numeric(),
                          count = as.numeric()))
    } else if (n == 1 | v == 0) {
        ## No choice if only one observation or no variance
        assert_that(m > 0)
        return(data.frame(W = m,
                          count = n))
    }
    ## Otherwise proceed to the iterative algorithm
    ## Test here to allow for n = 0, m = NaN = 0/0 case
    assert_that(m >  0)
    assert_that(v >= 0)

    ## Calculate the lower data point
    lower <- m - sqrt(v / n)

    ## Error if negative
    if (lower <= 0) {
        stop("Not enough n to achieve the required mean and variance")
    }

    upper <- n * m - (n - 1) * lower

    data.frame(W = c(lower, upper),
               count = c(n - 1, 1))
}


### Generate a desired NA data frame
##' Generate a desired NA data frame
##'
##' Generate a desired NA data frame having specified number of rows, column names, and column classes.
##'
##' @param labels character vector for column names
##' @param as_funcs vector transforming functions such as \code{as.integer}, etc
##' @param nrow scalar value to specify the number of rows
##'
##' @return a data frame filled with |code{NA}'s and \code{nrow} rows with column names as specified in \code{labels} and column classes as specified in \code{as_funcs}.
##'
##' @author Kazuki Yoshida
CreateNaDf <- function(labels, as_funcs, nrow) {
    assert_that(is.numeric(nrow))
    assert_that(length(nrow) == 1)
    assert_that(length(labels) == length(as_funcs))
    assert_that(is.character(labels))
    assert_that(all(sapply(as_funcs, is.function)))

    out <- lapply(as_funcs, function(f) {
        f(rep(NA, nrow))
    }) %>% as.data.frame

    names(out) <- labels
    out
}


### Validate data frame for variability in exposure and outcome
##' Validate data frame for variability in exposure and outcome
##'
##' Validate a given data frame for variability in the exposure variable and outcome variable. The exposure variable specified in \code{expo_var} should have two values to make regression meaningful. Only one of \code{out_var} and \code{event_var} is required although both can be validated at the same time. Because the error raised can stop the simulation, it should be used within \code{try} function.
##'
##' @param df data frame to be used in a regression analysis
##' @param expo_var string for the exposure variable name
##' @param bin_var string for the binary outcome variable name
##' @param event_var string for the survival event indicator variable name (1 for event)
##'
##' @return the validated data frame if it is valid. Otherwise, an error is raised to stop estimation on invalid data.
##'
##' @author Kazuki Yoshida
ValidateDf <- function(df, expo_var, bin_var = NULL, event_var = NULL) {
    assert_that(is.character(expo_var) & length(expo_var) == 1)

    ## Check for degenerate exposure variable
    if (length(unique(df[,expo_var])) <= 1) {
        stop("Degenerate exposure variable. Effect estimation not possible.")
    }

    ## Check for valid binary outcome data
    if (!is.null(bin_var)) {
        assert_that(is.character(bin_var) & length(bin_var) == 1)
        tab <- table(df[,expo_var], df[,bin_var])
        ## No zero cells should exist in exposure-outcome cross table for valid estimation
        if (!all(tab > 0)) {
            stop("Zero cell(s) in exposure binary outcome cross table.")
        }
    }

    ## Check for valid event data for survival outcome data
    if (!is.null(event_var)) {
        assert_that(is.character(event_var) & length(event_var) == 1)
        tab <- tapply(df[,event_var], df[,expo_var], FUN = mean)
        ## Both exposure groups must have some events
        if (!all(tab > 0)) {
            stop("No events in either or both exposure groups in survival data.")
        }
    }

    ## If none of the errors have been raised, return as is
    df
}


##' Enhanced try that can treat some warnings as errors
##'
##' Behaves like \code{try} to return \code{try-error} object on an error. If a warning is encountered it is checked against the specified regular expression \code{regexp}. If the warning message matches the regular expression, it is handled as an error and return \code{try-error} object.
##'
##' @param expr expression to be evaluated that may give an warning or error.
##' @param regexp regular expression to specify warning to be handled as an error.
##' @param silent whether to silence error/warning messages.
##'
##' @return \code{try-error} object if an error or a specified warning is encountered. Otherwise, the result of evaluation returns.
##'
##' @author Kazuki Yoshida
##'
##' @export
try_w <- function(expr, regexp = NULL, silent = FALSE) {
    ## http://adv-r.had.co.nz/Exceptions-Debugging.html#condition-handling
    if (is.null(regexp)) {
        ## Just use try if no regexp is given
        try(expr, silent = silent)

    } else {
        ## If regexp is given, warnings matching regexp are considered fatal errors.
        tryCatch(expr,
                 ## Handling of an error
                 error = function(cond) {
                     ## If not silenet give message
                     if (!silent) {message(conditionMessage(cond))}

                     ## Return try-error object
                     invisible(structure(conditionMessage(cond),
                                         class = "try-error"))

                 },
                 ## Handling of a warning
                 warning = function(cond) {
                     ## If not silenet give message
                     if (!silent) {message(conditionMessage(cond))}

                     ## If it is a specific warning matching regexp
                     if (grepl(regexp, conditionMessage(cond))) {
                         ## Return try-error object
                         invisible(structure(conditionMessage(cond),
                                             class = "try-error"))
                     } else {
                         ## Otherwise, return the execution results
                         suppressWarnings(expr)
                     }
                 })
    }
}
