################################################################################
### Test miscellaneous functions
##
## Created on: 2016-05-05
## Author: Kazuki Yoshida
################################################################################

library(testthat)

### Context (1 for each file)
context("### Test 00.Miscellaneous.R")


test_that("is.error works as expected", {

    ## A list with elements that may give failure
    lst <- list(1, 2, "hello")
    ## Third element fails
    out <- lapply(lst, function(x) {try(1 + x, silent = TRUE)})
    ## return objects are either correct return value or try-error
    expect_equal(sapply(out, is.numeric),
                 c(TRUE, TRUE, FALSE))
    expect_equal(sapply(out, is.error),
                 c(FALSE, FALSE, TRUE))

})


test_that("broken model is handled correctly", {

    data1 <- data.frame(x = c(0),
                        y = c(-1))

    out <- try(glm(formula = y ~ x,
                   family  = binomial(link = "logit"),
                   data    = data1))

    ## Expectations
    expect_equal(CoefVar(out, "x"),
                 c(coef = as.numeric(NA),
                   var  = as.numeric(NA)))
})


test_that("variance calculation is safe", {

    ## For no element defined as 0
    expect_equal(safe_var(x = as.numeric()),
                 0)
    ## For single element defined as 1 unless it is NA
    expect_equal(safe_var(x = NA),
                 as.numeric(NA))
    expect_equal(safe_var(x = 3),
                 0)
    expect_equal(safe_var(x = 7),
                 0)
    ## Otherwise same as var
    expect_equal(safe_var(x = c(1,2,3)),
                 var(c(1,2,3)))
})


test_that("two-point sample is generated correctly", {

    ## Set values
    m <- 3.3
    v <- 7.2
    sample_even <- TwoPointSample(n = 10, mean = m, var = v)
    sample_odd  <- TwoPointSample(n = 11, mean = m, var = v)
    ## Expectations
    expect_equal(mean(sample_even), m)
    expect_equal(mean(sample_odd), m)
    expect_equal(var(sample_even), v)
    expect_equal(var(sample_odd), v)

    ## Set values
    m <- 0.1
    v <- 0.02
    sample_even <- TwoPointSample(n = 2, mean = m, var = v)
    sample_odd  <- TwoPointSample(n = 3, mean = m, var = v)
    ## Expectations
    expect_equal(mean(sample_even), m)
    expect_equal(mean(sample_odd), m)
    expect_equal(var(sample_even), v)
    expect_equal(var(sample_odd), v)

    ## Set values
    m <- 0.1
    v <- 0.02
    sample_even <- TwoPointSample(n = 20, mean = m, var = v)
    sample_odd  <- TwoPointSample(n = 33, mean = m, var = v)
    ## Expectations
    expect_equal(mean(sample_even), m)
    expect_equal(mean(sample_odd), m)
    expect_equal(var(sample_even), v)
    expect_equal(var(sample_odd), v)

})


test_that("skewed three-point sample is generated correctly", {

    ## no data point -> empty data frame
    expect_equal(SkewedTwoPointSample(n = 0,
                                      m = 0.1,
                                      v = 0),
                 data.frame(W = as.numeric(), count = as.numeric()))
    ## Mean 0 is unreasonable
    expect_error(SkewedTwoPointSample(n = 1,
                                      m = 0,
                                      v = 0),
                 "m not greater than 0")
    ## one data point no choice, but mean
    expect_equal(SkewedTwoPointSample(n = 1,
                                      m = 0.33,
                                      v = 0),
                 data.frame(W = 0.33, count = 1))
    ## one data point; ignores variance
    expect_equal(SkewedTwoPointSample(n = 1,
                                      m = 0.33,
                                      v = 0.1),
                 data.frame(W = 0.33, count = 1))
    ## If variance is zero, put everyone at m.
    expect_equal(SkewedTwoPointSample(n = 100,
                                      m = 0.33,
                                      v = 0),
                 data.frame(W = 0.33, count = 100))
    ## Correct mean and variance
    sample1 <- SkewedTwoPointSample(n = 3,
                                    m = 0.33,
                                    v = 0.1)
    print(sample1)
    expect_equal(mean(rep(sample1$W, sample1$count)), 0.33)
    expect_equal(var(rep(sample1$W, sample1$count)), 0.1)
    ##
    sample2 <- SkewedTwoPointSample(n = 100,
                                    m = 0.33,
                                    v = 0.5)
    print(sample2)
    expect_equal(mean(rep(sample2$W, sample2$count)), 0.33)
    expect_equal(var(rep(sample2$W, sample2$count)), 0.5)
    ##
    sample3 <- SkewedTwoPointSample(n = 100,
                                    m = 0.25,
                                    v = 0.5)
    print(sample3)
    expect_equal(mean(rep(sample3$W, sample3$count)), 0.25)
    expect_equal(var(rep(sample3$W, sample3$count)), 0.5)
    ##
    sample4 <- SkewedTwoPointSample(n = 100,
                                    m = 0.1,
                                    v = 0.5,
                                    max_r = 100)
    print(sample4)
    expect_equal(mean(rep(sample4$W, sample4$count)), 0.1)
    expect_equal(var(rep(sample4$W, sample4$count)), 0.5)
    ## Minimum achievable with these n and v
    sample5 <- SkewedTwoPointSample(n = 100,
                                    m = 0.0708,
                                    v = 0.5,
                                    max_r = 100)
    print(sample5)
    expect_equal(mean(rep(sample5$W, sample5$count)), 0.0708)
    expect_equal(var(rep(sample5$W, sample5$count)), 0.5)
    ##
    expect_error(SkewedTwoPointSample(n = 100,
                                      m = 0.05,
                                      v = 0.5,
                                      max_r = 100),
                 "Not enough n")

})


test_that("extreme two-point sample is generated correctly", {

    ## no data point -> empty data frame
    expect_equal(ExtremeTwoPointSample(n = 0,
                                       m = 0.1,
                                       v = 0),
                 data.frame(W = as.numeric(), count = as.numeric()))
    ## Mean 0 is unreasonable
    expect_error(ExtremeTwoPointSample(n = 1,
                                       m = 0,
                                       v = 0),
                 "m not greater than 0")
    ## one data point no choice, but mean
    expect_equal(ExtremeTwoPointSample(n = 1,
                                       m = 0.33,
                                       v = 0),
                 data.frame(W = 0.33, count = 1))
    ## one data point; ignores variance
    expect_equal(ExtremeTwoPointSample(n = 1,
                                       m = 0.33,
                                       v = 0.1),
                 data.frame(W = 0.33, count = 1))
    ## If variance is zero, put everyone at m.
    expect_equal(ExtremeTwoPointSample(n = 100,
                                       m = 0.33,
                                       v = 0),
                 data.frame(W = 0.33, count = 100))
    ## Correct mean and variance
    sample1 <- ExtremeTwoPointSample(n = 3,
                                     m = 0.33,
                                     v = 0.1)
    print(sample1)
    expect_equal(nrow(sample1), 2)
    expect_equal(mean(rep(sample1$W, sample1$count)), 0.33)
    expect_equal(var(rep(sample1$W, sample1$count)), 0.1)
    ##
    sample2 <- ExtremeTwoPointSample(n = 100,
                                     m = 0.33,
                                     v = 0.5)
    print(sample2)
    expect_equal(nrow(sample2), 2)
    expect_equal(mean(rep(sample2$W, sample2$count)), 0.33)
    expect_equal(var(rep(sample2$W, sample2$count)), 0.5)
    ##
    sample3 <- ExtremeTwoPointSample(n = 100,
                                     m = 0.25,
                                     v = 0.5)
    print(sample3)
    expect_equal(nrow(sample3), 2)
    expect_equal(mean(rep(sample3$W, sample3$count)), 0.25)
    expect_equal(var(rep(sample3$W, sample3$count)), 0.5)
    ##
    sample4 <- ExtremeTwoPointSample(n = 100,
                                     m = 0.1,
                                     v = 0.5)
    print(sample4)
    expect_equal(nrow(sample4), 2)
    expect_equal(mean(rep(sample4$W, sample4$count)), 0.1)
    expect_equal(var(rep(sample4$W, sample4$count)), 0.5)
    ## Minimum achievable with these n and v
    sample5 <- ExtremeTwoPointSample(n = 100,
                                     m = 0.0708,
                                     v = 0.5)
    print(sample5)
    expect_equal(nrow(sample5), 2)
    expect_equal(mean(rep(sample5$W, sample5$count)), 0.0708)
    expect_equal(var(rep(sample5$W, sample5$count)), 0.5)
    ##
    expect_error(ExtremeTwoPointSample(n = 100,
                                       m = 0.05,
                                       v = 0.5),
                 "Not enough n")

})


test_that("NA data frames are constructed as desired", {

    ## Single row all integers
    res1 <- CreateNaDf(labels = c("A", "events", "denom"),
                       as_funcs = c(as.integer, as.integer, as.integer),
                       nrow = 1)
    expect_equal(res1,
                 data.frame(A      = as.integer(NA),
                            events = as.integer(NA),
                            denom  = as.integer(NA)))

    ## More complicated
    res2 <- CreateNaDf(labels = c("A", "events", "denom"),
                       as_funcs = c(as.double, as.integer, as.integer),
                       nrow = 2)
    expect_equal(res2,
                 data.frame(A      = as.double(rep(NA, 2)),
                            events = as.integer(rep(NA, 2)),
                            denom  = as.integer(rep(NA, 2))))

})


test_that("Exposure/outcome validation works", {

    df_valid <- data.frame(A = c(0,0,1,1),
                           Y = c(1,0,1,0),
                           event = c(1,1,1,1))

    df_degenerate_A <- data.frame(A = c(1,1),
                                  Y = c(1,0))

    df_zero_cell <- data.frame(A = c(0,0,1,1),
                               Y = c(1,0,0,0))

    df_no_events <- data.frame(A = c(0,0,1,1),
                               event = c(1,0,0,0))

    ## This is ok
    df_all_events <- data.frame(A = c(0,0,1,1),
                                event = c(1,0,1,1))


    ## Expectations
    expect_equal(ValidateDf(df = df_valid, expo_var = "A", bin_var = "Y", event_var = "event"),
                 df_valid)

    expect_error(ValidateDf(df = df_degenerate_A, expo_var = "A", bin_var = "Y", event_var = "event"),
                 "Degenerate exposure variable")

    expect_error(ValidateDf(df = df_zero_cell, expo_var = "A", bin_var = "Y"),
                 "Zero cell\\(s\\) in exposure binary outcome cross table")

    expect_error(ValidateDf(df = df_no_events, expo_var = "A", event_var = "event"),
                 "No events in either or both exposure groups in survival data")

    expect_equal(ValidateDf(df = df_all_events, expo_var = "A", event_var = "event"),
                 df_all_events)

})


test_that("warning handler works as expected", {

    ## Without regexp, just use try()
    ## If fatal error, return try-error object
    expect_equal(try_w((1:10 + "1:3"), silent = TRUE),
                try((1:10 + "1:3"), silent = TRUE))
    ## If fatal error, return try-error object
    expect_equal(try_w(1:10),
                 try(1:10))

    ## With regexp
    ## If not specified warning just pass value
    expect_equal(try_w((1:10 + 1:3), "Crazy Warning", silent = TRUE),
                 c(2L, 4L, 6L, 5L, 7L, 9L, 8L, 10L, 12L, 11L))
    expect_equal(try_w((1:10 + 1:3), "Crazy Warning", silent = FALSE),
                 c(2L, 4L, 6L, 5L, 7L, 9L, 8L, 10L, 12L, 11L))
    ## If specified warning, return try-error object
    expect_true(is.error(try_w((1:10 + 1:3), regexp = "not a multiple of", silent = TRUE)))
    expect_true(is.error(try_w((1:10 + 1:3), regexp = "not a multiple of", silent = FALSE)))
    ## If fatal error, return try-error object
    expect_true(is.error(try_w((1:10 + "1:3"), "not a multiple of", silent = TRUE)))
    expect_true(is.error(try_w((1:10 + "1:3"), "not a multiple of", silent = FALSE)))

})
