library(testthat)
library(wooldridge)  # Dataset source
library(AER)  # Required for ivreg
library(dplyr)

test_that("siv_reg runs without errors and produces expected outputs", {

  # Load dataset and preprocess
  df <- wooldridge::mroz
  data <- df[complete.cases(df), ]  # Remove missing values

  # Run the function with test parameters
  result <- siv_reg(data, "hours", "lwage",
                           c("educ", "age", "kidslt6", "kidsge6", "nwifeinc"), reps = 5)

  # Check if the function returns a list
  expect_type(result, "list")

  # Check if the required components exist
  expect_true("IV1" %in% names(result))
  expect_true("IV2" %in% names(result))
  expect_true("IV3" %in% names(result))
  expect_true("citable" %in% names(result))

  # Check if IV2, IV3, and IV4 are objects of class "ivreg" (instrumental variable regression)
  expect_s3_class(result$IV1, "ivreg")
  expect_s3_class(result$IV2, "ivreg")
  expect_s3_class(result$IV3, "ivreg")

  # Check if coefficient estimate of IV2 is within an expected range (example threshold)
  expect_true(abs(result$IV2$coefficients[2]) > 0)

  # Check if the confidence interval table contains the required columns
  expect_true(all(c("low beta", "mean b2", "high beta") %in% colnames(result$citable)))

})
