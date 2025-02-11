# Test process_all_replicates function
library(here)
library(testthat)
library(devtools)
library(usethis)
source(here("qPCR_western_Analysis", "qPCR analysis Script_updated.R"))

test_that("process_all_replicates functions correctly", {
  print("testing process_all_replicates")
  test_data <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Test_Data_Results/input_data_formatted.csv")
  expected_result <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Test_Data_Results/RQs.csv")
  
  result <- process_all_replicates(
    test_data, 
    analysistype = "qPCR",
    control_gene = "b-actin",
    control_to_check = "mir+vehcon"  
  )
  
  expect_equal(
    result[order(result$chemical, result$gene), ],
    expected_result[order(expected_result$chemical, expected_result$gene), ],
    tolerance = 1e-4
  )
  expect_true(all(c("chemical", "gene", "RQ") %in% colnames(result)))
})

# Test process_CTdf function
test_that("process_CTdf functions correctly", {
  test_data <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Test_Data_Results/RQs.csv")
  expected_result <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Test_Data_Results/Final_results_with_stats.csv")

  print("testing process_CTdf")
  sink(nullfile())
  result <- suppressWarnings(
    process_CTdf(
      test_data,
      stats_on = "RQ", 
      control_to_check = "mir+vehcon",
      analysistype = "qPCR",
      unit = "uM",
      gene_or_protein_label = "gene"
  ))
  sink()
  
  # Compare only the relevant columns with correct case
  result <- result[, c("chemical", "gene", "DeltaCT", "RQ", "average_RQ", "SEM", 
                      "concentration", "chemicalgroup", "pvalue", "sigcode", "graphing_sort")]
  expected_result <- expected_result[, c("chemical", "gene", "DeltaCT", "RQ", "average_RQ", "SEM",
                                       "concentration", "chemicalgroup", "pvalue", "sigcode", "graphing_sort")]
  
  # Convert gene names to lowercase and chemical to factor
  result$gene <- tolower(result$gene)
  expected_result$gene <- tolower(expected_result$gene)
  
  # Convert chemical to factor and concentration to numeric
  result$chemical <- as.factor(result$chemical)
  expected_result$chemical <- as.factor(expected_result$chemical)
  result$concentration <- as.numeric(as.character(result$concentration))
  expected_result$concentration <- as.numeric(as.character(expected_result$concentration))
  
  # Drop unused levels
  result <- droplevels(result)
  expected_result <- droplevels(expected_result)
  
  # Sort both dataframes using numeric concentration
  result <- result[order(result$graphing_sort, as.numeric(result$concentration), result$gene), ]
  expected_result <- expected_result[order(expected_result$graphing_sort, as.numeric(expected_result$concentration), expected_result$gene), ]
  
  # Original comparison
  expect_equal(
    result,
    expected_result,
    tolerance = 1e-4,
    check.attributes = FALSE
  )
})

# Test calculate_deltaCT function
test_that("calculate_RQ functions correctly", {
  test_data <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Test_Input/input_data_formatted.csv")
  expected_result <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Test_Input/deltaCT_RQs.csv")
  
  print("testing calculate_RQ")
  average_df_test <- calculate_average_CT(test_data)
  result <- calculate_RQ(average_df_test,"qPCR", "b-actin", "mie")
  # Test first row valus
  expect_equal(result$gene[1], "ADIPOQ")
  expect_equal(result$chemical[1], "0.01umdeha") 
  expect_equal(result$average_CT[1], 25.97, tolerance = 0.0001)
  expect_equal(result$control_CT[1], 18.4067, tolerance = 0.0001)
  expect_equal(result$DeltaCT[1], -7.5633, tolerance = 0.0001)
  expect_equal(result$RQ[1], 2.6386, tolerance = 0.0001)
})


