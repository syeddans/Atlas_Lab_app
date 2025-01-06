library(performance)
library(lmtest)
library(car)
library(glmnet)
library(lme4)
library(dplyr)
library(tm)
library(ggplot2)
library(forcats)
library(plotrix)
library(matrixStats)
library(remotes)
library("multcompView")
library("plotly")
library("crosstalk")
library("scatterD3")
library(tidyr)
#remotes::install_github("snandi/RFunctionsSN")
library(RFunctionsSN)
#devtools::install_github("jcheng5/d3scatter")
library(tibble)
options(dplyr.summarise.inform = FALSE)

library(dplyr)
library("DescTools")
library("CATT")
library("coin")
library(plotrix)
library(multcomp)
library("aod")
library(outliers)
library(lmPerm)

library(MASS)
library(glmnet)
library(lme4)
library(Metrics)
#rep_files_dir <- "~/storage/qPCR Analysis/Phenanthrene western"
#output_dir <- "~/storage/Atlas_lab_app"
#control_gene <- "b-actin"
#control_to_check <- "mir+vehcon"
#analysis_type <- "western"
unit <- "uM"
#chemicalsymbols <- c("phe","9p","bap") #if using full names already leave empty 
#chemicalnames<- c("phe","9p","bap")
#stats_on <- 'RQ' #RQ or DeltaCT

#output_dir = output_dir
plate <-0
CTdf <- data.frame()

# Load required packages
library(dplyr)
library(tidyr)

# Function to read and parse files in a plate directory
parse_plate_files <- function(plate) {
  plate_files <- list.files(plate, pattern = NULL, all.files = FALSE, full.names = TRUE)
  platesoutput <- data.frame(value = character(), well = character())
  
  for (file in plate_files) {
    #print(file)
    dftemp <- data.frame()
    table <- read.csv(file, sep = ",", header = FALSE, fill = TRUE)
    colnames(table) <- table[1,]
    
    for (row in 2:nrow(table)) {
      for (col in 2:ncol(table)) {
        dftemp <- rbind(dftemp, data.frame(value = table[row, col], well = paste0(table[row, 1], table[1, col])))
      }
    }
    
    if (grepl("-CTmap", file, ignore.case = TRUE)) {
      colnames(dftemp) <- c("CTvalue", "well")
    } else if (grepl("-wellmap", file, ignore.case = TRUE)) {
      colnames(dftemp) <- c("chemical", "well")
    } else if (grepl("-genemap", file, ignore.case = TRUE) || grepl("-proteinmap", file, ignore.case = TRUE)) {
      colnames(dftemp) <- c("gene", "well")
    }
    
    if (nrow(platesoutput) == 0) {
      platesoutput <- dftemp
    } else {
      platesoutput <- left_join(platesoutput, dftemp, by = "well")
    }
  }
  
  platesoutput <- platesoutput %>% distinct() %>% drop_na()
  platesoutput$chemical <- tolower(gsub(" ", "", platesoutput$chemical))
  platesoutput$CTvalue <- as.numeric(platesoutput$CTvalue)
  platesoutput
}

# Function to calculate average CT values for each gene and chemical
calculate_average_CT <- function(platesoutput) {
  platesoutput %>%
    group_by(gene, chemical) %>%
    summarize(average_CT = mean(CTvalue, na.rm = TRUE), .groups = "drop")
}

# Function to calculate Delta CT and RQ values
calculate_RQ <- function(average_df, analysistype, control_to_check) {
  bactin <- filter(average_df, gene == control_gene)
  gene_of_interest <- filter(average_df, gene != control_gene)
  joined <- left_join(bactin, gene_of_interest, by = 'chemical')
  joined <- joined %>%
    mutate(control = ifelse(chemical == control_to_check, "C", ""))
  if (analysistype == "western") {
    joined <- joined %>%
      mutate(DeltaCT = round(average_CT.y / average_CT.x, digits = 2)) %>%
      group_by(gene.y) %>%
      mutate(RQ = round(DeltaCT / DeltaCT[control == "C"], digits = 2))
  } else {
    joined <- joined %>%
      mutate(DeltaCT = round(average_CT.y - average_CT.x, digits = 2)) %>%
      group_by(gene.y) %>%
      mutate(
        control_DeltaCT = ifelse(any(control == "C"), DeltaCT[control == "C"], NA),
        DDCT = round(DeltaCT - control_DeltaCT, digits = 2),
        RQ = round(2^(-DDCT), digits = 2)
      ) %>%
      ungroup()
  }
  rename(joined, gene = gene.y)
}

# Function to process each replicate and update the results dataframe
process_replicate <- function(rep, analysistype, control_to_check) {
  plates_in_rep <- list.dirs(rep, recursive = FALSE)
  rep_result_df <- data.frame("chemical" = character(0), "gene" = character(0), "RQ" = numeric(0))
  for (plate in plates_in_rep) {
    platesoutput <- parse_plate_files(plate)
    average_df <- calculate_average_CT(platesoutput)
    joined <- calculate_RQ(average_df, analysistype, control_to_check)
    rep_result_df <- rbind(rep_result_df, joined[, c("chemical", "gene", "DeltaCT", "RQ")])
  }
  rep_result_df %>%
    group_by(chemical, gene) %>%
    summarize(RQ = mean(RQ, na.rm = TRUE), DeltaCT = mean(DeltaCT, na.rm = TRUE), .groups = "drop")
}

# Main loop to iterate over replicates
process_all_replicates <- function(reps_dir, analysistype, control_to_check) {
  CTdf <- data.frame()
  rep_num <- 0
  rep_folders <- list.dirs(reps_dir, recursive = FALSE)
  for (rep in rep_folders) {
    print(rep)
    rep_num <- rep_num + 1
    rep_result_df <- process_replicate(rep, analysistype, control_to_check)
    
    if (nrow(CTdf) == 0) {
      CTdf <- bind_rows(CTdf, rep_result_df)
    } else {
      RQ_rep <- paste0("RQ", as.character(rep_num))
      DeltaCT_rep <- paste0("DeltaCT", as.character(rep_num))
      names(rep_result_df)[names(rep_result_df) == "RQ"] <- RQ_rep
      names(rep_result_df)[names(rep_result_df) == "DeltaCT"] <- DeltaCT_rep
      
      nonmatching <- anti_join(rep_result_df, CTdf, by = c("chemical", "gene"))
      if (nrow(nonmatching) != 0) {
        CTdf <- bind_rows(CTdf, nonmatching)
        rep_result_df <- anti_join(rep_result_df, nonmatching, by = c("chemical", "gene"))
      }
      
      CTdf<- left_join(CTdf, rep_result_df, by= c("chemical", "gene"))
      columns_to_merge <- grep(paste0("^",RQ_rep), names(CTdf), value = TRUE)
      columns_to_merge_DeltaCT <- grep(paste0("^",DeltaCT_rep), names(CTdf), value = TRUE)
      # Merge selected columns into a new column
      # Remove the original columns if needed
      if(length(columns_to_merge)>1){ 
        RQ_toadd<- coalesce(!!!c(CTdf[columns_to_merge]))
        DeltaCT_toadd <- coalesce(!!!c(CTdf[columns_to_merge_DeltaCT]))
        RQ_tobind <- data.frame(value = RQ_toadd)
        DeltaCT_tobind <- data.frame(value = DeltaCT_toadd)
        colnames(RQ_tobind) <- RQ_rep
        colnames(DeltaCT_tobind) <- DeltaCT_rep
        CTdf <- CTdf[, !names(CTdf) %in% columns_to_merge]
        CTdf <- CTdf[, !names(CTdf) %in% columns_to_merge_DeltaCT]
        CTdf<-bind_cols(CTdf,RQ_tobind) 
        CTdf<-bind_cols(CTdf,DeltaCT_tobind)
      }

      
    }
  }
  
  CTdf
}

# Helper function to perform outlier test and set outliers to NA
perform_outlier_test <- function(row_data, stats_on_columns) {
  row_data[stats_on_columns] <- as.numeric(trimws(row_data[stats_on_columns]))
  values_to_check_outliers <- as.numeric(unlist(row_data[stats_on_columns]))
  grubbs_result <- tryCatch({
    grubbs_result <- grubbs.test(values_to_check_outliers)
    #only try and remove the data point if test is successful and to avoid removing values when reps are consistent which would show 
    # 0 variablitiy and therefore 0 pvalue like in controls and also dont error out trying to remove NA when there are less than 3 reps 
    if (grubbs_result$p.value > 0 & grubbs_result$p.value < 0.05 & all(!is.na(values_to_check_outliers))){ 
      matches <- gregexpr("\\b[-\\d]+\\.\\d+\\b", grubbs_result$alternative, perl = TRUE)
      numbers <- regmatches(grubbs_result$alternative, matches)[[1]]
      outlier_to_remove <- as.numeric(numbers)
      row_data[row_data == outlier_to_remove] <- NA
      }
  }, error = function(e) {
    # If an error occurs in the test, return the data as is
    print("error in grubbs test (Outlier test)")
    print(row_data["chemical"])
    return(row_data)
  })
  return(row_data)
}


# Apply outlier test to CTdf
apply_outlier_test <- function(CTdf, stats_prefix) {
  # Select columns based on stats_prefix
  stats_on_columns <- names(CTdf)[startsWith(names(CTdf), stats_prefix)]
  
  # Apply the outlier test to each row and ensure the result is a data frame
  CTdf <- as.data.frame(t(apply(CTdf, 1, function(row_data) {
    # Call the outlier test function for each row
    row_data <- perform_outlier_test(row_data, stats_on_columns)
    return(row_data)  # Return modified row
  })))
  
  return(CTdf)
}

# Round SEM values
calculate_sem <- function(df, stats_prefix) {
  # Select columns that match the stats_prefix and ensure they are numeric
  selected_columns <- df[, grepl(stats_prefix, names(df))]
  
  # Convert selected columns to numeric if they aren't already
  selected_columns <- data.frame(lapply(selected_columns, as.numeric))
  
  # Calculate average RQ (row means)
  df$average_RQ <- round(rowMeans(selected_columns, na.rm = TRUE), digits = 2)
  
  # Calculate SEM for each row
  df$SEM <- apply(selected_columns, 1, function(row) std.error(na.omit(row)))
  
  return(df)
}

# Function to extract concentration from a string
extract_concentration <- function(string) {
  result <- regmatches(string, regexpr("^\\d*\\.?\\d+", string))
  if (length(result) == 0) NA else result
}

#TODO: add symbol to full name functionality 
replace_symbols <- function(string) {
  if(length(chemicalsymbols)>0){ 
    for (i in seq_along(chemicalsymbols)) {
      if (grepl(chemicalsymbols[i], string)) {
        return(chemicalnames[i])
      }
    }
  }
  else{ 
    for (i in seq_along(chemicalnames)){ 
      if (grepl(chemicalnames[i], string, ignore.case = TRUE)) {
        return(chemicalnames[i])
        }
      }
    }
  
  return(string)
}

# Updated function with replace_symbols included
add_concentration_and_groups <- function(df) {
  df$concentration <- as.double(unlist(lapply(df$chemical, extract_concentration)))
  #df$chemicalgroup <- unlist(lapply(df$chemical, replace_symbols))
  df$chemicalgroup <- as.character(gsub("^[0-9.]+(um)?", "", df$chemical))
  return(df)
}

# Adjust gene or protein labeling based on analysis type
adjust_labels <- function(df, analysistype, gene_or_protein_label) {
  if (analysistype == "western") {
    gene_or_protein_label <- "protein"
  }
  names(df)[names(df) == "gene"] <- gene_or_protein_label
  df[[gene_or_protein_label]] <- tolower(df[[gene_or_protein_label]])
  df$graphing_sort <- paste(df$chemicalgroup, df[[gene_or_protein_label]], sep = "_")
  df
}

# Function to calculate p-values and significance codes
calculate_pvalues <- function(CTdf, stats_on, control_to_check) {
  CTdf$pvalue <- NA
  CTdf$model_name <- NA
  
  CTdf$chemical <- factor(
    CTdf$chemical,
    levels = c(
      CTdf %>%
        arrange(as.numeric(concentration)) %>%
        pull(chemical) %>%
        unique()
    )
  )
  CTdf <- CTdf %>% arrange(chemical)
  
  unique_chemicals <- unique(CTdf$chemicalgroup)
  # Create empty vectors to store p-values and model names
  pvalues <- vector("numeric", length = nrow(CTdf))
  model_names <- vector("character", length = nrow(CTdf))
  #TODO: currently if there is little to no variation like 0.07 vs 0 in expression, linear model 
  for (chemical in unique_chemicals) {
    if (chemical != control_to_check) {
      for (gene in unique(CTdf$gene)) {
        print(chemical)
        print(gene)
        x<-chemical
        z<-gene
        subset_df <- CTdf[CTdf$chemicalgroup == chemical & CTdf$gene == gene, ]
        control_subset <- CTdf[CTdf$chemicalgroup == control_to_check & CTdf$gene == gene, ]
        #dont do stats on control or chemicals with only one rep 
        if (nrow(subset_df) > 0) { 
          subset_df_with_control <- rbind(control_subset,subset_df)
          stats_columns <- names(subset_df_with_control)[startsWith(names(subset_df_with_control), stats_on)]
          subset_df_with_control <- subset_df_with_control[, c("chemical", "concentration", stats_columns)]
          
          subset_df_with_control$non_na_count <- rowSums(!is.na(subset_df_with_control[, grep("RQ", colnames(subset_df_with_control))]))
          # Step 2: Restructure data to have "Treatment" and "Expression" columns

           long_df <- subset_df_with_control %>%
            pivot_longer(
              cols = starts_with("RQ"),
              names_to = "Replicate",
              values_to = "Expression"
            ) %>%
            mutate(
              Replicate = if_else(Replicate == "RQ", 1L, as.integer(gsub("RQ", "", Replicate))),
              Expression = as.numeric(Expression)# Ensure Expression is numeric
            )
          # Step 3: Ensure equal rows per Treatment (fill with NA if needed)
          treatments <- unique(long_df$chemical)
          max_reps <- max(table(long_df$chemical))

          balanced_df <- long_df %>%
            group_by(chemical) %>%
            complete(Replicate = 1:max_reps, fill = list(Expression = NA))
          balanced_df <- balanced_df %>%
            rename(Treatment = chemical)
          # Ensure 'Treatment' is a factor and that control is first
          # Arrange Treatment levels based on concentration, keeping control as the first level
          balanced_df$Treatment <- factor(
            balanced_df$Treatment,
            levels = c(
              control_to_check, 
              levels(balanced_df$Treatment)[levels(balanced_df$Treatment) != control_to_check]
            )
          )
          balanced_df <- balanced_df %>% 
            arrange(Treatment)
          # Remove rows with NAs for modeling
          balanced_df <- balanced_df %>%
            filter(!is.na(Expression))
          
     
          check_assumptions <- function(model, model_type) {
            failed_tests <- c()
            # 1. Normality of residuals for Linear Models (LM) and Mixed Effects Models
            if (model_type == "LM" || model_type == "MixedEffects") {
              normality_result <- check_normality(model)
              # Extract the p-value from the result message
              p_value <- as.numeric(gsub(".*p = ([0-9.]+).*", "\\1", normality_result))
              if (p_value < 0.05) {
                failed_tests <- c(failed_tests, "Normality of residuals")
              }
            }
            
            # 2. Normality of residuals for GLMs (using simulation)
            if (model_type == "GLM") {
              # Simulate residuals and check uniformity
              simulated_residuals <- simulate_residuals(model)
              residual_check_result <- check_residuals(simulated_residuals)
              # Extract the p-value from the result message
              p_value <- as.numeric(gsub(".*p = ([0-9.]+).*", "\\1", residual_check_result))
              if (p_value < 0.05) {
                failed_tests <- c(failed_tests, "Normality of residuals")
              }
            }
            
            # 3. Homoscedasticity (for LM/GLM/MixedEffects, use appropriate method)
            if (model_type %in% c("LM", "MixedEffects")) {
              homoscedasticity_result <- check_heteroscedasticity(model)
              # Extract the p-value from the result message
              p_value <- as.numeric(gsub(".*p = ([0-9.]+).*", "\\1", homoscedasticity_result))
              if (p_value < 0.05) {
                failed_tests <- c(failed_tests, "Homoscedasticity")
              }
            }
            
            # For GLMs that are not Gaussian (Poisson, Gamma, etc.), use residuals for heteroscedasticity check
            if (model_type == "GLM" && !is.null(family(model)$family) && family(model)$family != "gaussian") {
              # For Poisson or Gamma, check residuals behavior
              residuals <- residuals(model, type = "pearson")
              # Use residual plots or other checks for heteroscedasticity
              # Here, you can replace this with your method to check residuals distribution/variance.
              if (sd(residuals) > 2) {
                failed_tests <- c(failed_tests, "Heteroscedasticity (Poisson/Gamma)")
              }
            }
            
            
            # 4. Independence (Durbin-Watson for LM/GLM)
            if (model_type %in% c("LM", "GLM", "MixedEffects")) {
              dw_test <- durbinWatsonTest(model)
              if (dw_test$p < 0.05) {
                failed_tests <- c(failed_tests, "Independence of residuals (Durbin-Watson test)")
              }
            }
            
            # 5. Overdispersion (for GLM)
            if (model_type == "GLM") {
              residual_deviance <- sum(residuals(model)^2)
              overdispersion_ratio <- residual_deviance / df.residual(model)
              if (overdispersion_ratio > 1) {
                failed_tests <- c(failed_tests, "Overdispersion")
              }
            }
            
            if (length(failed_tests) == 0) {
              return(TRUE)  # All assumptions passed
            } else {
              return(failed_tests)  # Return the names of failed tests
            }
          }
          
          # Define model families for GLM
          glm_families <- list(
            "Gaussian Identity" = gaussian(link = "identity"),
            "Binomial Logit" = binomial(link = "logit"),
            "Poisson Log" = poisson(link = "log"),
            "Gamma Log" = Gamma(link = "log"),
            "Inverse Gaussian 1/mu^2" = inverse.gaussian(link = "1/mu^2")
          )
          
          # Define model types, including testing all GLM families
          #TODO: be able to choose which models you want to test out of
          models <- list(
           # "Linear Model" = function() lm(Expression ~ Treatment, data = balanced_df),
            "Generalized Linear Model" = function() {
              glm_results <- list()
              for (family_name in names(glm_families)) {
                tryCatch({
                  model <- glm(Expression ~ Treatment, data = balanced_df, family = glm_families[[family_name]])
                  assumption_results <- check_assumptions(model, model_type = "GLM")
                  if (isTRUE(assumption_results)) {
                    glm_results[[family_name]] <- list(
                      model = model,
                      AIC = AIC(model)
                    )
                  } else {
                    glm_results[[family_name]] <- list(
                      model = model,
                      AIC = NA,
                      failed_tests = assumption_results
                    )
                  }
                }, error = function(e) {
                  message("Error fitting GLM with family ", family_name, ": ", e$message)
                  glm_results[[family_name]] <- list(
                    model = NULL,
                    AIC = NA,
                    error_message = e$message
                  )
                })
              }
              return(glm_results)
            }#,
            #"Mixed Effects Model" = function() lmer(Expression ~ Treatment + (1 | Treatment), data = balanced_df)
          )
          
          # Function to evaluate models
            results <- list()
            
            for (model_name in names(models)) {
              tryCatch({
                if (model_name == "Generalized Linear Model") {
                  glm_results <- models[[model_name]]()  # Get results for all GLM families
                  for (family_name in names(glm_results)) {
                    glm_family_result <- glm_results[[family_name]]
                    model <- glm_family_result$model
                    print(family_name)
                    print(summary(model))
                    
                    if (!is.null(model) && !is.na(glm_family_result$AIC)) {
                      results[[paste(model_name, family_name, sep = " - ")]] <- glm_family_result$AIC
                    } else {
                      message("Model failed assumptions or encountered error: ", family_name, " - Failed Tests: ", paste(assumption_results, collapse = ", "))
                    }
                  }
                } else {
                  # Fit other models
                  model <- models[[model_name]]()
                  assumption_results <- check_assumptions(model, model_type = ifelse(model_name == "Linear Model", "LM", "Mixed"))
                  if (isTRUE(assumption_results)) {
                    results[[model_name]] <- AIC(model)
                    
                  } else {
                    message("Model failed assumptions: ", model_name, " - Failed Tests: ", paste(assumption_results, collapse = ", "))
                  }
                }
              }, error = function(e) {
                message("Error fitting model: ", model_name, " - ", e$message)
              })
            }
            
            # Determine the best model based on AIC
            if (length(results) > 0) {
              valid_results <- unlist(results)
              best_model <- names(valid_results)[which.min(valid_results)]
              cat("Best model based on AIC:", best_model, "\n")
              browser()
              # Extract p-values for the treatment comparison
              if (startsWith(best_model, "Generalized Linear Model")) {
              
                family_name <- sub("Generalized Linear Model - ", "", best_model)
                model_fit <- models[["Generalized Linear Model"]]()[[family_name]]$model
                model_summary <- summary(model_fit)
                # Extract p-values from the correct column for GLM (-1 to avoid storing y intercdept p value)
                browser()
                l<-model_summary$coefficients[-1, 4] 
                browser()
                CTdf$pvalue[CTdf$gene == gene & CTdf$chemicalgroup == chemical] <- model_summary$coefficients[-1, 4]  # p-values from summary
               
                } else if (best_model == "Linear Model") {
                  browser()
                model_fit <- models[["Linear Model"]]()
                model_summary <- summary(model_fit)
                # Extract p-values for linear models
                CTdf$pvalue[CTdf$gene == gene & CTdf$chemicalgroup == chemical] <- model_summary$coefficients[-1, 4]  # p-values from summary
                browser()
                } else if (best_model == "Mixed Effects Model") {
                  browser()
                # Fit the mixed effects model
                model_fit <- models[["Mixed Effects Model"]]()
                # Use lmerTest to get p-values for the fixed effects
                library("lmerTest")
                model_summary <- summary(model_fit)
                # Extract p-values for fixed effects (assuming the 5th column for p-values)
                CTdf$pvalue[CTdf$gene == gene & CTdf$chemicalgroup == chemical] <- model_summary$coefficients[-1, 5]  # p-values for fixed effects
                browser()
                }
              browser()
              # Store the best model name for later
              CTdf$model_name[CTdf$gene == gene & CTdf$chemicalgroup == chemical] <- best_model
            } else {
              cat("All models failed assumption tests or encountered errors.\n")
            }
        }
      }
    }
  }
  browser()
  CTdf$pvalue <- round(CTdf$pvalue, 4)  
  CTdf$pvalue <- format(CTdf$pvalue, scientific = FALSE)
  
  get_sigcode <- function(pvalue) {
    if (is.na(pvalue)) {
      return(" ")  # Return a default significance code for missing values
    } else if (pvalue < 0.001) {
      return("***")
    } else if (pvalue < 0.01) {
      return("**")
    } else if (pvalue < 0.05) {
      return("*")
    } else if (pvalue < 0.1) {
      return("o")
    } else {
      return(" ")
    }
  }
  CTdf$sigcode <- sapply(CTdf$pvalue, get_sigcode)
  
  return(CTdf)
}



# Main function to process data
process_CTdf <- function(CTdf, stats_on, control_to_check, analysistype, unit, gene_or_protein_label) {
  CTdf <- CTdf[!duplicated(CTdf[, c("chemical", "gene")]), ]
  CTdf <- apply_outlier_test(CTdf, stats_on)
  CTdf <- calculate_sem(CTdf, "^RQ")
  CTdf <- add_concentration_and_groups(CTdf)
  browser()
  CTdf <- calculate_pvalues(CTdf, stats_on, control_to_check)
  CTdf <- adjust_labels(CTdf, analysistype, gene_or_protein_label)
  
  # Adjust concentration display and move columns
  #CTdf$chemical <- paste0(format(CTdf$concentration, nsmall = 1), unit, " ", as.character(CTdf$chemicalgroup), " ", CTdf[[gene_or_protein_label]])
  #CTdf <- CTdf[, c("chemical", "concentration", "chemicalgroup", setdiff(names(CTdf), c("chemical", "concentration", "chemicalgroup")))]
  
  CTdf
}

# Run processing
final_df <- process_all_replicates(rep_files_dir, analysis_type, control_to_check)
CTdf <- process_CTdf(final_df, stats_on = "RQ", control_to_check = control_to_check, analysistype = analysis_type, unit = "uM", gene_or_protein_label = "gene")

write.table(CTdf, paste0("~/storage/Atlas_lab_app/Columbus_Analysis/", "resall.txt"), sep = "\t", row.names = FALSE)
