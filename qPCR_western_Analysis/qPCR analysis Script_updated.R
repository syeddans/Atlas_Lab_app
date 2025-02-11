library(performance)
library(dplyr)
library(tidyr)
library(plotrix)
library(multcompView)
library(plotly)
library(crosstalk)
library(tm)
library(ggplot2)
library(forcats)
library(matrixStats)
library(outliers)

options(dplyr.summarise.inform = FALSE)

if (exists("output_dir")) {
  log_file <<- file.path(output_dir, "analysis_log.txt")
  cat("Analysis Log\n", "Started at:", Sys.time(), "\n\n", file = log_file)
}else{
  log_file <<- file.path("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/analysis_log.txt")
  cat("Analysis Log\n", "Started at:", Sys.time(), "\n\n", file = log_file)
}


debug_script <- function(){
  output_dir <<- "/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/"
  control_gene <<- "b-actin"
  control_to_check <<- "mir+vehcon"
  analysis_type <<- "qPCR"
  unit <- "uM"
  chemicalsymbols <<- c()
  chemicalnames <<- c()
  stats_on <<- 'RQ' #RQ or DeltaCT 
  
  #formatted_data <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/qPCR_formatted_data_20250108_115413.csv")
  formatted_data <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/Phe_formatted_test_data.csv")
  formatted_data$chemical <- tolower(gsub(" ", "", formatted_data$chemical))

  final_df <- process_all_replicates(formatted_data, analysis_type, control_gene, control_to_check)

  CTdf <- process_CTdf(final_df, stats_on = "RQ", control_to_check = control_to_check, 
                      analysistype = analysis_type, unit = "uM", 
                      gene_or_protein_label = "gene")
  View(CTdf)
}

# Function to calculate average CT values for each gene and chemical
calculate_average_CT <- function(formatted_data) {
  formatted_data %>%
    group_by(gene, chemical, replicate) %>%
    summarize(average_CT = mean(CTvalue, na.rm = TRUE), .groups = "drop")
}

# Function to calculate Delta CT and RQ values
calculate_RQ <- function(average_df, analysistype, control_gene, control_to_check) {
  # Get control gene data - one value per chemical
  control_gene_data <- average_df %>%
    filter(gene == control_gene) %>%
    dplyr::select(chemical, average_CT, replicate) %>%
    rename(control_CT = average_CT)
  
  # Calculate Delta CT for all experimental genes
  deltaCT_df <- average_df %>%
    group_by(replicate) %>%
    filter(gene != control_gene) %>%
    # Join with control gene values by chemical
    left_join(control_gene_data, by = c("chemical", "replicate")) %>%
    # Calculate DeltaCT (experimental - control)
    mutate(DeltaCT = average_CT - control_CT) 

  # Calculate DDCT and RQ with error checking
  ddct_df <- deltaCT_df %>%
    group_by(gene, replicate) %>%
    mutate(
      control_deltaCT = DeltaCT[chemical == control_to_check],
      DDCT = DeltaCT - control_deltaCT
    ) %>%
    ungroup()
  
  # Calculate RQ with validation
  rq_df <- ddct_df %>%
    mutate(
      RQ = case_when(
        is.na(DDCT) ~ NA_real_,
        is.infinite(DDCT) ~ NA_real_,
        is.nan(DDCT) ~ NA_real_,
        TRUE ~ 2^(-DDCT)
      )
    )
  
  # Validate RQ values
  invalid_rq <- rq_df %>%
    filter(!is.na(RQ) & (is.infinite(RQ) | is.nan(RQ) | RQ < 0))
  
  if (nrow(invalid_rq) > 0) {
    warning("Invalid RQ values found:")
    print(invalid_rq)
  }
  
  return(rq_df)
}

# Function to process all replicates from formatted data
process_all_replicates <- function(formatted_data, analysistype, control_gene, control_to_check) {
  CTdf <- data.frame()
  
  formatted_data <- formatted_data %>% filter(!is.na(gene), !is.na(chemical), !is.na(CTvalue), !is.na(replicate))
  average_df <- calculate_average_CT(formatted_data)
  write.csv(average_df, file.path(output_dir, "Data_averaged.csv"), row.names = FALSE)
  rq_result_df <- calculate_RQ(average_df, analysistype, control_gene, control_to_check)
  
  write.csv(rq_result_df, file.path(output_dir, "Calculated_RQ_by_rep.csv"), row.names = FALSE)
  
  # Split data by replicate and rename columns
  replicate_dfs <- split(rq_result_df, rq_result_df$replicate)
  renamed_dfs <- lapply(names(replicate_dfs), function(rep_num) {
    df <- replicate_dfs[[rep_num]]
    names(df)[names(df) == "RQ"] <- paste0("RQ", rep_num)
    df$replicate <- NULL # Remove replicate column as it's now in column names
    return(df)
  })
  
  # Merge all replicate dataframes together
  CTdf <- Reduce(function(x, y) {
    merge(x, y, by = c("chemical", "gene"), all = TRUE)[,c("chemical", "gene", grep("^RQ", names(merge(x, y, by = c("chemical", "gene"), all = TRUE)), value = TRUE))]
  }, renamed_dfs)
  return(CTdf)
}

write_log <- function(...) {
    cat(..., "\n", file = log_file, append = TRUE)
  }

# Helper function to perform outlier test and set outliers to NA
perform_outlier_test <- function(row_data, stats_on_columns) {
  row_data[stats_on_columns] <- as.numeric(trimws(row_data[stats_on_columns]))
  values_to_check_outliers <- as.numeric(unlist(row_data[stats_on_columns]))
  
  # Log the chemical and gene being tested
  write_log(paste("\nTesting for outliers:"))
  write_log(paste("Chemical:", row_data["chemical"]))
  write_log(paste("Gene:", row_data["gene"]))
  write_log(paste("Values being tested:", paste(values_to_check_outliers, collapse=", ")))
  
  grubbs_result <- tryCatch({
    grubbs_result <- grubbs.test(values_to_check_outliers)
    #only try and remove the data point if test is successful and to avoid removing values when reps are consistent which would show 
    # 0 variablitiy and therefore 0 pvalue like in controls and also dont error out trying to remove NA when there are less than 3 reps 
    if (grubbs_result$p.value > 0 & grubbs_result$p.value < 0.05 & all(!is.na(values_to_check_outliers))){ 
      matches <- gregexpr("\\b[-\\d]+\\.\\d+\\b", grubbs_result$alternative, perl = TRUE)
      numbers <- regmatches(grubbs_result$alternative, matches)[[1]]
      outlier_to_remove <- as.numeric(numbers)
      
      # Find which replicate contains the outlier
      outlier_rep <- which(values_to_check_outliers == outlier_to_remove)
      write_log(paste("Outlier detected:", outlier_to_remove, "in replicate", outlier_rep))
      
      row_data[row_data == outlier_to_remove] <- NA
    } else {
      write_log("No outliers detected")
    }
  }, error = function(e) {
    # If an error occurs in the test, return the data as is
    write_log(paste("Error in Grubbs test (Outlier test):", e$message))
    
    return(row_data)
  })
  return(row_data)
}

# Apply outlier test to CTdf
apply_outlier_test <- function(CTdf, stats_prefix) {
  # Select columns based on stats_prefix
  write_log("Applying outlier test to CTdf")
  stats_on_columns <- names(CTdf)[startsWith(names(CTdf), stats_prefix)]
  write_log(paste("Stats on columns:", paste(stats_on_columns, collapse=", ")))
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
  
  # Convert selected_columns to numeric if they aren't already
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
  df$graphing_sort <- paste(df$chemicalgroup, df[[gene_or_protein_label]])
  df
}

# Function to calculate p-values and significance codes
calculate_pvalues <- function(CTdf, stats_on, control_to_check) {
  write_log("Analysis Parameters:")
  write_log("Control condition:", control_to_check)
  write_log("\n-------------------\n")
  
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
        write_log("\nAnalyzing:")
        write_log("Chemical:", chemical)
        write_log("Gene:", gene)
        
        subset_df <- CTdf[CTdf$chemicalgroup == chemical & CTdf$gene == gene, ]
        control_subset <- CTdf[CTdf$chemicalgroup == control_to_check & CTdf$gene == gene, ]
        #dont do stats on control or chemicals with only one rep 
      
        # Remove rows with NA in SEM column, meaning it only have 1 valid rep 
        subset_df <- subset_df[!is.na(subset_df$SEM), ]
      
        write_log("\nData Summary:")
        write_log("Number of treatment samples:", nrow(subset_df))
        write_log("Number of control samples:", nrow(control_subset))
        
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
          
          if (model_type == "Mixed Effects Model") {
            # 1. Normality of residuals
            model_resid <- residuals(model, type = "pearson")
            shapiro_test <- shapiro.test(model_resid)
            if (shapiro_test$p.value < 0.05) {
              failed_tests <- c(failed_tests, "Normality of residuals")
            }
            
            # 2. Homogeneity of variance (homoscedasticity)
            fitted_vals <- fitted(model)
            groups <- cut(fitted_vals, breaks = 3)
            levene_result <- car::leveneTest(model_resid ~ groups)
            if (levene_result$`Pr(>F)`[1] < 0.05) {
              failed_tests <- c(failed_tests, "Homoscedasticity")
            }
            
            # 3. Linearity - check if residuals have non-linear patterns
            cor_test <- cor.test(fitted_vals, model_resid)
            if (cor_test$p.value < 0.05) {
              failed_tests <- c(failed_tests, "Linearity")
            }
            
            # 4. Random effects normality
            rand_effects <- ranef(model)[[1]]  # Extract random effects
            rand_shapiro <- shapiro.test(as.vector(rand_effects))
            if (rand_shapiro$p.value < 0.05) {
              failed_tests <- c(failed_tests, "Random effects normality")
            }
            
            # 5. Check for influential observations
            influence <- influence(model, groups="Replicate")
            cook_d <- cooks.distance(influence)
            if (any(cook_d > 4/length(cook_d))) {  # Common threshold
              failed_tests <- c(failed_tests, "Influential observations")
            }
          } else if (model_type == "LM") {
            # Original Linear Model checks
            normality_result <- check_normality(model)
            p_value <- as.numeric(gsub(".*p = ([0-9.]+).*", "\\1", normality_result))
            if (p_value < 0.05) {
              failed_tests <- c(failed_tests, "Normality of residuals")
            }
            
            homoscedasticity_result <- check_heteroscedasticity(model)
            p_value <- as.numeric(gsub(".*p = ([0-9.]+).*", "\\1", homoscedasticity_result))
            if (p_value < 0.05) {
              failed_tests <- c(failed_tests, "Homoscedasticity")
            }
            
          } else if (model_type == "Generalized Linear Model") {
            # Original GLM checks
            residual_deviance <- sum(residuals(model, type = "deviance")^2)
            df_residual <- df.residual(model)
            if (residual_deviance/df_residual > 2) {
              failed_tests <- c(failed_tests, "Residual deviance")
            }
            
            if (family(model)$family %in% c("poisson", "binomial")) {
              phi <- residual_deviance/df_residual
              if (phi > 1.5) {
                failed_tests <- c(failed_tests, "Overdispersion")
              }
            }
          }
          
          if (length(failed_tests) == 0) {
            return(TRUE)
          } else {
            return(failed_tests)
          }
        }
        
        # Define all models in a single list
        models <- list(
          "Linear Model" = list(
            fit = function() lm(Expression ~ Treatment + Replicate, data = balanced_df),
            type = "LM"
          ),
          "Mixed Effects Model" = list(
            fit = function() {
              lmer(Expression ~ Treatment + (1|Replicate), 
                    data = balanced_df,
                    REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa",
                                        optCtrl = list(maxfun = 100000),
                                        check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.nRE = "ignore"))
            },
            type = "Mixed Effects Model"
          ),
          "GLM Gaussian" = list(
            fit = function() glm(Expression ~ Treatment + Replicate, 
                                data = balanced_df, 
                                family = gaussian(link = "identity")),
            type = "Generalized Linear Model"
          ),
          "GLM Poisson" = list(
            fit = function() glm(Expression ~ Treatment + Replicate, 
                                data = balanced_df, 
                                family = poisson(link = "log")),
            type = "Generalized Linear Model"
          ),
          "GLM Gamma" = list(
            fit = function() glm(Expression ~ Treatment + Replicate, 
                                data = balanced_df, 
                                family = Gamma(link = "log")),
            type = "Generalized Linear Model"
          ),
          "GLM Inverse Gaussian" = list(
            fit = function() glm(Expression ~ Treatment + Replicate, 
                                data = balanced_df, 
                                family = inverse.gaussian(link = "1/mu^2")),
            type = "Generalized Linear Model"
          )
        )
        
        # Function to evaluate models
        evaluate_models <- function() {
          results <- list()
          valid_models <- numeric()
          
          for (model_name in names(models)) {
            write_log("\nEvaluating", model_name, "...")
            
            tryCatch({
              model <- models[[model_name]]$fit()
              assumption_results <- check_assumptions(model, model_type = models[[model_name]]$type)
              
              if (isTRUE(assumption_results)) {
                write_log("✓ Model passed all assumptions")
                valid_models[model_name] <- AIC(model)
                
                write_log("\nModel Summary:")
                # Capture summary output
                sum_output <- capture.output(print(summary(model)))
                write_log(paste(sum_output, collapse = "\n"))
                write_log("AIC:", AIC(model))
                
                results[[model_name]] <- model
              } else {
                write_log("✗ Model failed assumptions:", paste(assumption_results, collapse = ", "))
              }
              
            }, error = function(e) {
              write_log("✗ Error fitting model:", conditionMessage(e))
            })
          }
          
          if (length(valid_models) > 0) {
            best_model_name <- names(valid_models)[which.min(valid_models)]
            best_model <- results[[best_model_name]]
            
            if (best_model_name == "Mixed Effects Model") {
              library(lmerTest)
              p_values <- summary(best_model)$coefficients[-1, "Pr(>|t|)"]
            } else {
              p_values <- summary(best_model)$coefficients[-1, 4]
            }
            
            write_log("\nBest model:", best_model_name)
            write_log("P-values:", paste(p_values, collapse = ", "))
            
            return(list(
              model_name = best_model_name,
              p_values = p_values
            ))
          } else {
            write_log("\nNo valid models found.")
            return(NULL)
          }
        }
        
        model_results <- evaluate_models()
        
        if (!is.null(model_results)) {
          # Get the number of treatments tested (excluding control)
          num_treatments <- nrow(subset_df)
          
          # Only take the first num_treatments p-values from the model results
          treatment_p_values <- model_results$p_values[1:num_treatments]
          
          # Only assign p-values to rows that were actually tested (have valid SEM)
          CTdf$pvalue[CTdf$gene == gene & 
                     CTdf$chemicalgroup == chemical & 
                     !is.na(CTdf$SEM)] <- treatment_p_values
          
          # Similarly, only assign model name to rows that were tested
          CTdf$model_name[CTdf$gene == gene & 
                         CTdf$chemicalgroup == chemical & 
                         !is.na(CTdf$SEM)] <- model_results$model_name
        }
      }
    }
  }
  CTdf$pvalue <- as.numeric(CTdf$pvalue)
  CTdf$pvalue <- round(CTdf$pvalue, 4)  
  CTdf$pvalue <- format(CTdf$pvalue, scientific = FALSE)
  
  get_sigcode <- function(pvalue) {
    # Debugging output
    pvalue <- as.numeric(pvalue)
    # Check if pvalue is NA or a character representation of NA
    if (is.na(pvalue) || pvalue == "" || pvalue == "NA" || pvalue == "NaN") {
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
  CTdf <- calculate_pvalues(CTdf, stats_on, control_to_check)
  CTdf <- adjust_labels(CTdf, analysistype, gene_or_protein_label)
}


#debug_script()
