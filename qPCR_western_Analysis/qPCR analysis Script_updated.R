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
# Load required packages
library(dplyr)
library(tidyr)
#rep_files_dir <- "~/storage/qPCR Analysis/Phenanthrene western" #update path to a formatted file 

#output_dir <- "~/storage/Atlas_lab_app"
#control_gene <- "b-actin"
#control_to_check <- "mie"
#analysis_type <- "qPCR"
#unit <- "uM"
#chemicalsymbols <- c() #if using full names already leave empty 
#chemicalnames<- c()
#stats_on <- 'RQ' #RQ or DeltaCT 

#output_dir = output_dir
plate <-0
CTdf <- data.frame()



# Function to calculate average CT values for each gene and chemical
calculate_average_CT <- function(formatted_data) {
  formatted_data %>%
    group_by(gene, chemical) %>%
    summarize(average_CT = mean(CTvalue, na.rm = TRUE), .groups = "drop")
}

# Function to calculate Delta CT and RQ values
calculate_RQ <- function(average_df, analysistype, control_to_check) {
  # Get control gene data - one value per chemical
  control_gene_data <- average_df %>%
    filter(gene == control_gene) %>%
    dplyr::select(chemical, average_CT) %>%
    rename(control_CT = average_CT)
  
  # Calculate Delta CT for all experimental genes
  deltaCT_df <- average_df %>%
    filter(gene != control_gene) %>%
    # Join with control gene values by chemical``
    left_join(control_gene_data, by = "chemical") %>%
    # Calculate DeltaCT (experimental - control)
    mutate(DeltaCT = average_CT - control_CT)
  
  # Calculate DDCT and RQ
  ddct_df <- deltaCT_df %>%
    group_by(gene) %>%
    mutate(
      control_deltaCT = DeltaCT[chemical == control_to_check],
      DDCT = DeltaCT - control_deltaCT
    ) %>%
    #select(-control_deltaCT) %>%
    ungroup()
  
  rq_df <- ddct_df %>%
    mutate(RQ = 2^(-DDCT))
  
  return(rq_df)
}

# Function to process all replicates from formatted data
process_all_replicates <- function(formatted_data, analysistype, control_to_check) {
  CTdf <- data.frame()
  
  # Get unique rep numbers
  reps <- unique(formatted_data$rep)
  rep_num <- 0
  
  #cat("Starting process_all_replicates function\n")
  #cat("Found reps:", paste(reps, collapse = ", "), "\n")
  
  for (current_rep in reps) {
    #cat("\nProcessing rep:", current_rep, "\n")
    rep_num <- rep_num + 1
    
    # Filter data for current rep
    rep_data <- formatted_data[formatted_data$rep == current_rep, ]
    #cat("\nRep data:\n")
    #print(rep_data)
    
    
    # Calculate average CT and RQ for this rep
    average_df <- calculate_average_CT(rep_data)
    #cat("average_df:\n")
    #print(average_df)
    rep_result_df <- calculate_RQ(average_df, analysistype, control_to_check)
    rep_result_df <- rep_result_df[, c("chemical", "gene", "DeltaCT", "RQ")]
    #message("\nProcessed results for rep:")
    #print(rep_result_df)
    
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
      
      CTdf <- left_join(CTdf, rep_result_df, by = c("chemical", "gene"))
      columns_to_merge <- grep(paste0("^", RQ_rep), names(CTdf), value = TRUE)
      columns_to_merge_DeltaCT <- grep(paste0("^", DeltaCT_rep), names(CTdf), value = TRUE)
      
      if(length(columns_to_merge) > 1) { 
        RQ_toadd <- coalesce(!!!c(CTdf[columns_to_merge]))
        DeltaCT_toadd <- coalesce(!!!c(CTdf[columns_to_merge_DeltaCT]))
        RQ_tobind <- data.frame(value = RQ_toadd)
        DeltaCT_tobind <- data.frame(value = DeltaCT_toadd)
        colnames(RQ_tobind) <- RQ_rep
        colnames(DeltaCT_tobind) <- DeltaCT_rep
        CTdf <- CTdf[, !names(CTdf) %in% columns_to_merge]
        CTdf <- CTdf[, !names(CTdf) %in% columns_to_merge_DeltaCT]
        CTdf <- bind_cols(CTdf, RQ_tobind) 
        CTdf <- bind_cols(CTdf, DeltaCT_tobind)
      }
      #message("\nCurrent CTdf:")
      #message(paste(capture.output(CTdf), collapse = "\n"))
    }
  }
  
  #message("\nFinal CTdf:")
  #message(paste(capture.output(CTdf), collapse = "\n"))
  return(CTdf)
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
              cat("\nEvaluating", model_name, "...\n")
              
              tryCatch({
                model <- models[[model_name]]$fit()
                
                assumption_results <- check_assumptions(model, model_type = models[[model_name]]$type)
                
                if (isTRUE(assumption_results)) {
                  cat("✓ Model passed all assumptions\n")
                  valid_models[model_name] <- AIC(model)
                  
                  cat("\nModel Summary:\n")
                  print(summary(model))
                  cat("AIC:", AIC(model), "\n")
                  
                  results[[model_name]] <- model
                } else {
                  cat("✗ Model failed assumptions:", paste(assumption_results, collapse = ", "), "\n")
                }
                
              }, error = function(e) {
                cat("✗ Error fitting model:", conditionMessage(e), "\n")
              })
            }
            
            # Handle results
            if (length(valid_models) > 0) {
              # Find best model
              best_model_name <- names(valid_models)[which.min(valid_models)]
              best_model <- results[[best_model_name]]
              
              # Extract p-values based on model type
              if (best_model_name == "Mixed Effects Model") {
                library(lmerTest)
                p_values <- summary(best_model)$coefficients[-1, "Pr(>|t|)"]
              } else {
                p_values <- summary(best_model)$coefficients[-1, 4]
              }
              
              cat("\nBest model:", best_model_name, "\n")
              
              return(list(
                model_name = best_model_name,
                p_values = p_values
              ))
            } else {
              cat("\nNo valid models found.\n")
              return(NULL)
            }
          }
          
          # Run the evaluation and store results
          model_results <- evaluate_models()
          
          # Store results in CTdf if a valid model was found
          if (!is.null(model_results)) {
            CTdf$pvalue[CTdf$gene == gene & CTdf$chemicalgroup == chemical] <- model_results$p_values
            CTdf$model_name[CTdf$gene == gene & CTdf$chemicalgroup == chemical] <- model_results$model_name
          }
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
  #CTdf <- apply_outlier_test(CTdf, stats_on)
  CTdf <- calculate_sem(CTdf, "^RQ")
  CTdf <- add_concentration_and_groups(CTdf)  
  CTdf <- calculate_pvalues(CTdf, stats_on, control_to_check)
  CTdf <- adjust_labels(CTdf, analysistype, gene_or_protein_label)
  # Adjust concentration display and move columns
  #CTdf$chemical <- paste0(format(CTdf$concentration, nsmall = 1), unit, " ", as.character(CTdf$chemicalgroup), " ", CTdf[[gene_or_protein_label]])
  #CTdf <- CTdf[, c("chemical", "concentration", "chemicalgroup", setdiff(names(CTdf), c("chemical", "concentration", "chemicalgroup")))]
  
  CTdf
}
#formatted_data <- read.csv("/home/ssyeddan/storage/Cursor_Atlas_Lab/Atlas_Lab_app/qPCR_formatted_data_20250108_115413.csv")
#formatted_data$chemical <- tolower(gsub(" ", "", formatted_data$chemical))
#final_df <- process_all_replicates(formatted_data, analysis_type, control_to_check)
#head(final_df)

#CTdf <- process_CTdf(final_df, stats_on = "RQ", control_to_check = control_to_check, 
#                    analysistype = analysis_type, unit = "uM", 
#                    gene_or_protein_label = "gene")
#write.table(CTdf, paste0("~/storage/Atlas_lab_app/Columbus_Analysis/", "resall.txt"), sep = "\t", row.names = FALSE)

