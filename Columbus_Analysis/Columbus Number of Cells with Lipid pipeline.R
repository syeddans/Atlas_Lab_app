args <- commandArgs(trailingOnly = TRUE)
Columbus_files_dir <- args[1]
output_dir <- args[2]
wellmap_dir <- args[3]
data_analysis_type <- args[4]

library(dplyr)
library(forcats)
library("Rcpp")
library("matrixStats")
library(parallel)

#"~/storage/Atlas_lab_app/2023-10-10 BODIPY 3T3L1 bisphenols cellmask rep1-p2-june19-JC" #"~/storage/Atlas_lab_app/3d_testdata/3d_test_data/2024-04-03 liver spheroids BODIPY lean vs Tox"

#Reads in the wellmap provided checking for errors
read_well_map <- function(wellmap_dir_path, analysis_name) {
  wellmap_file <- paste0(analysis_name, "-wellmap.csv")
  wellmap_path <- file.path(wellmap_dir_path, wellmap_file)
  
  # Read the well map with error handling
  well_map <- tryCatch(
    read.csv(wellmap_path, sep = ",", header = FALSE, row.names = 1),
    error = function(e) {
      #TODO: I think because its a script ran in the shiny app I can just print this error message directly to the server
      stop("Error reading well map. Ensure headers and row names are complete.\n", e)
    }
  )
  
  # Set column names and remove header row
  colnames(well_map) <- well_map[1, ]
  well_map <- well_map[-1, ]
  
  return(well_map)
}

#Connects a lookup table that has the chemical with its Well based on the wellmap provided 
generate_lookup <- function(well_map) {
  # Create all combinations of rows and columns
  combinations <- character()
  
  for (row in rownames(well_map)) {
    for (col in colnames(well_map)) {
      combinations <- c(combinations, paste0(row, col))
    }
  }
  lookup <- data_frame(WellName =combinations) 
  chemical <- character() 
  
  for (row in seq_along(rownames(well_map))) { 
    for (col in seq_along(colnames(well_map))){ 
      chemical <- c(chemical, paste0(well_map[row,col])) 
    }
  }
  lookup$Chemical <- chemical
  lookup[lookup == "NA"] <- NA
  lookup <- na.omit(lookup)
  
  return(lookup)
}

extract_well_name <- function(filename) {
  # Use your regex to extract well names
  gsub(".*\\.([A-Z0-9]+)\\[.*", "\\1", basename(filename))
}

generate_cellID <- function(length = 8) {
  paste0(sample(c(0:9, letters, LETTERS), length, replace = TRUE), collapse = "")
}

#calc_distances 
cppFunction('
        NumericMatrix calc_distances(NumericMatrix current, NumericMatrix above) {
          int n_current = current.nrow();
          int n_above = above.nrow();
          NumericMatrix distances(n_current, n_above);
          
          for (int i = 0; i < n_current; i++) {
              for (int j = 0; j < n_above; j++) {
                  // Using the explicit formula for clarity
                  distances(i, j) = sqrt(pow(above(j, 0) - current(i, 0), 2) +pow(above(j, 1) - current(i, 1), 2)  // euclidean_distance =  sqrt((x2 - x1)^2 + (y2 - y1)^2)
                  );
              }
          }
          
          return distances;
      }')

# Get closest cells
get_closest_cell <- function(column) {
  min_distance <- min(column)  
  closest_index <- which.min(column)  
  return(c(min_distance, closest_index))  
}


process_file_through_planes <- function(table,well, threshold_distance = 5) {
  
  table$CellID <- sapply(1:nrow(table), function(x) generate_cellID(8))
  #table <- table[table$Plane <5,]
  # Split table by Plane
  plane_tables <- split(table, table$Plane)
  
  # Process distances between planes
  for (i in seq_along(plane_tables)[-length(plane_tables)]) {
    current_plane_table <- plane_tables[[i]]       
    above_plane_table <- plane_tables[[i + 1]]     
    
    cell_file_condition = "Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm." %in% colnames(current_plane_table)
    lipid_file_condition = "Lipid...Corrected.Spot.Intensity" %in% colnames(current_plane_table)
    # Convert coordinates to a matrix
    if (cell_file_condition) { 
      current_coords <- as.matrix(current_plane_table[, c("Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm.", 
                                                          "Nuclei.Selected...Cytoplasm.Centroid.Y.in.Image..µm.")])
      above_coords <- as.matrix(above_plane_table[, c("Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm.", 
                                                      "Nuclei.Selected...Cytoplasm.Centroid.Y.in.Image..µm.")])
    } else if (lipid_file_condition) { 
      current_coords <- as.matrix(current_plane_table[, c("Position.X..µm.", 
                                                          "Position.Y..µm.")])
      above_coords <- as.matrix(above_plane_table[, c("Position.X..µm.", 
                                                      "Position.Y..µm.")])
    }
    
    # Calculate distances using the C++ function
    distances <- calc_distances(current_coords, above_coords)
    rownames(distances) <- current_plane_table$CellID
    colnames(distances) <- above_plane_table$CellID
    
    #TODO: I thought this would be a dependent process and mclapply would not help, but testing it showed that it did, im not sure why 
    results <- mclapply(1:ncol(distances), function(j) {
      get_closest_cell(distances[, j])
    }, mc.cores = parallel::detectCores() - 1)
    
    results_matrix <- do.call(rbind, results)
    min_distances <- results_matrix[, 1]  
    closest_cells <- results_matrix[, 2]
    
    # Update CellID where closest distance is below threshold
    above_plane_table$CellID[min_distances <= threshold_distance] <- 
      rownames(distances)[closest_cells[min_distances <= threshold_distance]]
    
    # Update the plane table
    plane_tables[[i + 1]] <- above_plane_table
  }
  
  # Combine processed data
  processed_table <- bind_rows(plane_tables)
  
  # Summarize results
  if (cell_file_condition) {
    result <- processed_table %>%
      group_by(CellID) %>%
      summarize(
        Plane = paste(min(Plane), max(Plane), sep = "-"),
        Num_of_Lipids = sum(Nuclei.Selected...Cytoplasm.Overlap...., na.rm = TRUE),
        WellName = well,
        .groups = "drop"
      )
  } else if (lipid_file_condition) {
    result <- processed_table %>%
      group_by(CellID) %>%
      summarize(
        Plane = paste(min(Plane), max(Plane), sep = "-"),
        Avg_Lipid_Intensity = mean(Lipid...Corrected.Spot.Intensity, na.rm = TRUE),
        Avg_Lipid_Area = mean(Lipid...Spot.Area..px.., na.rm = TRUE),
        WellName = well,
        .groups = "drop"
      )
  }
  
  
}
read_columbus_files <- function(file, lookup,  Num_lipid_threshold = 15) {
  
  possible_columns <- c("WellName", 
                        "Plane", 
                        "Nuclei.Selected...Cytoplasm.Overlap....", 
                        "Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm.", 
                        "Nuclei.Selected...Cytoplasm.Centroid.Y.in.Image..µm.",
                        "Lipid...Spot.Area..px..",
                        "Lipid...Corrected.Spot.Intensity",
                        "Position.X..µm.", 
                        "Position.Y..µm."
  )  
  
  well <- extract_well_name(file)
  print(well)
  
  # Skip if well is not in lookup
  if (!(well %in% as.vector(lookup$WellName))) {
    print("skipped")
    return(NULL)  # Return NULL for skipped files
  }
  
  # Read file if not empty
  if (file.size(file) > 0) { 
    table <- read.table(file, sep = "\t", header = TRUE, fill = TRUE)
    
    #grab columns of interest 
    selected_columns <- possible_columns[possible_columns %in% colnames(table)]
    # Create a new data frame with selected columns, filling with NA where necessary
    table <- table[, selected_columns, drop = FALSE]  # Use drop = FALSE to prevent dropping to vector
    
    cell_table <- process_file_through_planes(table = table, well= well)
    
    if("Num_of_Lipids" %in% colnames(cell_table)) { 
      cell_table_summarized <- cell_table %>%
        group_by(WellName) %>%
        summarise(
          Number_of_Cells = n(),
          "Cells_with_Lipid" = sum(Num_of_Lipids > Num_lipid_threshold),
          "Lipid_Count" = sum(Num_of_Lipids)
        ) %>% 
        mutate( 
          "% Cells with Lipid" = `Cells_with_Lipid`/`Number_of_Cells`,
          "Num_Lipid_per_Cell" = `Lipid_Count`/`Number_of_Cells`,
        )
    } 
    else if("Avg_Lipid_Intensity" %in% colnames(cell_table)) { 
      cell_table_summarized <- cell_table %>%
        group_by(WellName) %>%
        summarise(
          Number_of_Lipid = n(),
          "Lipid_Intensity" = sum(Avg_Lipid_Intensity),
          "Lipid_Area" = sum(Avg_Lipid_Area)
        ) %>% 
        mutate( 
          "Avg_Lipid_Intensity" = `Lipid_Intensity`/`Number_of_Lipid`,
          "Avg_Lipid_Area" = `Lipid_Area`/`Number_of_Lipid`,
        ) 
    }
    result = cell_table_summarized
    
    
    
    
    return(result)  
  }
  
  return(NULL)  # Return NULL for empty files
}

analysis_name = basename(Columbus_files_dir)
w <- read_well_map(wellmap_dir, analysis_name)
l <- generate_lookup(w)
files <- list.files(Columbus_files_dir, full.names = TRUE)
results <- lapply(files, read_columbus_files, lookup = l)
combined_results <- bind_rows(results) 
combined_results_with_chemicals <- left_join(combined_results, l, by = "WellName")
combined_results_with_chemicals$chemicalgroup <- gsub("^.*? ", "", combined_results_with_chemicals$Chemical) 

# Formatting of Data
# concentration after 
#datamerged$chemicalgroup <- gsub("\\s.*", "", datamerged$chemical)

#This is to merge the two file type of data together into one row based on WellName
concat_unique <- function(x){paste(unique(x[!is.na(x)]),  collapse=',')}
combined_results_with_chemicals_merged <- combined_results_with_chemicals %>%
  group_by(WellName) %>%
  summarize(across(everything(), concat_unique), .groups = "drop")

outputFileName = paste0(analysis_name, ".txt") 
write.table(combined_results_with_chemicals_merged, paste0(output_dir, outputFileName), sep="\t", row.names=FALSE)





