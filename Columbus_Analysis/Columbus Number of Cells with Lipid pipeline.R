#required arguments
args <- commandArgs(trailingOnly = TRUE)
Columbus_files_dir <- args[1]
output_dir <- args[2]
wellmap_dir <- args[3]

library(dplyr)
library(forcats)
library("Rcpp")
library("matrixStats")
library(parallel)

path = Columbus_files_dir #"~/storage/Atlas_lab_app/2023-10-10 BODIPY 3T3L1 bisphenols cellmask rep1-p2-june19-JC" #"~/storage/Atlas_lab_app/3d_testdata/3d_test_data/2024-04-03 liver spheroids BODIPY lean vs Tox"
AnalysisName = basename(path) 
path2 = paste0(path,"/")
files = list.files(path, pattern=NULL, all.files=FALSE, 
           full.names=FALSE)

# Map the chemical to its well (moved to be processed first so program can crash rigth away if there is problem instead of waiting processing and wasting time)
wellMapFileName = paste0(AnalysisName, "-wellmap.csv") #wellmap file name should come in as analysis file name + -wellmap
#TODO: add error catching: if it says "row.names cannot be duplicated" it likely means the user forgot to fill out the wellmap headers fully like forgetting P or 12. 
wellMap <- read.csv(paste0(wellmap_dir, "/", wellMapFileName),sep =",", header = FALSE, row.names = 1)#read.csv(paste0("~/storage/Atlas_lab_app/2023-10-10 BODIPY 3T3L1 bisphenols cellmask rep1-p2-june19-JC-wellmap", "/", wellMapFileName),sep =",", header = FALSE, row.names = 1)

colnames(wellMap) <- wellMap[1, ]
wellMap <- wellMap[-1, ]
combinations <- character()

for (row in rownames(wellMap)) {
  for (col in colnames(wellMap)) {
    combinations <- c(combinations, paste0(row, col))
  }
}
lookup <- data_frame(WellName =combinations) 
chemical <- character() 

for (row in seq_along(rownames(wellMap))) { 
  for (col in seq_along(colnames(wellMap))){ 
    chemical <- c(chemical, paste0(wellMap[row,col])) 
  }
}
lookup$chemical <- chemical

#remove NA rows so that we do not need to process them later if they are not indicated in the wellmap
lookup[lookup == "NA"] <- NA
lookup <- na.omit(lookup)


generate_cellID <- function(length = 8) {
  paste0(sample(c(0:9, letters, LETTERS), length, replace = TRUE), collapse = "")
}

# first tried the above function in R but was too slow when using lapply, so opted to use c++ code for faster processing of euclidean distance which is a pairwise calculation
euclidean_distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}
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


full_table <- data.frame()
# Grab data from files
for (file in files) { 
  well <- gsub(".*\\.([A-Z0-9]+)\\[.*", "\\1", file)
  print(well)
  file_name <- paste0(path2,file)
  
  #if chemical not given in wellmap, skip its processing
  if(!(well %in% as.vector(lookup$WellName))) {
    print("skipped")
    next
  }
  
  #check to make sure the file isnt empty, the well could have no cells and thus would provide an empty file even though it was processed
  if(file.size(file_name)>0){ 
    table <- read.table(file_name,sep="\t",header=TRUE,fill = TRUE)
    #print(colnames(table))
    table <- table[,c("WellName","Plane", "Nuclei.Selected...Cytoplasm.Overlap....", "Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm.", "Nuclei.Selected...Cytoplasm.Centroid.Y.in.Image..µm.")]
    num_of_cells <- nrow(table)
    # Threshold of 5 lipids as buffer for analytical error
    # table <- table[table$`Nuclei.Selected...Cytoplasm.Overlap....` > 5, ]
    table$CellID <- sapply(1:nrow(table), function(x) generate_cellID(8))
    
    plane_tables <- split(table, table$Plane)
    #TODO: calculating through each plane increases time linearly, find a way to speed this up, cant use mclapply because plane calculations are dependent on one another 
    for (i in seq_along(plane_tables)[-length(plane_tables)]) {  # Skip the last plane as it has no plane above it
      # Define the C++ function
      
      # Inside your main function or processing loop
      current_plane_table <- plane_tables[[i]]       
      above_plane_table <- plane_tables[[i + 1]]     
      
      # Convert the relevant columns to a numeric matrix
      current_coords <- as.matrix(current_plane_table[, c("Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm.", "Nuclei.Selected...Cytoplasm.Centroid.Y.in.Image..µm.")])
      above_coords <- as.matrix(above_plane_table[, c("Nuclei.Selected...Cytoplasm.Centroid.X.in.Image..µm.", "Nuclei.Selected...Cytoplasm.Centroid.Y.in.Image..µm.")])
     
      #start <- Sys.time()
       # Calculate distances using the C++ function
      distances <- calc_distances(current_coords, above_coords)
      #print( Sys.time() - start )
       # Assign row and column names to the distance matrix based on CellID
      rownames(distances) <- current_plane_table$CellID
      colnames(distances) <- above_plane_table$CellID
      
      # Set a threshold for determining if cells are "close enough"
      threshold_distance <- 5
      
      get_closest_cell <- function(column, distances) {
        min_distance <- min(column)  # Get the minimum distance in the column
        closest_index <- which.min(column)  # Get the index of the closest cell
        return(c(min_distance, closest_index))  # Return both values as a vector
      }

      # Use mclapply to compute min distances and closest cells in parallel
      results <- mclapply(1:ncol(distances), function(j) {
        get_closest_cell(distances[, j], distances)  # Process each column
      }, mc.cores = parallel::detectCores() - 1)  # Use one less than the available cores
      
      # Combine the results into a matrix
      results_matrix <- do.call(rbind, results)
      min_distances <- results_matrix[, 1]  # Extract the first column for minimum distances
      closest_cells <- results_matrix[, 2]
      
      # Update CellID where closest distance is below threshold
      above_plane_table$CellID[min_distances <= threshold_distance] <- 
        rownames(distances)[closest_cells[min_distances <= threshold_distance]]
      
      # Update the plane table
      plane_tables[[i + 1]] <- above_plane_table
    }
    
    processed_table <- bind_rows(plane_tables)
    result <- processed_table %>%
      group_by(CellID) %>%                            # Group by cellID
      summarize(
        Plane = paste(min(Plane), max(Plane), sep = "-"),  # List Plane "min-max" as start-stop
        Num_of_Lipids = sum(Nuclei.Selected...Cytoplasm.Overlap....), # Sum number of lipids in cell accross all planes
        WellName = well
      ) %>%
      ungroup()
  }
  else{next}
  
  full_table <-bind_rows(full_table,result)
}

summary_table <- full_table %>% 
  group_by(WellName) %>%
  summarize(
    "Number of Cells"= n(),
    "Cells with Lipid" = sum(Num_of_Lipids >10),
  )
# Formatting of Data
summary_table$`Number of Cells` <- as.numeric(summary_table$`Number of Cells`)
summary_table$`Cells with Lipid` <- as.numeric(summary_table$`Cells with Lipid`)
summary_table$"% Cells with Lipid" <- summary_table$`Cells with Lipid`/summary_table$`Number of Cells`

data = summary_table
datamerged <- left_join(data, lookup, by = "WellName")
datamerged$chemicalgroup <- gsub("^.*? ", "", datamerged$chemical) 

# concentration after 
#datamerged$chemicalgroup <- gsub("\\s.*", "", datamerged$chemical)

outputFileName = paste0(AnalysisName, ".txt") 
write.table(datamerged, paste0(output_dir, outputFileName), sep="\t", row.names=FALSE)



