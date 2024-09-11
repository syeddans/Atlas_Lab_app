#required arguments
args <- commandArgs(trailingOnly = TRUE)
Columbus_files_dir <- args[1]
output_dir <- args[2]
wellmap_dir <- args[3]

library(dplyr)
library(forcats)

AnalysisName = basename(Columbus_files_dir) 
path = Columbus_files_dir
path2 = paste0(path,"/")
files = list.files(path, pattern=NULL, all.files=FALSE, 
           full.names=FALSE)
plane = 1
db = data.frame()
db2 = data.frame() 

# Grab data from files
for (file in files) { 
  well <- gsub(".*\\.([A-Z0-9]+)\\[.*", "\\1", file)
  print(well)
  file_name <- paste0(path2,file)
  table <- read.table(file_name,sep="\t",header=TRUE,fill = TRUE)
  table <- table[,c(7,11,13,19)]
  #colnames(table)
  if (plane >0){ 
    table <-table[table$Plane == plane, ]
  }
  
  # Threshold of 10 lipids as buffer for overlapping and analytical error
  nonZeroRows <- subset(db2, table$Nuclei...Lipid.In.Cell..Overlap.... > 10)
  
  # Get the number of rows in the subset which is the number of cells that have lipid in it
  numNonZeroRows <- nrow(nonZeroRows)
  db2 <- rbind(db2, c(nrow(table), numNonZeroRows, well))
  db <- rbind(db,table)
}

# Formatting of Data
colnames(db) <- colnames(table)
colnames(db2) <- c("Number of Cells", "Cells with Lipid", "WellName")
db2$`Number of Cells` <- as.numeric(db2$`Number of Cells`)
db2$`Cells with Lipid` <- as.numeric(db2$`Cells with Lipid`)
sum(db2$`Number of Cells`) 
db2$"% Cells with Lipid" <- db2$`Cells with Lipid`/db2$`Number of Cells`

# Map the chemical to its well 
wellMapFileName = paste0(AnalysisName, "-wellmap.csv") #wellmap file name should come in as analysis file name + -wellmap
wellMap <- read.csv(paste0(wellmap_dir, "/", wellMapFileName),sep =",", header = FALSE, row.names = 1)
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


data = db2
datamerged <- left_join(data, lookup, by = "WellName")
datamerged$chemicalgroup <- gsub("^.*? ", "", datamerged$chemical) 

# concentration after 
#datamerged$chemicalgroup <- gsub("\\s.*", "", datamerged$chemical)

outputFileName = paste0(AnalysisName, ".txt") 
write.table(datamerged, paste0(output_dir, outputFileName), sep="\t", row.names=FALSE)



