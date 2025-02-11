# qPCR Data Analysis Pipeline

This document explains how to prepare and analyze qPCR data using the Atlas Lab analysis tools.

## File Requirements

For each qPCR plate, you need 4 CSV files:
1. CT values file - Contains the raw CT values from your qPCR run
2. Genes tested file - Lists the genes analyzed in each well
3. Chemical treatments file - Specifies the treatment for each well
4. Replicate information file - Indicates which wells are replicates

## Data Organization

### Example Data Structure
You can download example data <a id="download_example1" href="#" class="shiny-download-link">Example1_qPCR_input_data.zip</a>, <a id="download_example2" href="#" class="shiny-download-link">Example2_qPCR_input_data.zip</a> to see the required format.

### File Organization Options
1. **Folder of Folders**: Create a main folder (make sure to compress to .zip if uploading a folder) containing subfolders for each plate, with the 4 required CSV files in each subfolder
2. **Individual Upload**: Upload the 4 required files for a single plate directly

## Data Processing Pipeline

1. **qPCR Formatter**
   - Takes your raw data files
   - Combines and summarizes the data into a single CSV file
   - This formatted file feeds into the analysis pipeline
Here's an example of the results format: <a id="download_formatted_results1" href="#" class="shiny-download-link">Example_Formatted_Results1.csv</a>, <a id="download_formatted_results2" href="#" class="shiny-download-link">Example_Formatted_Results2.csv</a>

2. **qPCR Analysis**
   - Calculates RQ values
   - Performs statistical analysis
   - Generates interactive graphs and tables
   - Produces a complete data summary webpage
<a id="download_final_results1" href="#" class="shiny-download-link">Example_Results_Final1.zip</a>, <a id="download_final_results2" href="#" class="shiny-download-link">Example_Results_Final2.zip</a>

## Stats
- Statistical analysis is done on RQ values using a generalized linear model (GLM), all model families are tested and the best model is selected based on AIC Values for each gene for each treatment to gather p-values. 

### Important Considerations

#### Naming Consistency
- Be consistent with chemical names (e.g., use either "Phenanthrene" or "Phe" consistently, not both)
- Include concentrations in treatment names when multiple concentrations are tested (e.g., "100uM Phe")
- Treatment names don't require concentrations for single-concentration experiments

### Control Selection
- When selecting control gene (e.g., β-actin) or control treatment:
  - Must match exactly as written in your files
  - Be consistent with naming (e.g., "beta-actin" vs "b-actin")
  - For treatments, be precise (e.g., "MIE" vs "MI ethanol")

## Output Files

The analysis pipeline generates:
1. Interactive web summary with graphs and tables
2. CSV files showing data at each analysis step
3. Log file documenting statistical tests and analysis parameters


## Individual stats 
You may also chose to use your own RQ values already generated and run stats on those. You may upload them in this format <a id="download_stats_format" href="#" class="shiny-download-link">Stats_Only_Format.csv</a> to get stats results. Note: Every new rep should have RQ+rep number as the header, for example RQ1, RQ2, RQ3 etc. 

For questions or issues, please contact the Atlas Lab team.

