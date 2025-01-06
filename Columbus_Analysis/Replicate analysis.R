library(dplyr)
library("DescTools")
library("CATT")
library("coin")
library(plotrix)
library(multcomp)
library(stringr)
library("aod")

# test path 
path = results_dir #"~/storage/Atlas_lab_app/test_output/"
files = list.files(path, pattern=NULL, all.files=FALSE, 
                   full.names=FALSE)

#default "DMSO" MAKE SURE THE CONTROL IS RIGHT OR ELSE IT WILL BREAK DUE TO STATS COMPARISONS 
control_dose_num = control_dose
control = control_name

print(files)
# Extract Data from files created from script 
data_analysis_type <-input$data_analysis_type
df <- data_frame()
for (file in files) { 
    dftemp <- read.table(paste0(path,"/", file),sep="\t",header=TRUE,fill = TRUE)
    dftemp<- dftemp %>%
      group_by(Chemical) %>% 
      reframe("Metric_full"= case_when( 
                data_analysis_type == "1" ~ X..Cells.with.Lipid),
                #data_analysis_type == "2" ~ Number_of_Lipid,
                #data_analysis_type == "3" ~ Avg_Lipid_Intensity, 
                #data_analysis_type == "4" ~ Avg_Lipid_Area), #technical rep average
                "chemicalgroup"= chemicalgroup,
                "Well"= WellName,
                "SEM" = std.error(X..Cells.with.Lipid))
    df <- rbind(df,dftemp)
}

#if concentration before chemical name
df$Chemical <- as.double(gsub("([0-9]*\\.?[0-9]+).*", "\\1", as.character(df$Chemical)))

#if concentration is after chemical name 
#df$chemical <- str_extract(df$chemical, "(?<=\\s)[0-9.]+")


# Formatting for graphing later on
names(df)[names(df) == "Chemical"] <- "Dose"
df$Dose <- replace(df$Dose, is.na(df$Dose), 0)
df$Dose<-as.double(df$Dose)
df$graphing_sort <- paste0(format(df$Dose, nsmall = 1), "_", format(df$chemicalgroup, nsmall = 1))
df3<- df %>% 
  group_by(graphing_sort, Dose, chemicalgroup) %>% 
  summarize("Metric"= mean(Metric_full), SEM = std.error(Metric_full))

browser() 
null_stats = TRUE
if(length(files)>1){
  null_stats = FALSE
    # Statisical Analysis
    # TODO: add other statisical tests
    # TODO: right now stats are assuming technical replications are biological replicates
    for(chemicals in unique(df3$chemicalgroup)) { 
      if (chemicals == control) {
        next 
      }
      # Have to compare a group of chemicals with their concentrations to control so put them together
      df2 <- df[df$chemicalgroup==chemicals, ]
      #if(control_dose_num !=0) {
      #  df2<- rbind(df2,df[(df$chemicalgroup==control&df$Dose==control_dose_num), ])
      #}
      #else{ 
        df2<- rbind(df2,df[df$chemicalgroup==control, ])
      #  }
      
      
      #TODO:for now set control dose to 0 as I figure out how to compare dose controls
      df2$Dose[df2$chemicalgroup==control] <- -1 
      df2 <- df2[order(df2$Dose), ]
      df2$Dose <- factor(df2$Dose, levels = sort(unique(df2$Dose)))
      
      # This is for comparing multiple concentrations using a quasibinomial distribution model 
      #browser()
        model <- glm(Metric_full ~ Dose, family = quasibinomial, data = df2)
        #anova(model, test= "F")
        model_summary <- summary(model)
        p_values <- model_summary$coefficients[-1, 4]
        graphing_sort_values <- df3[df3$chemicalgroup == chemicals, "graphing_sort"]
       
        #browser()
        # Create a data frame with p-values and graphing_sort values as separate rows
        temp <- data.frame(pvalues = p_values, graphing_sort = graphing_sort_values)
        df3<- left_join(df3, temp, by="graphing_sort")
      
      
      # This is for if there is only one concentration for example MID, to do a basic anova 
      # else if(chemicals != control&& length(unique(df2$Dose))==1) { 
      #   a<-aov(Metric~chemicalgroup, data = df2)
      #   p_value <- summary(a)[[1]]$`Pr(>F)`[1]
      #   graphing_sort_values <- df3[df3$chemicalgroup == chemicals, "graphing_sort"]
      #   
      #   # Create a data frame with p-values and graphing_sort values as separate rows
      #   temp <- data.frame(pvalues = p_value, graphing_sort = graphing_sort_values)
      #   df3<- left_join(df3, temp, by="graphing_sort")
      # }
    }
  }


# p value indicators based on analysis 
if(null_stats == FALSE){ 
  pvalue_columns <- names(df3)[startsWith(names(df3), "pvalues")]
  df3 <- df3 %>% 
    mutate(pvalue = coalesce(!!!syms(pvalue_columns)))
  df3 <- df3[, !names(df3) %in% pvalue_columns]
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
      return("")
    }
  }
  # Add a column to the dataframe with the significance codes
  df3$sigcode <- sapply(df3$pvalue, get_sigcode)
}else { 
  df3$sigcode <- ""
  }
# This is needed to bring the control to the end if needed 
# TODO: Figure out a more elegant way to do this 
#df3$chemicalgroup <- gsub("MI", "~MI", df3$chemicalgroup)
#df3$graphing_sort <- gsub(" 0.00_MI", "~MI", df3$graphing_sort)


df3$graphing_sort <- gsub(" ", "", df3$graphing_sort)
df3 <- df3 %>%
  arrange(chemicalgroup, Dose)

# Factor levels for chemicalgroup and Dose
df3$chemicalgroup <- factor(df3$chemicalgroup, levels = unique(df3$chemicalgroup))
df3$Dose <- factor(df3$Dose, levels = unique(df$Dose))
df <- df[order(df$chemicalgroup, df$Dose), ]

