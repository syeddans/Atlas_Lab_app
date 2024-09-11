library(dplyr)
library("DescTools")
library("CATT")
library("coin")
library(plotrix)
library(multcomp)
library(stringr)
library("aod")
path = "~/storage/Columbus Processing/Nilesh Results2"
files = list.files(path, pattern=NULL, all.files=FALSE, 
                   full.names=FALSE)

control = "DMSO"

filenum = 0
for (file in files) { 
  filenum <- filenum +1 
  if (filenum  == 1){ 
    df <- read.table(paste0(path,"/", file),sep="\t",header=TRUE,fill = TRUE)
  }
  else { 
    dftemp <- read.table(paste0(path,"/", file),sep="\t",header=TRUE,fill = TRUE)
    df <- rbind(df,dftemp)
  }
  
}

#if concentration 
#df$chemical <- as.double(gsub("([0-9]*\\.?[0-9]+).*", "\\1", as.character(df$chemical)))

#if concentration is after chemical name 
df$chemical <- str_extract(df$chemical, "(?<=\\s)[0-9.]+")



names(df)[names(df) == "chemical"] <- "Dose"
df$Dose <- replace(df$Dose, is.na(df$Dose), 0)
df$Dose<-as.double(df$Dose)
as.character(df$Dose)
df$graphing_sort <- paste0(format(df$Dose, nsmall = 1), "_", format(df$chemicalgroup, nsmall = 1))



df3<- df %>% 
  group_by(graphing_sort, Dose, chemicalgroup) %>% 
  summarize("cellswithlipid"= mean(X..Cells.with.Lipid), SEM = std.error(X..Cells.with.Lipid))



for(chemicals in unique(df3$chemicalgroup)) { 
  df2 <- df[df$chemicalgroup==chemicals, ]
  df2<- rbind(df2,df[df$chemicalgroup==control, ])
  if(length(unique(df2$Dose))>1) { 
    df2$Dose<- factor(df2$Dose)
    model <- glm(X..Cells.with.Lipid ~ Dose, family = quasibinomial, data = df2)
    anova(model, test= "F")
    summary<-summary(glht(model, linfct = mcp(Dose = "Dunnett")))
    p_values <- summary$test$pvalues
    graphing_sort_values <- df3[df3$chemicalgroup == chemicals, "graphing_sort"]
  
    # Create a data frame with p-values and graphing_sort values as separate rows
    temp <- data.frame(pvalues = p_values, graphing_sort = graphing_sort_values)
    
    df3<- left_join(df3, temp, by="graphing_sort")
  }
  else if(chemicals != control&& length(unique(df2$Dose))==1) { 
    a<-aov(X..Cells.with.Lipid~chemicalgroup, data = df2)
    p_value <- summary(a)[[1]]$`Pr(>F)`[1]
    graphing_sort_values <- df3[df3$chemicalgroup == chemicals, "graphing_sort"]
    
    # Create a data frame with p-values and graphing_sort values as separate rows
    temp <- data.frame(pvalues = p_value, graphing_sort = graphing_sort_values)
    
    df3<- left_join(df3, temp, by="graphing_sort")
    
    }
  }
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
    return(" ")
  }
}
df3$chemicalgroup <- gsub("MI", "~MI", df3$chemicalgroup)
df3$chemicalgroup <- gsub(" P201", "P201", df3$chemicalgroup)
df3$graphing_sort <- gsub(" 0.00_MI", "~MI", df3$graphing_sort)
# Add a column to the dataframe with the significance codes
df3$sigcode <- sapply(df3$pvalue, get_sigcode)
df3$graphing_sort <- gsub(" ", "", df3$graphing_sort)
df3$Dose <- as.double(df3$Dose)

df3 <- df3[order(df3$chemicalgroup), ]

rows_to_move <- c(1:4)

# Subset the rows to move
rows_to_move_df <- df3[rows_to_move, ]

# Remove the rows to move from the original data frame
df3 <- df3[!(1:nrow(df3) %in% rows_to_move), ]

# Append the subsetted rows to the end of the original data frame
df3 <- rbind(df3, rows_to_move_df)

