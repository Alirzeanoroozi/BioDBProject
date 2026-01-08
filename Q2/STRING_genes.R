library(readxl)  
library(dplyr)    
setwd("~/Desktop")
#To load data
string_data <- read_excel("INTERACTIONS.xlsx") 
#Keep columns with the genes
string_genes <- string_data[, 1:2]

#Combining columns into one vector
all_genes <- c(string_genes[[1]], string_genes[[2]])
all_genes_df <- data.frame(Gene = all_genes)

#To view enteries
head(all_genes_df)
#Checking total number of genes in combined vector
length(all_genes)

#Count unique genes
length(unique(all_genes))

# Remove duplicates
unique_genes <- unique(all_genes)
gene_vector <- as.vector(unique_genes)


head(gene_vector)
length(gene_vector)

# Export to txt file
write.table(gene_vector,
            file = "~/Desktop/STRING_gene_list.txt",  
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
