# Author: Angela Sofia Garcia
# For Agrosavia
# Master researcher call A-772020 – Technical Test
# Code for Exercise 2 - parsing data

# First load the files
ariba <- read.csv("ariba_amr_output.csv", header = T)
metadata <- read.csv("ncbi_acquired_genes_metadata.csv", header = T)

## To find only those clusters genes of data for CEPHALOSPORIN and CARBAPENEM
## ID column of metadata to use for gene name
summary(is.na(metadata$allele)) # Too many NA
summary(is.na(metadata$gene_family)) # Only 13 NA, column to use

## Extract all gene names per antimicrobial subclass
g_cephalosporin <- unique(as.vector(metadata[(metadata$subclass == "CEPHALOSPORIN"), 'gene_family']))
g_carbapenem <- unique(as.vector(metadata[(metadata$subclass == "CARBAPENEM"), 'gene_family']))
g_ceph_cabap <- unique(as.vector(metadata[(metadata$subclass == "CEPHALOSPORIN" | metadata$subclass == "CARBAPENEM"), 'gene_family']))

# Note: There are shared gene family names between those two subclasses

# Subset data-frame with only the genes of interest by their prefix
new <- data.frame(matrix(ncol = 0, nrow = nrow(ariba)))
new <- cbind(new, name= ariba$name)
myGenes <- c()

for (i in g_ceph_cabap) {
  temp <- ariba[,grepl(paste0("^", i), colnames(ariba))]
  if (ncol(temp)!=0) {myGenes <- c(myGenes, i)}
  new <- cbind(new, temp)
}


# Get only present gen names and positions in my new data frame
gen_assembled_names <- colnames(new)[c(endsWith(colnames(new), ".assembled"))]
gen_coverage_names <- colnames(new)[c(endsWith(colnames(new), ".ctg_cov"))]
gen_assembled_col <- which(colnames(new) %in% gen_assembled_names)

samples <- data.frame(matrix(ncol = 4, nrow = 0))

for (j in 1:length(gen_assembled_col)) {
  for (k in nrow(new)) {
    # conditional to be present in the allele
    if (new[k,gen_assembled_names[j]]!="no") {
      # conditional to have at least 10 coverage
      if (new[k,gen_coverage_names[j]]>10) {
        samples <- rbind(samples, data.frame(sampleID=new[k,'name'],
                                             gene=gsub( ".assembled*", " ", gen_assembled_names[j]),
                                             allele=new[k,gen_assembled_names[j]],
                                             coverage=new[k,gen_coverage_names[j]]))
      }
    }
  }
}

write.csv(samples, file = "samples_antimicrobial_genes.csv")
