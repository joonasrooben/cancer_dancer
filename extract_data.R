#setwd("C:/Users/Pauline Baur/Documents/STUDIUM/Roma/Digital_Epidemiology/Run2")


## I think you can start running this code in the ANALYSIS part (and just load 
## data from the zip folder i put on github)


### Load packages ###
library(TCGAbiolinks)
library(SummarizedExperiment)


## Code, if you want to isntall the the TCGAbiolinks and SummarizedExperiment
## (I was confused at first, because you can't get them from CRAN)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("SummarizedExperiment")


## Name of our cancer data
proj <- "TCGA-KIRC"

## Creates Folder to Store the Data in
dir.create(file.path(proj))

## Run1 - get Count Data
rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", 
                        sample.type = "Primary Tumor")
## Here I changed "Changed Primary Solid Tumor" to "Primary Tumor"

GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C) #doi:10.1038/nature12222
rna.expr.data.C <- assay(rna.data.C)
rna.genes.info.C <- rowRanges(rna.data.C)
rna.sample.info.C <- colData(rna.data.C)
write.table(rna.expr.data.C, 
            file=file.path(proj, paste(proj, "_rna_expr_data_C.txt",sep = "")), 
            row.names=TRUE, col.names=TRUE, quote = FALSE)
write.table(rna.sample.info.C@listData$patient, 
            file=file.path(proj,paste(proj, "_rna_patients_C.txt",sep = "")), 
            row.names=FALSE, col.names=FALSE, quote = FALSE)
write.table(rna.genes.info.C@ranges@NAMES, 
            file=file.path(proj,paste(proj, "_rna_genes_C.txt",sep = "")), 
            row.names=FALSE, col.names=FALSE, quote = FALSE)

save(rna.sample.info.C, file = "fullInfoC.RData")


rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "HTSeq - Counts", sample.type = "Solid Tissue Normal")
GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N)
rna.expr.data.N <- assay(rna.data.N)
rna.genes.info.N <- rowRanges(rna.data.N)
rna.sample.info.N <- colData(rna.data.N)
write.table(rna.expr.data.N, file=file.path(proj,paste(proj, "_rna_expr_data_N.txt",sep="")), row.names=TRUE, col.names=TRUE, quote = FALSE)
write.table(rna.sample.info.N@listData$patient, file=file.path(proj,paste(proj, "_rna_patients_N.txt",sep = "")), row.names=FALSE, col.names=FALSE, quote = FALSE)
write.table(rna.genes.info.N@ranges@NAMES, file=file.path(proj,paste(proj, "_rna_genes_N.txt",sep = "")), row.names=FALSE, col.names=FALSE, quote = FALSE)

#I modified the next couple lines of code, because the original one did not work
clinical.query<-GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
#write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.csv",sep="")), 
#          row.names = FALSE, quote = FALSE)
#changed sep = ";"
write.table(clinical.query, 
            file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), 
          row.names = FALSE, quote = FALSE, sep = ";")

save(rna.sample.info.N, file = "fullInfoN.RData")

## New Code 

load("fullInfoN.RData")
load("fullInfoC.RData")


#Create Summarized Experiments with the files extracted with the professor's code

SumExpN <- SummarizedExperiment(assays=rna.expr.data.N,
                     rowData=NULL, 
                     rowRanges=rna.genes.info.N,
                     colData=rna.sample.info.N,
                     metadata=list(),
                     checkDimnames=TRUE)

SumExpC <- SummarizedExperiment(assays=rna.expr.data.C,
                                rowData=NULL, 
                                rowRanges=rna.genes.info.C,
                                colData=rna.sample.info.C,
                                metadata=list(),
                                checkDimnames=TRUE)
## Save the summarized experiment files
save(list = c("SumExpN", "SumExpC"), file = "SumExp.RData")


### ANALYSIS ###

## load the saved .txt files
setwd("C:/Users/Pauline Baur/Documents/STUDIUM/Roma/Digital_Epidemiology/Run2/TCGA-KIRC")

rna_expr_data_C <- read.table("TCGA-KIRC_rna_expr_data_C.txt", header = TRUE)
rna_expr_data_N <- read.table("TCGA-KIRC_rna_expr_data_N.txt") #all different colnames

## these two data frames have the same rownames
all(rownames(rna_expr_data_C) == rownames(rna_expr_data_N))

## this object just contains the rownames of the two data frames above
rna_genes_C <- read.table("TCGA-KIRC_rna_genes_C.txt")
all(rownames(rna_expr_data_C) == rna_genes_C )

## same thing again (for normal cells)
rna_genes_N <- read.table("TCGA-KIRC_rna_genes_N.txt")
all(rna_genes_N == rna_genes_C )

## Also load the files containing a list of all the patient IDs
## Save them as a vector instead of a one-column data frame
rna_genes_C <- read.table("TCGA-KIRC_rna_patients_C.txt")
rna_genes_C  <- rna_genes_C[,1]
rna_genes_N <- read.table("TCGA-KIRC_rna_patients_N.txt")
rna_genes_N  <- rna_genes_N[,1]

## Load the clinical data
clinical_data <- read.table("TCGA-KIRC_clinical_data.txt", na = "NA",
                            heade = TRUE, sep = ";")



## I noticed, that some patients appear in the study multiple times
## but they have different values each time
doppelt <- which(duplicated(rna_genes_C))
#view an extract of the data of the people appearing multiple times
rna_expr_data_C[1:5, doppelt]



## Extract patients, that have normal an cancer tissue information

## Get a logical vector, indicating for each patient with cancer tissue,
## if normal tissue is available as well


both <- sapply(rna_genes_C, function(x) x %in% rna_genes_N)
sum(both)

## 72 people have both cancer and normal tissue

## use full rna_expr_data_N, but only partial rna_expr_data_C
rna_expr_data_C_part <- rna_expr_data_C[,both]
## this contains only patients with cancer and normal tissue
## this dataframe now has 72 columns

## modify colnames of our two data frames in a way that only the part
## identifying the person (e.g "TCGA.A3.3358") remains, in order to be able
## to compara the data frames


colnames(rna_expr_data_C_part) <- sapply(colnames(rna_expr_data_C_part), function(x)
gsub("\\.[0-9]{2}.\\.[0-9]{2}.\\.[0-9]{4}\\.[0-9]{2}", "",x))

colnames(rna_expr_data_N) <- sapply(colnames(rna_expr_data_N), function(x)
  gsub("\\.[0-9]{2}.\\.[0-9]{2}.\\.[0-9]{4}\\.[0-9]{2}", "",x))

## reorder rna_expr_data_C_part to match the order of rna_expr_data_N
Match <- sapply(colnames(rna_expr_data_C_part), 
                function(x)which(x == colnames(rna_expr_data_N)))


## reorder the columns in rna_expr_data_N
rna_expr_data_N <- rna_expr_data_N[, Match]
all(colnames(rna_expr_data_C_part) == colnames(rna_expr_data_N))
## it worked

## Merge the two data frames with count data of normal and cancer tissue into one
## data frame


## merge rna_expr_data_C_part and rna_expr_data_N

## Ad row to the bottom of the data frame indicating the label
## (do you know a prettier solution for this??)

dataC <- rbind(rna_expr_data_C_part, "cancer")
dataN <- rbind(rna_expr_data_N, "normal")

## Combine the two data frames
dataFull <- cbind(dataC, dataN)


## Remove all rows, that contain at least one zero
## are there any zeros in rows? (except for last row of the data frame)
tmp <- apply(dataFull[-nrow(dataFull),], 1, function(x)!any(x == 0))
## this is a logical vector, indicating if there are zeros in a row

##also, we want to keep the last row keep last row
tmp <- c(tmp, TRUE)

## kick out all the rows with zeros
dataFull1 = dataFull[tmp,]


## save this
save(dataFull1, file = "CountDataWithoutZeros.RData")
