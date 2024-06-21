

# Therefore, we use 8 different but well-established datasets: METABRIC56, GSE989357, GSE739058, GSE9605859, GSE1112160, GSE492261, NKI62,63,

# Package installation and library imports
#install.packages("BiocManager")
BiocManager::install("Biobase")
BiocManager::install("preprocessCore")
BiocManager::install("affy")
BiocManager::install("affyPLM")
BiocManager::install("leukemiasEset")
BiocManager::install("gcrma")
BiocManager::install("biomaRt")
BiocManager::install("broom")


library(preprocessCore)
library(Biobase)
library(gcrma)
library(affy)
library(affyPLM) # for normalization


# Set the path to the GSM_data folder
gsmDataFolder <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/GSM_data"

# Get the list of files in the GSM_data folder
fileList <- list.files(gsmDataFolder, full.names = TRUE)


# Read in the data so that the samples are the columns and the genes are the rows
# We will assume that F635.Mean represents the mean fluorescence intensity of the red channel for each spot which quantifies gene expression

# Set the path to the GSM_data folder
gsmDataFolder <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/GSM_data"

# Get the list of .gpr files in the GSM_data folder
fileList <- list.files(gsmDataFolder, pattern = "\\.gpr$", full.names = TRUE)

# Initialize an empty list to store the data from each file
dataList <- list()

# Loop over the files
for (i in seq_along(fileList)) {
  # Read the file
  gprDataRaw <- read.table(fileList[i], sep = "\t", header = TRUE, skip = 31)
  
  # Extract the F635.Mean column and store it in the list
  dataList[[i]] <- gprDataRaw$F635.Mean
}

# Combine the data from all files into a single matrix
dataMatrix <- do.call(cbind, dataList)


# Create a second data matrix to preserve it
dataMatrix2 <- dataMatrix
# Set the column names of the matrix to the filenames (without the path and extension)
colnames(dataMatrix2) <- basename(fileList)
colnames(dataMatrix2) <- sub("\\.gpr$", "", colnames(dataMatrix2))


# now read in the phenotype data (GSE9893_clinical_data.csv)

# Set the path to the .txt file
txtFilePath <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/GSM_data/GSE9893_clinicalData.txt"

# Read the .txt file with fill = TRUE
phenotypeData <- read.table(txtFilePath, sep = "\t", header = TRUE, fill = TRUE)

# Transpose phenotypeData so that columns are samples and rows are variables
transposedPhenoData <- as.data.frame(t(phenotypeData))

# Set the column names to be the 'Tumor.sample' column from the original data
colnames(transposedPhenoData) <- phenotypeData$"Tumor.sample"

# Remove the 'Tumor.sample' row from the transposed data
transposedPhenoData <- transposedPhenoData[-which(rownames(transposedPhenoData) == "Tumor.sample"), ]

# View the first few rows of the data
head(transposedPhenoData)


phenotypeData <- transposedPhenoData


# Now lets change the column name from dataMatrix from GSM number to the EB number

identifierPath <- "C:/Users/patri/OneDrive/Desktop/snpdr_update/snpdr_update/data/GSM_data/identifier.csv"
identifierData <- read.csv(identifierPath, header = FALSE, stringsAsFactors = FALSE)

# Set the column names
colnames(identifierData) <- c("GSM", "EB")

# Replace the GSM numbers with the corresponding EB numbers in the data matrix
colnames(dataMatrix2) <- identifierData$EB[match(colnames(dataMatrix2), identifierData$GSM)]

# Now order the columns of data matrix to match the order of the rows in the phenotype data
dataMatrix2 <- dataMatrix2[, match(colnames(phenotypeData), colnames(dataMatrix2))]

# Transpose phenoData back so that the rows of phenoData equal the columns of dataMatrix2
phenotypeData <- as.data.frame(t(phenotypeData))

# Now lets convert this into an expression set
assayDataMatrix <- dataMatrix2
phenoData <- new("AnnotatedDataFrame", data = phenotypeData)


# Create ExpressionSet
eset <- ExpressionSet(assayData = assayDataMatrix,
                      phenoData = phenoData)



# Save the ExpressionSet object to a file
save(eset, file = "eset.RData")

# Load the ExpressionSet object from a file
load("eset.RData")


# objects in the data
str(eset)
slotNames(eset) # phenoData, etc

dim(eset) # 22680 genes, 155 samples
sampleNames(eset) # cell file names

allPheno <- pData(eset)
head(allPheno)  # look at phenotype info
censorStatus <- allPheno$"State.of.health"  # censor status as state of health (alive or deceased)
summary(censorStatus)
# ALL AML CLL CML NoL 
# 12  12  12  12  12 
featureNames(eset)[1:5] # first 5 gene ensemble ids

ExprData <- exprs(eset) # exprs is an affy function to extract expression data from eset
colnames(ExprData) <- censorStatus  # add phenotype names to matrix


# quantiles function needs eset to operate on
ExprData_quantile <- normalize.quantiles(ExprData)
boxplot(ExprData_quantile, range=0, ylab="raw intensity", main="Quantile Normalized")
#ExprData_quantileLog2 <- log2(exprs(ExprData_quantile))
ExprData_quantileLog2 <- log2(ExprData_quantile)
colnames(ExprData_quantileLog2) <- censorStatus  # add phenotype names to matrix
boxplot(ExprData_quantileLog2, range=0, ylab="log2 intensity", main="Quantile Normalized Log2")
