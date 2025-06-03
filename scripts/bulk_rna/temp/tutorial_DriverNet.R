### R code from vignette source 'DriverNet-Overview.Rnw'
# BiocManager::install("DriverNet")
###################################################
### code chunk number 1: options
###################################################
options(width=60)


###################################################
### code chunk number 2: preliminaries
###################################################
library(DriverNet)

data(samplePatientMutationMatrix)
data(samplePatientOutlierMatrix)
data(sampleInfluenceGraph)
data(sampleGeneNames)
length(sampleGeneNames)
dim(sampleGeneNames)
sampleGeneNames[1:5]
colnames(sampleInfluenceGraph)[1:5]
dim(sampleInfluenceGraph)
sampleInfluenceGraph[1:10,1:10]

sum(sampleGeneNames %in% colnames(samplePatientMutationMatrix))
pids <- rownames(samplePatientMutationMatrix)
pids[1:10]
ext_pids <- sample(pids, 10)
samplePatientMutationMatrix <- samplePatientMutationMatrix[ext_pids,]
samplePatientOutlierMatrix <- samplePatientOutlierMatrix[ext_pids,]



dim(samplePatientMutationMatrix)
head(samplePatientMutationMatrix[1:3,1:3]) # patients x genes, values 0, 1
# with row names and col names are pt and genes

## checking how to define outliers
head(samplePatientOutlierMatrix[1:3,1:3])  # patients x genes, values false, true

dim(sampleInfluenceGraph) # a square mtx with 1 is connected, and 0 is not. Same vertex is 1
sampleInfluenceGraph[1:3,1:3]


cna_sids <- rownames(samplePatientMutationMatrix)
cna_gids <- colnames(samplePatientMutationMatrix)

dim(samplePatientOutlierMatrix)
head(samplePatientOutlierMatrix[1:3,1:3])
sum(samplePatientOutlierMatrix[2,]==T)
exp_sids <- rownames(samplePatientOutlierMatrix)
exp_gids <- colnames(samplePatientOutlierMatrix)

dim(sampleInfluenceGraph)
length(sampleGeneNames)
head(sampleInfluenceGraph[1:3,1:3])

length(exp_gids)
length(cna_gids)
sum(cna_sids %in% exp_sids)
###################################################
### code chunk number 3: options
###################################################
options(width=60)


###################################################
### code chunk number 4: preliminaries
###################################################
# The main function to compute drivers

driversList = DriverNet::computeDrivers(samplePatientMutationMatrix, 
                                        samplePatientOutlierMatrix,
                                        sampleInfluenceGraph,
                                        outputFolder=NULL, printToConsole=FALSE)

drivers(driversList)[1:10]
class(driversList)

input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/test_DriverNet/"
driversList1 <- readRDS(paste0(input_dir, 'drivernet_driversList.rds'))

length(driversList@actualEvents$PTEN)

head(driversList@actualEvents$EGFR)

###################################################
### code chunk number 5: options
###################################################
options(width=60)


###################################################
### code chunk number 6: preliminaries
###################################################
# random permute the gene labels to compute p-values
dim(samplePatientOutlierMatrix)
dim(samplePatientMutationMatrix)
dim(sampleInfluenceGraph)   

length(sampleGeneNames)
randomDriversResult = computeRandomizedResult(
  patMutMatrix=samplePatientMutationMatrix,
  patOutMatrix=samplePatientOutlierMatrix, 
  influenceGraph=sampleInfluenceGraph,
  geneNameList= sampleGeneNames, outputFolder=NULL, 
  printToConsole=FALSE,numberOfRandomTests=20, weight=FALSE, 
  purturbGraph=FALSE, purturbData=TRUE)

randomDriversResult <- readRDS(paste0(input_dir, 'drivernet_randomDriversResult.rds'))

###################################################
### code chunk number 7: options
###################################################
options(width=60)


###################################################
### code chunk number 8: preliminaries
###################################################
# Summarize the results
res = resultSummary(driversList, randomDriversResult, 
                    samplePatientMutationMatrix,sampleInfluenceGraph,
                    outputFolder=NULL, printToConsole=FALSE)


res <- readRDS(paste0(input_dir, 'drivernet_tuto_res.rds'))
res[1:2,1:10]

input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/test_DriverNet/"
input_dir <- "/Users/hoatran/Documents/BCCRC_projects/hakwoo_project/testing/DriverNet/"
# dir.create(input_dir)
saveRDS(res, paste0(input_dir, 'drivernet_tuto_res.rds'))
saveRDS(driversList, paste0(input_dir, 'drivernet_driversList.rds'))
saveRDS(randomDriversResult, paste0(input_dir, 'drivernet_randomDriversResult.rds'))

dim(samplePatientMutationMatrix)
