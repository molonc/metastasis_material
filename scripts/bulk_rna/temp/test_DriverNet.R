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


dim(samplePatientMutationMatrix)
head(samplePatientMutationMatrix[1:3,1:3])
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
driversList = computeDrivers(samplePatientMutationMatrix, 
                             samplePatientOutlierMatrix,
                             sampleInfluenceGraph,
                             outputFolder=NULL, printToConsole=FALSE)

drivers(driversList)[1:10]

input_dir <- "/Users/miu/Documents/workspace/projects_BCCRC/hakwoo_project/metastasis_material/materials/bulkRNAseq/test_DriverNet/"
driversList <- readRDS(paste0(input_dir, 'drivernet_driversList.rds'))

length(driversList@actualEvents$PTEN[1:3,2])

driversList@actualEvents$EGFR

###################################################
### code chunk number 5: options
###################################################
options(width=60)


###################################################
### code chunk number 6: preliminaries
###################################################
# random permute the gene labels to compute p-values
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
# dir.create(input_dir)
saveRDS(res, paste0(input_dir, 'drivernet_tuto_res.rds'))
saveRDS(driversList, paste0(input_dir, 'drivernet_driversList.rds'))
saveRDS(randomDriversResult, paste0(input_dir, 'drivernet_randomDriversResult.rds'))



# influence graph: using pathway ref to create gene-gene connection
# binary matrix, 1: connected, 0: not connected. Same gene - connected 
# samplePatientMutationMatrix: using cnv gene mapping matrix, 
# test with only cis genes + genes in pathway first to see if algo works? 
# how many cis genes are in pathway sets?
# binary form: if change in CN between C versus B, true else false, check again binary format 
# samplePatientOutlierMatrix: 2 samples first, if up-reg true, else false

get_influence_graph <- function(pathway_fn=NULL){
  if(is.null(pathway_fn)){
    # pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/c2.cp.kegg.v7.1.symbols.gmt'  
    pathway_fn = '/home/htran/storage/datasets/drug_resistance/rna_results/biodatabase/pathway_set/h.all.v7.0.symbols.gmt'  
  }
  
  ref_set <- fgsea::gmtPathways(pathway_fn)
}
