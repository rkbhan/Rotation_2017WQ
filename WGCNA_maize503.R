##load data from online 
##knit to html looks cool
setwd("/Users/apple/Desktop/UC_Davis/Lab/2017WQ_Rotation/")
expdata <- read.table("./data/maize_503genotype_cleaned_expression_values_FPKM.txt", sep = "\t", header = TRUE)
library(WGCNA)
options(stringsAsFactors = FALSE)
rm(missing)
dim(expdata)
names(expdata)
expdatawork <- as.data.frame(t(expdata[,-c(1:4)]))
names(expdatawork) <- expdata$gene
View(expdatawork)

##take out genes with <0.5 expression in any lines
expdatawork[is.na(expdatawork)] <- 0
obs.index <- apply(expdatawork,2, function(x) all(x>0.5))
summary(obs.index)
expdatanew <- expdatawork[,obs.index]
dim(expdatanew)

#log transforming the expression data
expdatalog <- log(expdatanew)
range(expdatalog)

save(expdatalog, file="expdata_logtrans.RData")
#load(file = "expdata_logtrans.RData")
#View(expdatalog)

gsg = goodSamplesGenes(expdatalog, verbose = 3)
gsg$allOK
##if gsg#allOK returns false: 
#if (!gsg$allOK)
#{
  # Optionally, print the gene and sample names that were removed:
  #if (sum(!gsg$goodGenes)>0)
    #printFlush(paste("Removing genes:", paste(names(expdata0)[!gsg$goodGenes], collapse = ", ")));
  #if (sum(!gsg$goodSamples)>0)
    #printFlush(paste("Removing samples:", paste(rownames(expdata0)[!gsg$goodSamples], collapse = ", ")));
  ## Remove the offending genes and samples from the data:
  #expdata0 = expdata0[gsg$goodSamples, gsg$goodGenes]
#}
Tree1 <- hclust(dist(expdatalog), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(Tree1, main = "expression data clustering to detect outliers", sub ="", xlab = "",cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
##If there are obvious outliers, remove them with:
#clust <- cutreeStatic(Tree1, cutHeight = XXX, minSize = 10)
#table(clust)
#keeplines <- (clust==1)
#expdatalog <- expdatalogkeeplines,]
nGenes <- ncol(expdatalog)
nLines <- nrow(expdatalog)


##add a new trait: differece between years
## Loading phenotype data
td <- read.csv("~/Desktop/UC_Davis/Lab/2017WQ_Rotation/data/maize_503_phenotypic_data.csv", header = TRUE)
View(td)

## model for year effect (line~year*rep(use 1,2 in year1 and 3,4 in year2), with lsmeans function in package lsmeans) instead of removing na's. 
trait2008 <- td[td$Year=="2008",]
trait2009 <- td[td$Year=="2009",]
trait2009$Rep[trait2009$Rep==1] <- 3
trait2009$Rep[trait2009$Rep==2] <- 4
td <- rbind(trait2008, trait2009)

# change the class 
td[,7:20] = sapply(7:20, function(i) as.numeric(as.character(td[,i])))
View(td)
td[,3] = as.factor(td[,3])

#Remove trait "kernel300"
td = td[,-8]

#removing big gaps (can't think of a smarter way...)
td.clean <- td[-c(2073:2142),]
td.clean <- td.clean[-c(1360:1428),]
td.clean <- td.clean[-c(646:714),]
td.clean <- td.clean[-c(2580:2648),]

#Calculating weighted mean with lsmeans
library(lsmeans)

pred.td = as.data.frame(td.clean[c(1:646),3])
colnames(pred.td) <- "Entry"
for(i in 7:19){
  y <- td.clean[,i]
  trait.lm <- lm(y ~ 0 + Entry + Year*Rep, data = td.clean, na.action = na.omit)
  #grid <- ref.grid(trait.lm)
  #grid.sum = summary(grid)
  pred.test <- lsmeans(trait.lm, "Entry")
  pred.test.t = summary(pred.test)[,c(1:2)]
  colnames(pred.test.t)[2] <- colnames(td.clean)[i]
  pred.td = merge(pred.td,pred.test.t,by = "Entry")
  #pos = is.na(y) 
  #pos = (1:nrow(td))[pos]
  #y[pos] = pred.test.t[match(as.character(td[pos,3]),as.character(pred.test.t[,1])),2]
  #y[pos] = unlist(sapply(pos,function(x) pred.test.t[as.character(pred.test.t[,1])==as.character(td[x,3]),2]))
  #td[,i] = y
}

#removed duplicated lines
pos <- duplicated(pred.td)
pred.td.nodup <- pred.td[!pos,]

save(pred.td.nodup, file="trait_impute.RData")
write.csv(pred.td.nodup, file="trait_impute.csv")

##IMPORTANT: some of the line names are in different formats between the phenotype data and the gene expression data
##I have manually changed the line names in the trait_impute.csv file (probably not the best thing to do...). Therefore here we need to load the "trait_impute.csv"file again.
#####write in the script what got changed

pred.td.nodup <- read.csv("trait_impute.csv", header = TRUE)
View(pred.td.nodup)

#Form a data frame analogous to expression data that will hold the phenotype traits

Entry = rownames(expdatalog)
phenoRows = match(Entry, pred.td.nodup$Entry)
table(is.na(phenoRows))
pheno = pred.td.nodup[phenoRows,-1]
apply(pred.td.nodup,1,function(x) length(is.na(x)[is.na(x)==TRUE]))
pheno = pheno[complete.cases(pheno),]
dim(pheno) #We have 361 lines with complete 14 phenotype data and complete expression data
rownames(pheno) = pheno[,1]
pheno = pheno[,-1]
collectGarbage()

#Here we are going to create a gene expression dataset with just the lines that we have complete phenutyping data for
expRows = match(pred.td.nodup$Entry,Entry)
table(is.na(expRows))
expPheno = expdatalog[expRows,-1]
View(expPheno)
expPheno = expPheno[complete.cases(expPheno),]
dim(expPheno)
collectGarbage()

#Visualizing how the phentypes related to the dendrogram
sampleTree2 = hclust(dist(expPheno), method = "average")
traitColors = numbers2colors(pheno, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(pheno),
                    main = "Preliminaty dendrogram and trait heatmap_361lines")

save(expPheno, pheno, file = "361line_phenotype_and_expression.RData")
#load("361line_phenotype_and_expression.RData")

##At this point, we have:
#490 lines with complete phenotypes ("pred.td.nodup", saved in "trait_impute.RData");
#503 lines with complete expression data ("expdatalog", saved in "expdata_logtrans.RData");
#361 lines (at least right now) with complete both phenotype and expression dta ("expPheno" and "pheno" saved in "361line_phenotype_and_expression.RData")


##Constructing separate models for the 490-line expression dataset and the 361-line dataset. 
##The 361-line dataset will be used for correlation with phenotype data later.

#490_Step-by-step network construction and module detection
power <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(expdatalog, powerVector = power, verbose = 5)
##pdf("home/rongkui/topology")
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=power,cex=cex1,col="red")
#Save plot here
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1,col="red")

#361_Step-by-step network construction and module detection
power <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(expPheno, powerVector = power, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="361 Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("361 Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=power,cex=cex1,col="red")
#Save plot here
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="361 Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("361 Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=power, cex=cex1,col="red")

#490 Co-expression similarity and adjacency
softPower = 3
adjacency_490 = adjacency(expdatalog, power = softPower);
save(adjacency_490, file="490_adjacency.RData")

#361 Co-expression similarity and adjacency
softPower = 3
adjacency_361 = adjacency(expPheno, power = softPower);
save(adjacency_361, file="361_adjacency.RData")

##In case of crashing:
#load("expdata_logtrans.RData")
#load("trait_impute.RData")
#load("adjacency.RData")
#library(WGCNA)

#490_Topological Overlap Matrix (TOM)
TOM_490 = TOMsimilarity(adjacency_490)
dissTOM_490 = 1-TOM_490

#361_Topological Overlap Matrix (TOM)
TOM_361 = TOMsimilarity(adjacency_361)
dissTOM_361 = 1-TOM_361

save(dissTOM_361,dissTOM_490,file = "dissTOMs.RData")
#load("disTOMs.RData")

############490############

##490_Clustering using TOM
geneTree_490 = hclust(as.dist(dissTOM_490), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree_490, xlab="", sub="", main = "490-line Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize_490 = 30
dynamicMods_490 = cutreeDynamic(dendro = geneTree_490, distM = dissTOM_490,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize_490)
table(dynamicMods_490)
dynamicColors_490 = labels2colors(dynamicMods_490)
table(dynamicColors_490)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_490, dynamicColors_490, "490-line Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "490-line Gene dendrogram and module colors")

###########361###########

##361_Clustering using TOM
geneTree_361 = hclust(as.dist(dissTOM_361), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree_361, xlab="", sub="", main = "361-line Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize_361 = 30
dynamicMods_361 = cutreeDynamic(dendro = geneTree_361, distM = dissTOM_361,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize_361)
table(dynamicMods_361)
dynamicColors_361 = labels2colors(dynamicMods_361)
table(dynamicColors_361)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_361, dynamicColors_361, "361-line Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "361-line Gene dendrogram and module colors")

###########490############

##490 Merging modules whose expression profiles are very similar
MEList_490 = moduleEigengenes(expdatalog, colors = dynamicColors_490)
MEs_490 = MEList_490$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss_490 = 1-cor(MEs_490)
# Cluster module eigengenes
METree_490 = hclust(as.dist(MEDiss_490), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree_490, main = "490-linr Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres_490 = 0.4
# Plot the cut line into the dendrogram
abline(h=MEDissThres_490, col = "red")
# Call an automatic merging function
merge_490 = mergeCloseModules(expdatalog, dynamicColors_490, cutHeight = MEDissThres_490, verbose = 3)
# The merged module colors
mergedColors_490 = merge_490$colors
# Eigengenes of the new merged modules:
mergedMEs_490 = merge_490$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree_490, cbind(dynamicColors_490, mergedColors_490),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors_490 = mergedColors_490
# Construct numerical labels corresponding to the colors
colorOrder_490 = c("grey", standardColors(50));
moduleLabels_490 = match(moduleColors_490, colorOrder_490)-1;
MEs_490 = mergedMEs_490;
# Save module colors and labels for use in subsequent parts
save(MEs_490, moduleLabels_490, moduleColors_490, geneTree_490, file = "maize490-networkConstruction-stepByStep.RData")

##########361############

##361 Merging modules whose expression profiles are very similar
MEList_361 = moduleEigengenes(expPheno, colors = dynamicColors_361)
MEs_361 = MEList_361$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss_361 = 1-cor(MEs_361)
# Cluster module eigengenes
METree_361 = hclust(as.dist(MEDiss_361), method = "average")
# Plot the result
sizeGrWindow(7, 6)
plot(METree_361, main = "361-line Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres_361 = 0.4
# Plot the cut line into the dendrogram
abline(h=MEDissThres_361, col = "red")
# Call an automatic merging function
merge_361 = mergeCloseModules(expPheno, dynamicColors_361, cutHeight = MEDissThres_361, verbose = 3)
# The merged module colors
mergedColors_361 = merge_361$colors
# Eigengenes of the new merged modules:
mergedMEs_361 = merge_361$newMEs

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree_361, cbind(dynamicColors_361, mergedColors_361),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors_361 = mergedColors_361
# Construct numerical labels corresponding to the colors
colorOrder_361 = c("grey", standardColors(50));
moduleLabels_361 = match(moduleColors_361, colorOrder_361)-1;
MEs_361 = mergedMEs_361;
# Save module colors and labels for use in subsequent parts
save(MEs_361, moduleLabels_361, moduleColors_361, geneTree_361, file = "maize361-networkConstruction-stepByStep.RData")

#############################

##Relating model to phenotype
##The 361-line datasets are used here
nGenes = ncol(expPheno)
nSamples = nrow(expPheno)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(expPheno, moduleColors_361)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, pheno, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                   signif(moduleTraitPvalue,1),")",sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6,8.5,3,3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(pheno),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

##Gene significance of module membership
##Only analyzing GDD, stalk, leaf-number, stover_yield, ear_internode, and precentage_adult. 
## Define variable GDD containing the GDD column of pheno data
GDD = as.data.frame(pheno$GDD)
names(GDD) = "GDD"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(expPheno, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(expPheno, GDD, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(GDD), sep="");
names(GSPvalue) = paste("p.GS.", names(GDD), sep="")

# identifying genes with high MM and GS in lightyellow module
module_GDD = "lightyellow"
column_GDD = match(module_GDD, modNames)
modulesGenes_GDD = moduleColors_361==module_GDD
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes, column_GDD]),
                   abs(geneTraitSignificance[modulesGenes, 1]),
                   xlab = paste("Module Membership in", module_GDD, "module_GDD"),
                   ylab = "Gene significance GDD",
                   main = paste("GDD Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "khaki1")
##cor = 0.5, p= 0.032
# identifying genes with high MM and GS in black module
module_GDD2 = "black"
column_GDD2 = match(module_GDD2, modNames)
modulesGenes_GDD2 = moduleColors_361==module_GDD2
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_GDD2, column_GDD2]),
                   abs(geneTraitSignificance[modulesGenes_GDD2, 1]),
                   xlab = paste("Module Membership in", module_GDD2, "module_GDD2"),
                   ylab = "Gene significance GDD",
                   main = paste("GDD Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_GDD2)
#nothing significant I guess

# Define variable stalk containing the stalk column of pheno data
stalk = as.data.frame(pheno$stalk)
names(stalk) = "stalk"
# names (colors) of the modules
# modNames = substring(names(MEs), 3)
geneModuleMembership_stalk = as.data.frame(cor(expPheno, MEs, use = "p"))
MMPvalue_stalk = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_stalk), nSamples))
names(geneModuleMembership_stalk) = paste("MM", modNames, sep="");
names(MMPvalue_stalk) = paste("p.MM", modNames, sep="");
geneTraitSignificance_stalk = as.data.frame(cor(expPheno, stalk, use = "p"));
GSPvalue_stalk = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_stalk), nSamples));
names(geneTraitSignificance_stalk) = paste("GS.", names(stalk), sep="");
names(GSPvalue_stalk) = paste("p.GS.", names(stalk), sep="")

# identifying genes with high MM and GS in green module
module_stalk = "green"
column_stalk = match(module_stalk, modNames)
modulesGenes_stalk = moduleColors_361==module_stalk
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_stalk, column_stalk]),
                   abs(geneTraitSignificance[modulesGenes_stalk, 1]),
                   xlab = paste("Module Membership in", module_stalk, "module_stalk"),
                   ylab = "Gene significance stalk",
                   main = paste("Stalk Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_stalk)

# identifying genes with high MM and GS in cyan module
module_stalk2 = "cyan"
column_stalk2 = match(module_stalk2, modNames)
modulesGenes_stalk2 = moduleColors_361==module_stalk2
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_stalk2, column_stalk2]),
                   abs(geneTraitSignificance[modulesGenes_stalk2, 1]),
                   xlab = paste("Module Membership in", module_stalk2, "module_stalk2"),
                   ylab = "Gene significance stalk2",
                   main = paste("Stalk2 Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_stalk2)

# identifying genes with high MM and GS in grey60 module
module_stalk3 = "grey60"
column_stalk3 = match(module_stalk3, modNames)
modulesGenes_stalk3 = moduleColors_361==module_stalk3
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_stalk3, column_stalk3]),
                   abs(geneTraitSignificance[modulesGenes_stalk3, 1]),
                   xlab = paste("Module Membership in", module_stalk, "module_stalk3"),
                   ylab = "Gene significance stalk3",
                   main = paste("Stalk3 Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_stalk3)

# Define variable leafNo containing the leaf_number column of pheno data
leafNo = as.data.frame(pheno$leaf_number)
names(leafNo) = "Leaf_Number"
# names (colors) of the modules
# modNames = substring(names(MEs), 3)
geneModuleMembership_leafNo = as.data.frame(cor(expPheno, MEs, use = "p"))
MMPvalue_leafNo = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_leafNo), nSamples))
names(geneModuleMembership_leafNo) = paste("MM", modNames, sep="");
names(MMPvalue_leafNo) = paste("p.MM", modNames, sep="");
geneTraitSignificance_leafNo = as.data.frame(cor(expPheno, leafNo, use = "p"));
GSPvalue_leafNo = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_leafNo), nSamples));
names(geneTraitSignificance_leafNo) = paste("GS.", names(leafNo), sep="");
names(GSPvalue_leafNo) = paste("p.GS.", names(leafNo), sep="")

# identifying genes with high MM and GS in cyan module
module_leafNo = "cyan"
column_leafNo = match(module_leafNo, modNames)
modulesGenes_leafNo = moduleColors_361==module_leafNo
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_leafNo, column_leafNo]),
                   abs(geneTraitSignificance[modulesGenes_leafNo, 1]),
                   xlab = paste("Module Membership in", module_leafNo, "module_leaf_number"),
                   ylab = "Gene significance leaf number",
                   main = paste("Leaf number Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_leafNo)

# identifying genes with high MM and GS in black module
module_leafNo2 = "black"
column_leafNo2 = match(module_leafNo2, modNames)
modulesGenes_leafNo2 = moduleColors_361==module_leafNo2
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_leafNo2, column_leafNo2]),
                   abs(geneTraitSignificance[modulesGenes_leafNo2, 1]),
                   xlab = paste("Module Membership in", module_leafNo2, "module_leaf_number 2"),
                   ylab = "Gene significance leaf number 2",
                   main = paste("Leaf number Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_leafNo2)

# Define variable stover containing the stover_yield column of pheno data
stover = as.data.frame(pheno$stover_yield)
names(stover) = "Stover_yield"
# names (colors) of the modules
# modNames = substring(names(MEs), 3)
geneModuleMembership_stover = as.data.frame(cor(expPheno, MEs, use = "p"))
MMPvalue-stover = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_stover), nSamples))
names(geneModuleMembership_stover) = paste("MM", modNames, sep="");
names(MMPvalue_stover) = paste("p.MM", modNames, sep="");
geneTraitSignificance_stover = as.data.frame(cor(expPheno, stover, use = "p"));
GSPvalue_stover = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_stover), nSamples));
names(geneTraitSignificance_stover) = paste("GS.", names(stover), sep="");
names(GSPvalue_stover) = paste("p.GS.", names(stover), sep="")

# identifying genes with high MM and GS in grey60 module
module_stover = "grey60"
column_stover = match(module_stover, modNames)
modulesGenes_stover = moduleColors_361==module_stover
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_stover, column_stover]),
                   abs(geneTraitSignificance[modulesGenes_stover, 1]),
                   xlab = paste("Module Membership in", module_stover, "module_stover_yield"),
                   ylab = "Gene significance stover yield",
                   main = paste("Stover Yield Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_stover)

# identifying genes with high MM and GS in blue module
module_stover2 = "blue"
column_stover2 = match(module_stover2, modNames)
modulesGenes_stover2 = moduleColors_361==module_stover2
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_stover2, column_stover2]),
                   abs(geneTraitSignificance[modulesGenes_stover2, 1]),
                   xlab = paste("Module Membership in", module_stover2, "module_stover2"),
                   ylab = "Gene significance stover yield 2",
                   main = paste("Stover yield Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_stover2)

# Define variable internode containing the ear_internode column of pheno data
internode = as.data.frame(pheno$ear_internode)
names(internode) = "Ear-internode"
# names (colors) of the modules
# modNames = substring(names(MEs), 3)
geneModuleMembership_internode = as.data.frame(cor(expPheno, MEs, use = "p"))
MMPvalue_internode = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_internode), nSamples))
names(geneModuleMembership_internode) = paste("MM", modNames, sep="");
names(MMPvalue_internode) = paste("p.MM", modNames, sep="");
geneTraitSignificance_internode = as.data.frame(cor(expPheno, internode, use = "p"));
GSPvalue_internode = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_internode), nSamples));
names(geneTraitSignificance_internode) = paste("GS.", names(internode), sep="");
names(GSPvalue_internode) = paste("p.GS.", names(internode), sep="")

# identifying genes with high MM and GS in lightgreen module
module_internode = "lightgreen"
column_internode = match(module_internode, modNames)
modulesGenes_internode = moduleColors_361==module_internode
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_internode, column_internode]),
                   abs(geneTraitSignificance[modulesGenes_internode, 1]),
                   xlab = paste("Module Membership in", module_internode, "module_ear_internode"),
                   ylab = "Gene significance ear internode",
                   main = paste("Ear internode Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_internode)

# identifying genes with high MM and GS in black module
module_internode2 = "black"
column_internode2 = match(module_internode2, modNames)
modulesGenes_internode2 = moduleColors_361==module_internode2
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_internode2, column_internode2]),
                   abs(geneTraitSignificance[modulesGenes_internode2, 1]),
                   xlab = paste("Module Membership in", module_internode2, "module_ear_internode2"),
                   ylab = "Gene significance ear internode2",
                   main = paste("Ear internode2 Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_internode2)

# Define variable adult containing the percentage_adult column of pheno data
adult = as.data.frame(pheno$percentage_adult)
names(adult) = "Percentage_adult"
# names (colors) of the modules
# modNames = substring(names(MEs), 3)
geneModuleMembership_adult = as.data.frame(cor(expPheno, MEs, use = "p"))
MMPvalue_adult = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_adult), nSamples))
names(geneModuleMembership_adult) = paste("MM", modNames, sep="");
names(MMPvalue_adult) = paste("p.MM", modNames, sep="");
geneTraitSignificance_adult = as.data.frame(cor(expPheno, adult, use = "p"));
GSPvalue_adult = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_adult), nSamples));
names(geneTraitSignificance_adult) = paste("GS.", names(adult), sep="");
names(GSPvalue_adult) = paste("p.GS.", names(adult), sep="")

# identifying genes with high MM and GS in midlightblue module
module_adult = "midnightblue"
column_adult = match(module_adult, modNames)
modulesGenes_adult = moduleColors_361==module_adult
sizeGrWindow(7,7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[modulesGenes_adult, column_adult]),
                   abs(geneTraitSignificance[modulesGenes_adult, 1]),
                   xlab = paste("Module Membership in", module_adult, "module_adult_percentage"),
                   ylab = "Gene significance adult percentage",
                   main = paste("Adult percentage Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module_adult)

##Here we are on WGCNA tutor 3: need gene annotation file. 

###########Step 5: visualization of networks ##########

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
# dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM_361 = dissTOM_361^7
# Set diagonal to NA for a nicer plot
diag(plotTOM_361) = NA;
# Call the plot function
sizeGrWindow(9,9) 
pdf(heatmap.pdf) 
TOMplot(plotTOM_361, geneTree_361, moduleColors_361, main = "Network heatmap plot, 361 lines, all genes")
dev.off()

###For GWAS: cleaning genotyping data
cat(readLines("maize_503genotypes_cleaned_SNP.txt"),sep = "\n")
geno = read.table("maize_503genotypes_cleaned_SNP.txt")

##1/31/17 ToDo for the week
##Finish WGCNA: Need gene annotation file

##Find out Eigengene values of mudules (there should be a value for each module for each individual)
##Try downloading the genotype data from the Hirsch paper (check for name differences)
##Do PCA on the genotype data and put the first couple PC's into TASSEL for GWAS (as a subsitution for kinship matrix) /generate kinship matrix of the genotype data with genotype data
##GWAS