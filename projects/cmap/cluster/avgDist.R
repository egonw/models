#options(java.parameters = c("-Xmx4000m"))
library(rcdk)

compounds = read.csv("../data/Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

mols <- load.molecules(c("../data/CMAP_ALL_names.sdf"))
print(paste("number of molecules: ", length(mols)))

smallMols = c()
smallNames = c()

for (mol in mols) {
  name <- get.property(mol, "cdk:Title")
  if (name %in% names) {
    smallMols = c(smallMols, mol)
    smallNames = c(smallNames, name)
  }
}
print(paste("number of molecules: ", length(smallMols)))

fps <- lapply(smallMols, get.fingerprint, type='standard')
sim.matrix = fp.sim.matrix(fps)
fp.dist = 1 - sim.matrix

clustering = hclust(as.dist(fp.dist), method="complete")
plot(clustering, hang=-1)
identified = identify(clustering)
view.molecule.2d(smallMols[identified[[1]]])

png(filename="histSimilarities.png")
hist(sim.matrix, main="Compound Similarities", xlab="Tanimoto similarities")
dev.off()

# new code to calculate the tanimoto between the two interesting structures
interestingMols = c()
interestingNames = c()
selectedMols = c("tanespimycin","alvespimycin","geldanamycin")

for (mol in mols) {
  name <- get.property(mol, "cdk:Title")
  if (name %in% selectedMols) {
    interestingMols = c(interestingMols, mol)
    interestingNames = c(interestingNames, name)
    view.molecule.2d(mol)
  }
}
print(paste("number of interesting molecules: ", length(interestingMols)))

fps <- lapply(interestingMols, get.fingerprint, type='standard')
sim.matrix = fp.sim.matrix(fps)
rownames(sim.matrix) = interestingNames
colnames(sim.matrix) = interestingNames
sim.matrix
fp.dist = 1 - sim.matrix
