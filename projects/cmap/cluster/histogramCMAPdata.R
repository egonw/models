#options(java.parameters = c("-Xmx4000m"))
library(rcdk)

compounds = read.csv("../data/Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

mols <- load.molecules(c("../data/CMAP_ALL_names.sdf"))
print(paste("number of molecules: ", length(mols)))

smallMols = c()

for (mol in mols) {
  name <- get.property(mol, "cdk:Title")
  if (name %in% names) {
    smallMols = c(smallMols, mol)
  }
}
print(paste("number of molecules: ", length(smallMols)))

clustering = hclust(as.dist(fp.dist), method="complete")
plot(clustering, hang=-1)
identified = identify(clustering)
view.molecule.2d(smallMols[identified[[1]]])
