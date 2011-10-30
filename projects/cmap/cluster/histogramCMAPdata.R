#options(java.parameters = c("-Xmx4000m"))
library(rcdk)

compounds = read.csv("../data/Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

load("../data/GI50_20110415.RData")
yAll = apply(GI50, 1,  function(x) mean(x, na.rm=TRUE))

mols <- load.molecules(c("../data/CMAP_ALL_names.sdf"))
print(paste("number of molecules: ", length(mols)))

load("top20_chems_for_comps_20110915.RData")
for (i in 1:length(top20)) {
  names = top20[[i]]
  
  smallMols = c()
  y = c()
  for (mol in mols) {
    name <- get.property(mol, "cdk:Title")
    if (name %in% names) {
      smallMols = c(smallMols, mol)
      y = c(y, yAll[name])
    }
  }
  print(paste("number of molecules: ", length(smallMols)))
  fps <- lapply(smallMols, get.fingerprint, type = "extended")
  fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
  fp.dist <- 1 - fp.sim

  clustering = hclust(
    as.dist(fp.dist), method="complete"
  )
  hist(
    as.dist(fp.dist), breaks=50, xlab="Tanimoto Distance",
    main="Frequency of Tanimoto Distances",
    sub=paste("Component", names(top20)[i])
  )
} 
  
text(0.10,400, "Few similar\nmolecules", col="red")
text(0.6,1200, "Most are\ndissimilar", col="red")
memb <- cutree(clustering, h=0.5)

identified = identify(clustering)
view.molecule.2d(smallMols[identified[[1]]])

hist(
  fp.dist[identified[[1]],identified[[1]]], breaks=50,
  xlab="Tanimoto Distance",
  main="Frequency of Tanimoto Distances\n(next cluster)", sub="",
  xlim=c(0,1)
)

hist(
  y[identified[[1]]], breaks=30,
  xlab="log GI_50", main="Log GI_50\n(steroid cluster)",
  xlim=c(-10,-2)
)

hist(
  y, breaks=30,
  xlab="log GI_50", main="Log GI_50 (all)"
)
