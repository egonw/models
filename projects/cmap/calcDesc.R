# options(java.parameters = c("-Xmx4000m"))
library(rcdk)

compounds = read.csv("Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

descNames <- unique(unlist(sapply(
  get.desc.categories(), get.desc.names
)))

iter <- iload.molecules("CMAP_ALL_names.sdf", type = "sdf")
x = c()
xNames = c()
while (hasNext(iter)) {
  mol <- nextElem(iter)
  name <- get.property(mol, "cdk:Title")
  if (name %in% names) {
    allDescs <- eval.desc(mol, descNames)
    x = rbind(x, allDescs)
    xNames = c(xNames, name)
    print(name)
  }
}

yUnclean = compounds[
  compounds[,c("CMAP_chemical_name")] %in% xNames,
  c("CMAP_chemical_name", "Response_MCF7", "Response_PC3", "Response_HL60")
]

