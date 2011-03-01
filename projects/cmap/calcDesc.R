# options(java.parameters = c("-Xmx4000m"))
library(rcdk)

descNames <- unique(unlist(sapply(
  get.desc.categories(), get.desc.names
)))

iter <- iload.molecules("CMAP_ALL.sdf", type = "sdf")
while (hasNext(iter)) {
  mol <- nextElem(iter)
  allDescs <- eval.desc(mol, descNames)
  x = rbind(x, allDescs)
  print(get.property(mol, "cdk:Title"))
}
