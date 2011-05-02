#options(java.parameters = c("-Xmx4000m"))
library(rcdk)

compounds = read.csv("../data/Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

descNames <- unique(unlist(sapply(
  get.desc.categories(), get.desc.names
)))

fixNitrogenCharge <- function(mol) {
  atoms = get.atoms(mol)
  for (atom in atoms) {
    if (get.symbol(atom) == "N") {
      neighborCount = .jcall(mol, "I", "getConnectedBondsCount", atom) 
      if (neighborCount == 4) {
        # set formal charge to +1
        .jcall(atom, "V", "setFormalCharge", .jnew("java/lang/Integer", "1"))
      }
    }
  }
}

iter <- iload.molecules("../data/CMAP_ALL_names.sdf", type = "sdf")
x = c()
xNames = c()
blacklist = c("Wgamma1.unity", "Wgamma2.unity", "Wgamma3.unity", "WG.unity")
bads = c()
while (hasNext(iter)) {
  mol <- nextElem(iter)
  fixNitrogenCharge(mol)
  name <- get.property(mol, "cdk:Title")
  if (name %in% names) {
    print(name)
    allDescs <- eval.desc(mol, descNames)
    allDescNames = colnames(allDescs)
    allDescNames = allDescNames[!(allDescNames %in% blacklist)]
    allDescs = allDescs[,allDescNames]
    nanDescs = is.na(allDescs)
    # print(rbind(rbind(allDescs, nanDescs), allDescNames))
    if (!all(!nanDescs)) {
      print(paste("", "descriptor with NaN:"))
      print(as.vector(allDescNames[nanDescs]))
      bads = c(name, bads)
    }
    x = rbind(x, allDescs)
    xNames = c(xNames, name)
  }
}

print(paste("number of bad molecules: ", length(bads)))
print(bads)

write.csv(file="../data/molDesc.csv", x)
write.csv(file="../data/molDescNames.csv", xNames)

