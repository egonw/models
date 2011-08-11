#options(java.parameters = c("-Xmx4000m"))
library(rcdk)
load("../rdf/cbmMols.R")
source("../../cdkpropdist/plot.propDist.R")

alldata = getChemBioModData()

query = "
SELECT ?title ?PC3 WHERE {
  ?compound <http://purl.org/dc/terms/title> ?title ;
    <http://egonw.imm.ki.se/chembiomod/nci60-prop/hasMeasurement> ?m1 .
  ?m1 <http://egonw.imm.ki.se/chembiomod/nci60-prop/hasValue> ?PC3 ;
    <http://egonw.imm.ki.se/chembiomod/nci60-prop/inCellLine> ?line .
  FILTER regex(?line, \"PC3\")
}
"

pc3 = sparql.rdf(alldata, query, rowvarname="title")

mols <- load.molecules(molfiles=c("../data/CMAP_ALL_names.sdf"))

molNames = lapply(mols, get.title)
  
toxic = which(molNames %in% names(pc3[pc3[,"PC3"] < -5,]))
nontoxic = which(molNames %in% names(pc3[pc3[,"PC3"] >= -5,]))
all = 1:length(mols)

plot.propdist(
  mols,
  selections=list(toxic, nontoxic, all),
  plotColors=c("red", "green", "black"),
  descriptor="org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
  main="", xlab="Weight"
)
