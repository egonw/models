#options(java.parameters = c("-Xmx4000m"))
source("../rdf/cbmMols.R")
source("../../cdkpropdist/plot.propDist.R")

alldata = getChemBioModData("../rdf")
summarize.rdf(alldata)

query = "
SELECT ?title ?PC3 WHERE {
  ?compound <http://purl.org/dc/terms/title> ?title ;
    <http://egonw.imm.ki.se/chembiomod/nci60-prop/hasMeasurement> ?measurement .
  ?measurement
    <http://egonw.imm.ki.se/chembiomod/nci60-prop/inCellLine> \"PC3\" ;
    <http://egonw.imm.ki.se/chembiomod/nci60-prop/hasValue> ?PC3 .
}
"

pc3 = sparql.rdf(alldata, query, rowvarname="title")

library(rcdk)
mols <- load.molecules(molfiles=c("../data/CMAP_ALL_names.sdf"))
print(paste("number of molecules: ", length(mols)))

compounds = read.csv("../data/Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

molNames = lapply(mols, get.title)

smallMols = mols[which(molNames %in% names)]
print(paste("number of molecules:", length(smallMols)))
smallNames = unlist(lapply(smallMols, get.title))

chemInGroups = read.csv("Chemicals_in_Groups_20110919.csv")

classa = which(smallNames %in% chemInGroups[chemInGroups[,"Group"] == "A","Chemical"])
classb = which(smallNames %in% chemInGroups[chemInGroups[,"Group"] == "B","Chemical"])
classc = which(smallNames %in% chemInGroups[chemInGroups[,"Group"] == "C","Chemical"])
classd = which(smallNames %in% chemInGroups[chemInGroups[,"Group"] == "D","Chemical"])
#previously:
#toxic = which(smallNames %in% names(pc3[pc3[,"PC3"] < -5,]))
#nontoxic = which(smallNames %in% names(pc3[pc3[,"PC3"] >= -5,]))
all = 1:length(mols)

plot.propdist(
  mols,
  selections=list(classa, classb, classc, classd, all),
  plotColors=c("red", "green", "yellow", "orange", "black"),
  descriptor="org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
  main="", xlab="Molecular Weight",
  ylab="Relative Distribution"
)


