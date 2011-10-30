source("http://www.bioconductor.org/biocLite.R")
biocLite("ALL")
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

# select some interesting genes

load("/home/egonw/data/data/cmap/xpression.3062.RData")

diffexpr180 = differential.expression.3062[unlist(lapply(smallNames, FUN = function(x) { sample(which(rownames(differential.expression.3062) == x),1) } )),]
# saved as file

pvals = apply(diffexpr180, MARGIN=2, function(x) {t.test(x[steroids],x[-steroids])$p.value} )
interesting.probes = which(pvals<0.05)

heatmap(diffexpr180[steroids,interesting.probes], labCol=rep("",length(pvals)))
