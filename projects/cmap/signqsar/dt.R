library("rpart")
source("signatures.R")

# read the file content, but save it as separate objects not to be modified
mat = signature.read.to.matrix('signatures.txt', size=4629, header=TRUE)
xNamesCSV = read.csv("compoundNames.txt", sep="#")[,1]
compoundsCSV = read.csv("../data/Compounds.csv")

# create initial working objects
x <- mat
xNames <- xNamesCSV

cellline = "HL60"

cat("dim(x): ", dim(x), "\n")
cat("length(xNames): ", length(xNames), "\n")

compounds <- compoundsCSV
compounds <- compounds[-6,]
compounds <- compounds[1:223,]
compounds <- compounds[-222,] # duplicate 
compounds <- compounds[-199,] # duplicate 
y = compounds[,
  c("Response_MCF7", "Response_PC3", "Response_HL60")
]
yNames = compounds[,c("ChemBank_ID.1.")]
cat("dim(y): ", dim(y), "\n")
cat("length(yNames): ", length(yNames), "\n")

# OK, now find the common compounds
commonNames = yNames[yNames %in% xNames]
x = x[xNames %in% commonNames,]
xNames = xNames[xNames %in% commonNames]
y = y[yNames %in% xNames,]
yNames = yNames[yNames %in% commonNames]
cat("dim(x): ", dim(x), "\n")
cat("length(xNames): ", length(xNames), "\n")
cat("dim(y): ", dim(y), "\n")
cat("length(yNames): ", length(yNames), "\n")

# select the cell line
response = paste("Response_", cellline, sep="")
y = y[,response]

# remove descriptor columns with NAs (not applicable to signatures)
#removeNADesc = function(x) {
#  TRUE %in% is.na(x)
#}
#binaryFilter = apply(x, removeNADesc, MARGIN=2)
#x = x[,!binaryFilter]

# remove descriptor columns with zeros (not applicable to signatures)
#removeZero = function(x) {
#  length(unique(x)) == 1
#}
#binaryFilter = apply(x, removeZero, MARGIN=2)
#x = x[,!binaryFilter]

# remove rows with no end point
ynaFilter = is.na(y) == FALSE
x = x[ynaFilter,]
xNames = xNames[ynaFilter]
y = y[ynaFilter]
yNames = yNames[ynaFilter]

x = as.matrix(x)

# make it a classification problem


cat("dim(x): ", dim(x), "\n")
cat("length(xNames): ", length(xNames), "\n")
cat("length(y): ", length(y), "\n")
cat("length(yNames): ", length(yNames), "\n")

# all vars
for (i in 1:5) {
  test.set = sample(nrow(x),ceiling(length(y)/5)) # ~20% test set
  cat(paste("#mols:", length(test.set), "\n"))
  train.x = x[-test.set,]
  train.y = y[-test.set]
  test.x = x[test.set,]
  test.y = y[test.set]
  tree.model = rpart(train.y ~ train.x, minsplit=10)
  test.x = x[test.set,]
  test.y = y[test.set]
  predicted.test.y = predict(tree.model, data.frame(test.x));

  tp = length(which(((predicted.test.y == 1) & (test.y == 1)) == TRUE))
  fp = length(which(((predicted.test.y == 1) & (test.y == 0)) == TRUE))
  tn = length(which(((predicted.test.y == 0) & (test.y == 0)) == TRUE))
  fn = length(which(((predicted.test.y == 0) & (test.y == 1)) == TRUE))

  cat(paste(
     "tp:", tp,
     "tn:", tn,
     "fp:", fp,
     "fn:", fn,
     "\n"
  ))

  a1 = tp / (tp + fn)
  a2 = tn / (tn + fp)

  cat(paste(
     "a1:", a1,
     "a2:", a2,
     "\n"
  ))
}
