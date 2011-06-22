library("pls")

# load descriptor data
x = read.csv("../data/molDesc.csv")
xNames = read.csv("../data/molDescNames.csv")[,2]

# load GI_50 values
load("../data/GI50_20110415.RData")
yUnclean = GI50[rownames(GI50) %in% xNames,]

# select cell line MCF7, PC3, HL60
#response = "HL60"
#y = yUnclean[,response]

# take the mean GI_50 over all cell lines
response = "mean log GI50"
y = apply(yUnclean, 1,  function(x) mean(x, na.rm=TRUE))

# remove molecules which have a NA x value
xnaFilter = !apply(x, 1, function(x) any(is.na(x)))
x = x[xnaFilter,]
y = y[xnaFilter]

print(dim(x))
print(length(y))

# remove descriptors which have NA values
removeNADesc = function(x) {
  TRUE %in% is.na(x)
}
binaryFilter = apply(x, removeNADesc, MARGIN=2)
x = x[,!binaryFilter]

# remove descriptors which have only one value
removeZero = function(x) {
  length(unique(x)) == 1
}
binaryFilter = apply(x, removeZero, MARGIN=2)
x = x[,!binaryFilter]

# remove molecules which have a NA y value
ynaFilter = is.na(y) == FALSE
x = x[ynaFilter,]
y = y[ynaFilter]

# remove all molecules with a GI_50 of -4.0
ynaFilter = which(y == -4.0)
x = x[-ynaFilter,]
y = y[-ynaFilter]

x = as.matrix(x)

print(dim(x))
print(length(y))

y[which(y > -5)] = 0  # non-toxic
y[which(y <= -5)] = 1 # toxic

# all vars
for (i in 1:20) {
  test.set = sample(nrow(x),ceiling(length(y)/5)) # ~20% test set
  cat(paste("#mols:", length(test.set), "\n"))
  train.x = x[-test.set,]
  train.y = y[-test.set]
  test.x = x[test.set,]
  test.y = y[test.set]
  tree.model = rpart(train.y ~ train.x)
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
