library("kohonen")

replaceNaNs <- function(x) {
  x[which(is.na(x))] = mean(x, na.rm=T)
  x
}

# load descriptor data
x = read.csv("../data/molDesc.csv")
xNames = read.csv("../data/molDescNames.csv")[,2]

# load GI_50 values
load("../data/GI50_20110415.RData")
yUnclean = GI50[rownames(GI50) %in% xNames,]

# select cell line
#response = "MCF7"
#y = yUnclean[,response]

# use all data
# y = replaceNaNs(yUnclean)

# take the mean GI_50 over all cell lines
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

print(dim(x))
print(length(y))

# remove descriptors which have only one value
removeZero = function(x) {
  length(unique(x)) == 1
}
binaryFilter = apply(x, removeZero, MARGIN=2)
x = x[,!binaryFilter]

# remove molecules which have a NA y value
#ynaFilter = is.na(y) == FALSE
#x = x[ynaFilter,]
#y = y[ynaFilter,]

x = as.matrix(x)

print(dim(x))
print(length(y))

# all vars
bestR2 = 0.0
for (i in 1:5) {
  test.set = sample(nrow(x), ceiling(length(y)/10)) # ~25% test set
  train.x = x[-test.set,]
  train.y = y[-test.set]
  test.x = x[test.set,]
  test.y = y[test.set]
  xyf.map = xyf(train.x, train.y, grid = somgrid(7, 7, "hexagonal"), rlen=500)
  xyf.prediction.train = predict(xyf.map, newdata=train.x)
  predicted.train.y = xyf.prediction.train$prediction
  xyf.prediction = predict(xyf.map, newdata=test.x)
  predicted.test.y = xyf.prediction$prediction
  RMSEP <- sqrt(sum((predicted.test.y-test.y)^2)/length(test.y)); RMSEP
  R2 = cor(train.y, predicted.train.y)^2
  if (R2 > bestR2) {
      bestR2 = R2
      best.model = xyf.map
      best.trainy = train.y
      best.predtrainy = predicted.train.y
      best.testy = test.y
      best.predy = predicted.test.y
  }

  cat(paste("   rank=", round(qr(train.x)$rank, digit=0),
            "     R2=", round(R2, digit=3),
            "  RMSEP=", round(RMSEP, digit=3), "\n", sep=""));
  save.image(file=paste("model.", i, ".RData", sep=""))
  plot(best.trainy, best.predtrainy)
  points(best.testy, best.predy, col="red")
  abline(0,1)
  dev.print(file=paste("xyfmodel.", i, ".prediction.ps", sep=""), width=6, height=6)
}
save(best.model, file="bestXYFModel.RData")
plot(best.trainy, best.predtrainy)
points(best.testy, best.predy, col="red")
abline(0,1)
dev.print(file=paste("bestXYFModel.prediction.ps", sep=""), width=6, height=6)

