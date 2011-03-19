library("kohonen")

compounds = read.csv("Compounds.csv")
names = compounds[1:224,c("CMAP_chemical_name")]

x = read.csv("molDesc.csv")
xNames = read.csv("molDescNames.csv")[,2]

yUnclean = compounds[
  compounds[,c("CMAP_chemical_name")] %in% xNames,
  c("CMAP_chemical_name", "Response_MCF7", "Response_PC3", "Response_HL60")
]
yUnclean = yUnclean[-6,]

response = "Response_PC3"
y = yUnclean[,response]

removeNADesc = function(x) {
  TRUE %in% is.na(x)
}
binaryFilter = apply(x, removeNADesc, MARGIN=2)
x = x[,!binaryFilter]

removeZero = function(x) {
  length(unique(x)) == 1
}
binaryFilter = apply(x, removeZero, MARGIN=2)
x = x[,!binaryFilter]

ynaFilter = is.na(y) == FALSE
x = x[ynaFilter,]
y = y[ynaFilter]

x = as.matrix(x)

# all vars
bestR2 = 0.0
for (i in 1:5) {
  test.set = sample(nrow(x),42)
  train.x = x[-test.set,]
  train.y = y[-test.set]
  test.x = x[test.set,]
  test.y = y[test.set]
  xyf.map = xyf(train.x, train.y, grid = somgrid(5, 5, "hexagonal"), rlen=500)
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
  plot(best.trainy, best.predtrainy, main=response)
  points(best.testy, best.predy, col="red")
  abline(0,1)
  dev.print(file=paste("xyfmodel.", i, ".prediction.ps", sep=""), width=6, height=6)
}
save(best.model, file="bestXYFModel.RData")
plot(best.trainy, best.predtrainy, main=response)
points(best.testy, best.predy, col="red")
abline(0,1)
dev.print(file=paste("bestXYFModel.prediction.ps", sep=""), width=6, height=6)

