library("pls")

x.raw = read.csv("../data/chemicals2components.csv")
cell.select = x.raw[,"Cell.line"] == "PC3" 

duplicate = 166
x = x.raw[cell.select, 4:103]
x = x[-duplicate,]
xNames = x.raw[cell.select, "Chemical"]
xNames = xNames[-duplicate]

compounds = read.csv("../data/Compounds.csv")
names.raw = compounds[1:224,c("CMAP_chemical_name")]
names = names.raw[names.raw %in% xNames]

yUnclean = compounds[
  compounds[,c("CMAP_chemical_name")] %in% xNames,
  c("CMAP_chemical_name", "Response_MCF7", "Response_PC3", "Response_HL60")
]

response = "Response_HL60"
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
bestQ2 = 0.0
for (i in 1:10) {
  test.set = sample(nrow(x),21)
  train.x = x[-test.set,]
  train.y = y[-test.set]
  test.x = x[test.set,]
  test.y = y[test.set]
  maxLV = min(nrow(test.x), ceiling(length(train.y)/5))
  pls.model = plsr(train.y ~ train.x, ncomp=maxLV, validation="CV")
  lv = pls.model$ncomp
  # lv = 15
  test.x = x[test.set,]
  test.y = y[test.set]
  predicted.test.y = predict(pls.model, test.x, ncomp=lv);
  RMSEP <- sqrt(sum((predicted.test.y-test.y)^2)/length(test.y)); RMSEP
  R2 = cor(train.y, pls.model$fit[,,lv])^2
  Q2 = cor(train.y, pls.model$validation$pred[,,lv])^2
  if (Q2 > bestQ2) {
      bestQ2 = Q2
      bestPLS.model = pls.model
      bestPLS.testy = test.y
      bestPLS.predy = predicted.test.y
  }

  cat(paste("   rank=", round(qr(train.x)$rank, digit=0),
            "     LV=", round(lv, digit=3),
            "     R2=", round(R2, digit=3),
            "     Q2=", round(Q2, digit=3),
            "  RMSEP=", round(RMSEP, digit=3), "\n", sep=""));
  save.image(file=paste("model.", i, ".RData", sep=""))
  plot(pls.model, plottype="prediction", ncomp=lv, main=response)
  points(bestPLS.testy, bestPLS.predy, col="red")
  abline(0,1)
  dev.print(file=paste("plsmodel.", i, ".prediction.ps", sep=""), width=6, height=6)
  plot(pls.model, plottype="coeff", ncomp=lv)
  dev.print(file=paste("plsmodel.", i, ".coeff.ps", sep=""), width=6, height=6)
}
save(bestPLS.model, file="bestPLSModel.RData")
plot(bestPLS.model, plottype="prediction", ncomp=lv)
points(bestPLS.testy, bestPLS.predy, col="red")
abline(0,1)
dev.print(file=paste("bestPLSModel.prediction.ps", sep=""), width=6, height=6)
plot(bestPLS.model, plottype="coeff", type="l")
dev.print(file="bestPLSModel.regression.ps", width=6, height=6)

