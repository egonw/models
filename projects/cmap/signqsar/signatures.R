signature.lf <- function (line, size=128) {
    p <- regexpr("\\[([:0-9,\\s]*)\\]", line, perl = T)
    s <- gsub(",", "", substr(line, p + 1, p + attr(p, "match.length") - 2))
    list <- strsplit(s, " ")[[1]]
    row <- rep(0, size)
    for (listBit in list) {
      bit = strsplit(listBit, ":")[[1]]
      if (as.numeric(bit[1]) > size) {
        cat(paste("Bit doesn't fit matrix:", as.numeric(bit[1]), "\n"));
      }
      row[as.numeric(bit[1])] <- as.numeric(bit[2])
    }
    row
}

signature.read.to.matrix <-
function(f="signatures.txt", size=128, header=FALSE) {
    fcon <- file(description = f, open = "r")
    lines = readLines(fcon, n = -1)
    if (header)
        lines = lines[-1]
    c = 1
    mat = c()
    for (line in lines) {
        row <- signature.lf(line, size)
        mat = rbind(mat, row)
        c <- c + 1
    }
    close(fcon)
    mat
}

