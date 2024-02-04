library(rrBLUP)

#####reading M
M <- read.csv("M.csv", sep = ',', header = FALSE)
colnames(M) <- 1:ncol(M)
rownames(M) <- 1:nrow(M)
M <- as.matrix(M)

#####reading y
y1 <- read.csv2("Input.csv", sep = ',', header = FALSE)
rownames(y1) <- 1:nrow(y1)
y1 <- y1[,1]
y1 <- as.numeric(as.character(y1))

##### writing answers
ans <- mixed.solve(y1,Z=M)
results <- M %*% as.matrix(ans$u) + ans$beta[1]
write.csv(results,"GEBV.csv", row.names=F)
