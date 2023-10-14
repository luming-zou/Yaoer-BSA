#!/usr/bin/env Rscript
# Filtering depth and GT
cat("Reading depth file\n")
Sys.time()
dp<-read.table("04.depth.csv", head = F, stringsAsFactors = F)

cat("Reading GT file\n")
Sys.time()
gt<-read.table("gt.GT.FORMAT", head = T, stringsAsFactors = F)

dp$dp_s01<-dp[,5] + dp[,6]
dp$dp_s02<-dp[,7] + dp[,8]
dp$dp_s03<-dp[,9] + dp[,10]
dp$dp_s04<-dp[,11] + dp[,12]
dp<-dp[dp$dp_s01 >= 3 & dp$dp_s01 <= 200 &
	dp$dp_s02 >= 3 & dp$dp_s02 <= 200 &
	dp$dp_s03 >= 3 & dp$dp_s03 <= 50 &
	dp$dp_s04 >= 3 & dp$dp_s04 <= 50&
	gt$S03 != gt$S04&
	gt$S01 !=  './.' &
	gt$S02 !=  './.' &
	gt$S03 !=  './.' &
	gt$S04 !=  './.', 
	]
dp<-dp[,1:12]


## S01
cat("S01 index\n")
Sys.time()
dp_mut<-ifelse(dp[,6]>dp[,5], dp[,6], dp[,5])
dp$index_s01<-dp_mut/(dp[,5] + dp[,6])

## S02
cat("S02 index\n")
Sys.time()
dp_mut<-ifelse(dp[,6]>dp[,5], dp[,8], dp[,7])
dp$index_s02<-dp_mut/(dp[,7] + dp[,8])

dp<-dp[dp$index_s01 > 0.9 & dp$index_s02 < 0.9, ]
dp<-dp[,1:12]
write.table(dp, "filtered_depth.txt", row.names = F, col.names = F, quote = F, sep = "\t")
