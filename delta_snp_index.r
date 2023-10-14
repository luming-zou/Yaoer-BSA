#!/usr/bin/env Rscript
options(scipen = 20)
dp<-read.table("filtered_depth.txt", head = F, sep = "\t", stringsAsFactors = F)

dp_mut<-ifelse(dp[,6]>dp[,5], dp[,6], dp[,5])
index_s01<-dp_mut/(dp[,5] + dp[,6])

dp_mut<-ifelse(dp[,6]>dp[,5], dp[,8], dp[,7])
index_s02<-dp_mut/(dp[,7] + dp[,8])

d_index<-index_s01 - index_s02

p_value<-NULL
depth<-as.matrix(dp[, 5:8])
for (i in 1:dim(dp)[1]){
	x<-matrix(depth[i,], ncol = 2, byrow = T)
	p<-fisher.test(x)$p.value
	p_value<-c(p_value, p)
}

index<-data.frame(CHROM = dp[,1], POS = dp[,2], INDEX_r = index_s01, INDEX_S = index_s02, 
	D_INDEX = d_index)

delta_index<-function(chr, start, end){
	index_temp1<-index_temp[index_temp[,2] > start & index_temp[,2] <= end,]
	if(dim(index_temp1)[1] >= 0){
		m<-sum(index_temp1$D_INDEX)/(end - start)
		return(m)
	}else{
		return(0)
	}
}

fai<-read.table("REFERENCE.fa.fai", head = F, stringsAsFactors = F)

index_all<-NULL
for (i in 1:dim(fai)[1]){
	D_INDEX_CHR<-NULL
	chr<-fai[i, 1]
	print(chr)
	index_temp<-index[index$CHROM == chr,]

	window<-1e3
	step<-7e2

	## windows
	length<-fai[i,2]
	if (length <= window){
		windows<-matrix(c(1, length), ncol = 2)
	}else{
		starts<-seq(1, length, step)
		ends<-starts + window
		windows<-cbind(starts, ends)
		windows<-windows[ends<=length, ]
		windows<-matrix(windows, ncol = 2)
		windows<-rbind(windows, c(max(windows[,1]) + step, length))
	}

	for (j in 1:dim(windows)[1]){
		D_INDEX_CHR<-rbind(D_INDEX_CHR, delta_index(chr, windows[j,1], windows[j,2]))
	}

	index_chr<-data.frame(CHROM = chr, STARTS = windows[,1], ENDS = windows[,2], D_INDEX = D_INDEX_CHR[,1])
	index_all<-rbind(index_all, index_chr)
}
write.table(index_all, "delta_snp_index.csv", row.names = F, sep = "\t", quote = F)
