#install.packages('pcadapt')
library(pcadapt)

#read data
readalcon <- read.pcadapt('/home/jonas/Documents/Masterthesis/data/Analyses/PCadapt/Alcon.lfmm', type = "lfmm")
readgentian <- read.pcadapt('/home/jonas/Documents/Masterthesis/data/Analyses/PCadapt/Gentian.lfmm', type = "lfmm")
#new lfmm for Alcon with 127 SNPs ('Alcon.lfmm' has only 122 SNPs)
A <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/alconSNP', h = T)
popA <- A$PopID


#create lfmm
A2 <- A[,3:ncol(A)]
refSNPs <- vector()
for (i in 1:254){
  refSNPs <- append(refSNPs, tail(names(sort(table(A2[,i]))), 1))
}

AlcLFMM <- data.frame(matrix(ncol=127, nrow = 661))
for (i in 1:127){
  for (j in 1:661){
    if (A2[j,2*i-1] == 0 && A2[j, 2*i] == 0){
      AlcLFMM[j,i] <- 9
      #if both values are 0, value becomes 9
    }
    if (A2[j,2*i-1] != 0 && A2[j,2*i-1] == refSNPs[2*i-1] && A2[j, 2*i]== refSNPs[2*i]){
      AlcLFMM[j,i] <- 0
      #if both values are same as reference, value becomes 0
    }
    if (A2[j,2*i-1] != 0 && A2[j,2*i-1] != refSNPs[2*i-1] && A2[j, 2*i] != refSNPs[2*i]){
      AlcLFMM[j,i] <- 2
      #if both values are different as reference, value becomes 2
    }
    if (A2[j,2*i-1] != 0 && A2[j,2*i-1] != refSNPs[2*i-1] && A2[j, 2*i] == refSNPs[2*i]){
      AlcLFMM[j,i] <- 1
      #if first value is different as reference, value becomes 1
    }
    if (A2[j,2*i-1] != 0 && A2[j,2*i-1] == refSNPs[2*i-1] && A2[j, 2*i] != refSNPs[2*i]){
      AlcLFMM[j,i] <- 1
      #if second value is different as reference, value becomes 1
    }
  }
}

colnames(AlcLFMM) <- coln
write.table(AlcLFMM, 'AlcLFMM.lfmm', row.names = F, col.names = F, quote = F)

readalcon2 <- read.pcadapt('/home/jonas/Documents/Masterthesis/data/Analyses/PCadapt/AlcLFMM.lfmm', type = "lfmm")

#perform pcadapt
alcon <- pcadapt(input = readalcon, K = 20)
alcon2 <- pcadapt(input = readalcon2, K = 20)
gentian <- pcadapt(input = readgentian, K = 20)


#screeplots
plot(x = alcon, option = 'screeplot', K = 10) #4 PCA's retained
plot(x = alcon2, option = 'screeplot', K = 10) #same
plot(x = gentian, option = 'screeplot', K = 10) #not so clear, 6 or 7?

#adapt K number (see screeplots)
alcon <- pcadapt(input = readalcon, K = 4)
alcon2 <- pcadapt(input = readalcon2, K = 4, min.maf = 0.01)
gentian <- pcadapt(input = readgentian, K = 6, min.maf = 0.01)

?pcadapt
length(alcon2$pvalues[is.na(alcon2$pvalues)]) #25 NAs (0.05), 7 (0.01), 
length(gentian$pvalues[is.na(gentian$pvalues)]) #3 NAs (0.05), 1 (O.O1)

#PCA's
plot(x = alcon, option = 'scores') #3 groups
plot(x = gentian, option = 'scores') #less clear
plot(x = alcon2, option = 'scores')

alcon$gif
gentian$gif
alcon2$gif

#graphical checks
#manhatten
plot(alcon, option ='manhattan') #1 highly significat, 1 suggestive
plot(alcon2, option = 'manhattan')
plot(gentian, option='manhattan') #few very very significant SNPs
#QQplot
plot(alcon, option ='qqplot') ##qqs looks ok, no systematic population structure
plot(alcon2, option = 'qqplot') #quite OK
plot(gentian, option='qqplot')

#histogram
hist(alcon$pvalues, xlab = "p-values", main = 'Distribution of P-values - Alcon', breaks = 50, col = "orange")
#doesnt look to good, low p-values might be false positives
hist(alcon2$pvalues, xlab = "p-values", main = 'Distribution of P-values - Alcon', breaks = 50, col = "orange") 
#quite good
hist(gentian$pvalues, xlab = "p-values", main = 'Distribution of P-values - Gentian', breaks = 50, col = "orange") #high frequency of low p-value SNPs:: outliers!
plot(alcon, option = "stat.distribution")
plot(alcon2, option = "stat.distribution")
plot(gentian, option = "stat.distribution")

#outliers
#bonferonni corrected
padja <- p.adjust(alcon$pvalues, method='bonferroni')
padja2 <- p.adjust(alcon2$pvalues, method = 'bonferroni')
padjg <- p.adjust(gentian$pvalues, method='bonferroni')
alpha <- 0.05 #false discovery rate of 5%
outlieralcon <- which(padja < alpha)
outlieralcon2 <- which(padja2 < alpha)
outliergentian <- which(padjg < alpha)
length(outlieralcon2)
length(outliergentian)

snpsG <- read.table(file = "/home/jonas/Documents/Masterthesis/data/snps.gentian", header = F, sep =" ")
snpsG[, 77]

snpsA <- read.table(file = "/home/jonas/Documents/Masterthesis/data/snps.alcon", header = F, sep =" ")
snpsA[, 108]

#plot mahanolobis distance
#alcon
?pcadapt
library(ggplot2)
library(ggrepel)
dev.off()
a.m <- alcon2$stat
?ifelse
a.p <- ifelse(alcon2$pvalues)
alcon2$pvalues
ggplot()+
  geom_point(data = NULL, aes(x= 1:length(a.m), y = a.m))+
  geom_point(data = NULL, aes(x= outlieralcon2, y = a.m[outlieralcon2]), 
             color ="#116E8A")+
  theme_classic()+
  geom_text_repel(mapping = aes(x = outlieralcon2, y = a.m[outlieralcon2]), 
                  label = outlieralcon2, nudge_y = 5, nudge_x = 2)+
  ylab('Mahalanobis Distance') + xlab('SNP Number')+ggtitle('Alcon')
#warning: 25 NA because of 5% minor allele freq -> 6 at 0.01MAF

g.m <- gentian$stat
ggplot()+
  geom_point(data = NULL, aes(x= 1:length(g.m), y = g.m))+
  geom_point(data = NULL, aes(x= outliergentian, y = g.m[outliergentian]), 
             color ="#116E8A")+
  theme_classic()+
  geom_text_repel(mapping = aes(x = outliergentian, y = g.m[outliergentian]), 
                  label = outliergentian, nudge_y = 5, nudge_x = 2)+
  ylab('Mahalanobis Distance') + xlab('SNP Number')+ggtitle('Gentian')
#Warning: 3 NA -> 1 at 0.01 MAF

##GGPlot
library(ggplot2)
library(ggthemes)
library(ggrepel)

#alcon
dev.off()
alcloadings <- as.data.frame(alcon2$loadings)
ggplot(data = alcloadings)+
  geom_point(aes(x = V1, y = V2))+
  xlab('PCA1')+ylab('PCA2') +ggtitle('Alcon')+
  theme_classic()+
  geom_point(data = alcloadings[outlieralcon2,],aes(x = V1, y = V2), color = "#116E8A")+
  geom_label_repel(data = alcloadings[outlieralcon2,],aes(x = V1, y = V2, label = rownames(alcloadings[outlieralcon2,])), 
                   color = "#116E8A")
  
gentian
genloadings <- as.data.frame(gentian$loadings)
ggplot(data = genloadings)+
  geom_point(aes(x = V1, y = V2))+
  xlab('PCA1')+ylab('PCA2') +ggtitle('Gentian')+
  theme_classic()+
  geom_point(data = genloadings[outliergentian,],aes(x = V1, y = V2), color = "#116E8A")+
  geom_label_repel(data = genloadings[outliergentian,],aes(x = V1, y = V2, label = outliergentian), 
                   color = "#116E8A")
