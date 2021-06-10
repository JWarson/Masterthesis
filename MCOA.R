#install.packages('adegenet')
#install.packages('ade4')
library(adegenet)
alcon <- read.delim("/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/alconSNP", h = T)

row.names(alcon) <- alcon$SampleID
alcon <- alcon[!(alcon$PopID =='Ker'),] #remove ker individuals
popsA <- alcon$PopID
alcon <- subset(alcon, select = -c(SampleID, PopID))

#function to set allele data ready for conversion to genind object
#initialization
dfA <- data.frame(matrix(nrow = 655))
#fill dfA with combination of the 2 columns per snp
for (i in (0:125)){
  dfA <- cbind(dfA, (paste(alcon[,2*i+1], alcon[,2*i+2])))
}
#remove first row (was needed for initialization)
dfA <- dfA[,2:ncol(dfA)]
#change colnames to SNP
c.names <- c()
for (i in 0:125){
  c.names <- append(c.names, colnames(alcon[2*i+1])) #selects every other column name from 'alcon'
}
colnames(dfA) <- c.names
View(dfA) #looks nice!

#conversion to genind
genindA <- df2genind(dfA, sep=" ", ploidy = 2, pop = popsA) #works now
par(mfrow = c(1,1))
#some statistics on the genind object
Asmry <- summary(genindA)
plot(Asmry$Hexp, Asmry$Hobs, main = "Observed vs expected heterozygosity", 
     xlab = "Expected Heterozygosity", ylab = "Observed Heterozygosity")
abline(0, 1, col = "red")
text(Asmry$Hexp[Asmry$Hobs > 0.8], Asmry$Hobs[Asmry$Hobs > 0.8], 
     rownames(Asmry$Hexp[Asmry$Hobs > 0.8]), pos = 2, cex = 0.7, 
     col = '#116E8A') #label outstanding SNPs (SNP1626b, SNP7633)
t.test(Asmry$Hexp, Asmry$Hobs, paired = TRUE, var.equal = TRUE) #is significantly differant -> deviation from HWE

##CO-INERTIA#
#############
#convert to genpop
genpopA <- genind2genpop(genindA)
seplocA <- seploc(genpopA)
AX <- lapply(seplocA, tab)
class(AX)
names(AX)
AX$SNP1043
AX <- lapply(AX, prop.table, 1) #turns allele counts in percentages
length(AX) #126 separate tables, 1 for each locus

#pca on each locus separately
#Error with fixed loci -> remove fixed SNP1600
AX <- AX[names(AX) != "SNP1600"] #pca doesnt work on this fixed allele (is also irrelavant)

A.pca <- list()
for (i in 1:length(AX)){
  A.pca[[length(A.pca)+1]] <- dudi.pca(AX[i],center = TRUE, scale = FALSE, 
                                       scannf = FALSE, nf = 3)
}
names(A.pca) <- names(AX) #add SNP name to each list

#visualize#
###########
#install.packages("remotes")
#remotes::install_github("JosephCrispell/basicPlotteR")
library(basicPlotteR)
pca.snp <- function(pca, SNP){
  x <- pca[[SNP]]
  if (length(x$li)>1){
    s.label(x$li)
    add.scatter.eig(x$eig, 3, 1, 2)
    title(paste("PCA 1 & 2 of", SNP))
  }
  else{
    plot(A.pca$SNP7633$li, main = paste("PCA 1 of",SNP))
    addTextLabels(x = x$li$Axis1, y= c(rep(1, times = nrow(x$li))), 
                    rownames(x$li), cex.label = 0.7, col.label = "#116E8A")
  }
}

pca.snp(pca = A.pca, SNP = 'SNP1626b')
pca.snp(pca = A.pca, SNP = 'SNP7633')

#MCOA#
######
A.ktab <- ktab.list.dudi(A.pca)
A.mcoa1 <- mcoa(A.ktab)


mcoa.snp <- function(SNP) {
  newCoord <- split(A.mcoa1$Tli, A.mcoa1$TL[, 1])
  names(newCoord) <- names(AX)
  s.label(newCoord[SNP], label = gsub(".SNP.*","",rownames(newCoord[[SNP]])), 
          xax = 1, yax = 2, sub = SNP,csub = 1.2)
  }
mcoa.snp(SNP = "SNP4388")
mcoa.snp(SNP = "SNP8192")

#mcoa of alcon SNPs identified as outlier in all analysis
mcoa.snp(SNP ="SNP7058") #98
mcoa.snp(SNP ="SNP7648") #109 #definitely interesting one, SW6, SW5, SW3 show high distinction from others for this SNP

#nicer ggplot
library(ggrepel)
newCoord <- split(A.mcoa1$Tli, A.mcoa1$TL[, 1])
mcoa109 <- newCoord[["SNP7648"]]
ggplot(mcoa109, aes(x=Axis1, y = Axis2))+
  geom_label(label = gsub(".SNP.*","",rownames(mcoa109)), 
             label.r = unit(0.3, "lines"), 
             label.padding = unit(0.2, "lines"),
             color = '#00407A', size = 3.5)+
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) +
  theme_classic()

contrib <- A.mcoa1$cov2[,1]/sum(A.mcoa1$cov2[,1])
names(contrib) <- rownames(A.mcoa1$cov2)
barplot(contrib)
contrib[contrib>0.025] #1937, 3745, 4855, 7381 and 6215 contribute most to differentiation
contrib[["SNP7648"]] #1.3 percent for 109

#
##
### Gentian
##  -------
#

gentian <- read.delim("/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/gentianSNP", h = T)

row.names(gentian) <- gentian$SampleID
popsA <- gentian$PopID
gentian <- subset(gentian, select = -c(SampleID, PopID))

#function to set allele data ready for conversion to genind object
#initialization
dfG <- data.frame(matrix(nrow = 652))
#fill dfG with combination of the 2 columns per snp
for (i in (0:103)){
  dfG <- cbind(dfG, (paste(gentian[,2*i+1], gentian[,2*i+2])))
}
#remove first row (was needed for initialization)
dfG <- dfG[,2:ncol(dfG)]
#change colnames to SNP
c.names <- c()
for (i in 0:103){
  c.names <- append(c.names, colnames(gentian[2*i+1])) #selects every other column name from 'gentian'
}
colnames(dfG) <- c.names
View(dfG) #looks good

#conversion to genind

genindG <- df2genind(dfG, sep=" ", ploidy = 2, pop = popsA) #works now

#some statistics on the genind object
Gsmry <- summary(genindG)
plot(Gsmry$Hexp, Gsmry$Hobs, main = "Gentian: Observed vs expected heterozygosity", 
     xlab = "Expected Heterozygosity", ylab = "Observed Heterozygosity")
abline(0, 1, col = "red")
text(Gsmry$Hexp[Gsmry$Hobs > 0.8], Gsmry$Hobs[Gsmry$Hobs > 0.8], 
     rownames(Gsmry$Hexp[Gsmry$Hobs > 0.8]), pos = 2, cex = 0.7, 
     col = 'darkblue') #label outstanding SNPs (SNP452, SNP3460, SNP 2099)
t.test(Gsmry$Hexp, Gsmry$Hobs, paired = TRUE, var.equal = TRUE) #is significantly differant -> deviation from HWE

##CO-INERTIA#
#############
#convert to genpop
genpopG <- genind2genpop(genindG)
seplocG <- seploc(genpopG)
GX <- lapply(seplocG, tab)
class(GX)
names(GX)
GX$SNP762
GX <- lapply(GX, prop.table, 1) #turns allele counts in percentages
length(GX) #104 separate tables, 1 for each locus

#pca on each locus separately
#Error with fixed loci -> remove fixed SNP1600
#GX <- GX[names(GX) != "SNP1600"] #pca doesnt work on this fixed allele (is also irrelavant)

G.pca <- list()
for (i in 1:length(GX)){
  G.pca[[length(G.pca)+1]] <- dudi.pca(GX[i],center = TRUE, scale = FALSE, 
                                       scannf = FALSE, nf = 3)
}
names(G.pca) <- names(GX) #add SNP name to each list


#visualize#
###########
#install.packages("remotes")
#remotes::install_github("JosephCrispell/basicPlotteR")
library(basicPlotteR)

pca.snp <- function(pca, SNP){
  x <- pca[[SNP]]
  if (length(x$li)>1){
    s.label(x$li)
    add.scatter.eig(x$eig, 3, 1, 2)
    title(paste("First 2 PC-GXes of", SNP))
  }
  else{
    plot(G.pca$SNP7633$li, main = paste("Single PC-GXis of",SNP))
    addTextLabels(x = x$li$GXis1, y= c(rep(1, times = nrow(x$li))), 
                  rownames(x$li), cex.label = 0.7, col.label = "cornsilk4")
  }
}

pca.snp(pca = G.pca, SNP = 'SNP1010')
pca.snp(pca = G.pca, SNP = 'SNP1142')

#MCOA#
######
G.ktab <- ktab.list.dudi(G.pca)
G.mcoa1 <- mcoa(G.ktab)
G.mcoa1

mcoa.snp <- function(SNP) {
  newCoord <- split(G.mcoa1$Tli, G.mcoa1$TL[, 1])
  names(newCoord) <- names(GX)
  s.label(newCoord[SNP], xGX = 1, yGX = 2, sub = SNP,csub = 1.5)
}
mcoa.snp(SNP = "SNP1142")#not working

contrib <- G.mcoa1$cov2[,1]/sum(G.mcoa1$cov2[,1])
names(contrib) <- rownames(G.mcoa1$cov2)
barplot(contrib)
contrib[contrib>0.030] #SNP3080  &  SNP3932  contribute most to differentiation
