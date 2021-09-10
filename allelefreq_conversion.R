alcon <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Arlequin/AlconAlleleFreqs.txt', h = T, sep = "\t")
alcon <- subset(alcon[,1:(ncol(alcon)-1)]) #last column are all NAs
popnamesA <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Arlequin/Alcon4group.res/PopNamesA.txt', h =F)
popnamesA$V3 <- 1:nrow(popnamesA) #numeric popnames, for RDA later
colnames(alcon) <- c("Locus", "Allele", popnamesA[,2]) #for data handling

#remove '?' allele
alconMAF <- subset(alcon, alcon$Allele != "?") #removes 108 '?' allele (zero anyway)

#remove 'allele' column and 'ker' location (last column = ncol(alconMAF))
alconMAF <- subset(alconMAF, select= -c(Allele, Ker))

#remove major allele
library(data.table)
alconMAF$sumAF <- rowSums(alconMAF[,2:ncol(alconMAF)]) #creates column with sum of the allele frequencies at each row
alconMAF <- setDT(alconMAF)[, .SD[which.min(sumAF)], by = Locus] #selects the row with smalles sumAF

#remove sumAF column
alconMAF <- subset(alconMAF, select = -c(sumAF))
#turn around rows and column, 'locus' is header
AlconAF <- as.data.frame(t(alconMAF))
names(AlconAF)<- AlconAF[1,]
AlconAF <- AlconAF[2:nrow(AlconAF),]
str(AlconAF)
AlconAF <- type.convert(AlconAF, as.is =T) #convert to numerical

write.table(AlconAF, file ='/home/jonas/Documents/Masterthesis/manuscript/data/alconMAF.txt',
            row.names = T, col.names = T, sep = '\t')

#gentian#
#########
gentian <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Arlequin/GentianAlleleFreqs.txt', h = F)
gentian <- subset(gentian[,1:(ncol(gentian)-1)]) #last column NA
popnamesG <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Arlequin/Gentian.res/PopNamesG.txt', h =F)
popnamesG$V3 <- 1:nrow(popnamesG) #numeric popnames
colnames(gentian) <- c("Locus", "Allele", popnamesG[,2]) #colnames for data handling

#remove '?' allele
gentianMAF <- subset(gentian, gentian$Allele != "?") #removes 108 '?' alleles
#remove 'allele' column
gentianMAF <- subset(gentianMAF, select= -c(Allele))

#remove major allele
gentianMAF$sumAF <- rowSums(gentianMAF[,2:ncol(gentianMAF)])
gentianMAF <- setDT(gentianMAF)[, .SD[which.min(sumAF)], by = Locus]

#remove sumAF 
gentianMAF <- subset(gentianMAF, select = -c(sumAF))
#turn around rows and column, 'locus' is header
GentianAF <- as.data.frame(t(gentianMAF))
names(GentianAF)<- GentianAF[1,]
GentianAF <- GentianAF[2:nrow(GentianAF),]
GentianAF <- type.convert(GentianAF, as.is =T)
str(GentianAF)

write.table(GentianAF, file ='/home/jonas/Documents/Masterthesis/manuscript/data/gentianMAF.txt',
            row.names = T, col.names = T, sep = '\t')
