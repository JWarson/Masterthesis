if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
devtools::install_github("gaospecial/ggVennDiagram")
library(ggVennDiagram)

AlcOutl <- list(
  PCA = c(30,34,60,108),
  Bayescan = c(7,25,29,98,99,108,109,116)
)
ggVennDiagram(AlcOutl) #108 is common

GenOutl <- list(
  PCA = c(14,25,27,34,36,55,63,76,77,94,98,105),
  Bayescan = c(5,21,23,29,32,53,59,72,73,77,80,84,87,88,97)
)
ggVennDiagram(GenOutl) #77 is common

##LFMM
lfmm.alc <- read.csv('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/alcon.lfmm.signSNPs', h = T, sep = '\t')
lfmm.gen <- read.csv('/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/gentian.lfmm.signSNPs', h = T, sep = '\t')
lfmm.alc2 <- read.csv('/home/jonas/Documents/Masterthesis/data/Analyses/LFMM/alconlfmm127', h = T, sep ='\t' )

dfA <- list(Altitude = lfmm.alc2[,1],Connectivity = lfmm.alc2[,2], Suitabliltiy = lfmm.alc2[,3], PatchSize = lfmm.alc2[,4])
dfA <- lapply(dfA, function(x) x[!is.na(x)])

dfG <- list(Altitude = lfmm.gen[,1],Connectivity = lfmm.gen[,2], Suitabliltiy = lfmm.gen[,3], PatchSize = lfmm.gen[,4])
dfG <- lapply(dfG, function(x) x[!is.na(x)])
?ggVennDiagram
ggVennDiagram(dfA, label = "percent", label_alpha = 0, label_size = 2.5, set_size = 3,
              category.names= c('Altitude', 'Connectivity', 'Suitability', 'Patch Size   '))+
  ggplot2::scale_fill_gradient(low="white",high = "#116E8A")+
  ggtitle('Alcon')+
  theme(legend.position = 'none',  plot.margin=unit(c(0,0.3,0,0),"cm"))

ggVennDiagram(dfG[c(1,3,4)], label = "both", label_alpha = 0, na.rm = T)+
  ggplot2::scale_fill_gradient(low="white",high = "#116E8A")+
  ggtitle('Gentian')+
  theme(legend.position = 'none')


##Search shared adaptive SNPs across methods

#combine to 1 vector 
aUL <- unique(unname(unlist(dfA)))
gUL <- unique(unname(unlist(dfG)))

outl <- readxl::read_excel('/home/jonas/Documents/Masterthesis/data/results.xlsx', 
                   sheet= 1, col_names = T)

compAlc <- list(Bayescan = as.integer(outl$Bayescan[2:9]), 
                PCA = as.integer(outl$PCA[2:10]), 
                LFMM = aUL)

compGen <- list(Bayescan = as.integer(outl$...3[2:16]), 
                PCA = as.integer(outl$...5[2:13]), 
                LFMM = gUL)

#52BDEC  #kul light blue
#00407A #kul dark blue
#116E8A #darker kul theme color



#Venn diagram
ggVennDiagram(compAlc, label = "both", label_alpha = 0, fill = 'black')+
  ggplot2::scale_fill_gradient(low="white",high = "#116E8A")+
  theme(legend.position = 'none',  plot.margin=unit(c(0,0.3,0,0),"cm"))

ggVennDiagram(compGen, label = "both", label_alpha = 0, )+
  ggplot2::scale_fill_gradient(low="white",high = "#116E8A")+
  theme(legend.position = 'none')

Reduce(intersect, compAlc) #98 SNP7058 and 109 SNP7648
Reduce(intersect, compGen) #0
