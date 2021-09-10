#RDA van SNPs (response matrix) vs. MEM variabele(n) --> R²adjusted om 
#isolation-by-distance (ibd) te bepalen, en een RDA biplot van de populaties 
#op te bekijken of er populaties zijn die genetisch geïsoleerd zijn door ibd

#load files
setwd('/home/jonas/Documents/Masterthesis/manuscript/data/rdaobjects')
alc.mems <- read.delim(file = 'alcmems', sep = ' ')
gen.mems <- read.delim(file = 'genmems', sep = ' ')
alconSNP <- read.delim('/home/jonas/Documents/Masterthesis/manuscript/data/alconMAF.txt', sep = '\t')
gentianSNP <- read.delim('/home/jonas/Documents/Masterthesis/manuscript/data/gentianMAF.txt', sep = '\t')

alc.mems <- alc.mems[rownames(alc.mems) %in% rownames(alconSNP),]
gen.mems.sampled <- gen.mems[rownames(gen.mems) %in% rownames(gentianSNP),]
library(vegan)

#hellinger-transformed snp data
alconSNP.hel <- decostand(alconSNP, 'hellinger')
gentianSNP.hel <- decostand(gentianSNP, 'hellinger')

#perform rda
alcon.rda <- rda(alconSNP.hel ~alc.mems,
    scale = F) #https://cran.r-project.org/web/packages/vegan/vignettes/decision-vegan.pdf
plot(alcon.rda, main = 'Alcon RDA', 
     xlab = paste('RDA1', '(',round(100*(alcon.rda$CCA$eig/alcon.rda$tot.chi), 
                                    digits = 2),'%)'))
RsquareAdj(alcon.rda) #14,1%
anova.cca(alcon.rda, step = 1000) #model is significant
anova.cca(alcon.rda, by = 'axis', step =1000)
alcon.rda$CCA$wa

sites.alc <- as.data.frame(scores(alcon.rda, display = 'wa'))
plot(sites.alc,  pch =19, col = 'lightblue', cex=0.5,
     main = 'Alcon Sites Scores')
text(sites.alc[,2] ~sites.alc[,1], labels=rownames(sites.alc), cex=0.6, font=2)

#create plot RDA plot with total connectivity
library(ggplot2)
library(ggrepel)
env <- read.csv('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Vegan/environment', h = T, sep = '\t')
sites.alc$T.C. <- env$TotalConnectivity[env$PopID %in% rownames(sites.alc)]
ggplot(data = as.data.frame(sites.alc))+
  geom_point(aes(x= RDA1, y = PC1, colour = T.C.))+
  scale_color_gradient(low = 'tomato1', high = 'springgreen')+
  theme_bw()+ggtitle('Alcon Sites with Total Connectivity')+
  geom_label_repel(aes(x= RDA1, y = PC1), label = rownames(sites.alc), size = 2)

#gentian.rda <- rda(gentianSNP ~ .,data = gen.mems.sampled[,1:10], scale = T)
#plot(gentian.rda, xlab = 'RDA1 (23.7%)', ylab = 'RDA2 (12.5%)')
#RsquareAdj(gentian.rda) #r² = 1 with all mems, only first 10 MEMs: R²adj = 41%

#forward selection on all mems
#ord.spa = rda(gentianSNP.hel ~ ., data=gen.mems.sampled, scale= FALSE)
#(R2.all.spa = RsquareAdj(ord.spa)$adj.r.squared)
#stp.spa = ordistep(rda(gentianSNP.hel ~ 1, data=gen.mems.sampled), scope = formula(ord.spa), scale= FALSE, direction="forward", pstep = 1000)
#(selected.spa = attributes(stp.spa$terms)$term.labels) #7,6,2,18,3,5,29,13

sel.mems <-gen.mems.sampled[,c(7,6,2,18,3,5,29,13)]

rda_test = rda(gentianSNP.hel~., data = sel.mems,scale=F) #better for testing
overall.test = anova.cca(rda_test, step=1000)
(axis.test = anova.cca(rda_test, by="axis", step= 1000))
(RsquareAdj(rda_test)) #53,4%

#try ordinate all sites (495) according to rda ordination
rda_plot = rda(gentianSNP.hel~., data = sel.mems,scale=T)
plot(rda_plot, scaling = 2)
plot(rda_plot$CCA$wa[,1:2])

?scores
reg <- scores(rda_plot, choices = c(1,2), display = c("reg"))
sites <- scores(rda_plot, choices = c(1,2), display = c("wa"))
gen.mems.sel<- gen.mems[,c(7,6,2,18,3,5,29,13)]

gen.mems.plotted <- as.matrix(gen.mems.sel) %*% reg

plot(sites, txt = rownames(sites), pch =19, col = 'lightblue', cex=0.5,
     main = 'Gentian Sites Scores (sampled only)')
text(sites[,2] ~sites[,1], labels=rownames(sites), cex=0.6, font=2)
?plot

library(ggplot2)
library(ggrepel)
ggplot(aes(x = RDA1, y = RDA2), data = as.data.frame(gen.mems.plotted))+
  geom_point(pch = 4)+
  geom_point(data = as.data.frame(gen.mems.plotted[(495-25):495,]), 
             col = 'darkorange1', pch =4)+
  geom_label_repel(data = as.data.frame(gen.mems.plotted[(495-25):495,]), 
             col = 'darkorange1', label = rownames(gen.mems.sampled), size = 2, max.overlaps = 19)+
  theme_bw()+ggtitle('Gentian Sites Scores (all pops)')

ggplot(aes(x = RDA1, y = RDA2), data = as.data.frame(gen.mems.plotted))+
  geom_point(pch = 4)+
  geom_point(data = as.data.frame(gen.mems.plotted[(495-25):495,]), 
                   col = 'darkorange1', pch =4)+
  theme_bw()+ggtitle('All Gentian Populations plotted on RDA axes')

#calculate overall fst
bstats = basic.stats(dat)
global_fst = bstats$overall[8]
print(bstats$overall)
#set some plotting variables
perc.rda = sum(rda_test$CCA$eig) / rda_test$tot.chi
perc.tot = global_fst * perc.rda
all_perc = rda_test$CCA$eig/rda_test$tot.chi

main =sprintf("Gentian AF ~ MEM (%.1f %%; Fst = %.2f; p = %s)", 100*perc.rda, perc.tot, overall.test[1,4]);
xlab=sprintf("RDA1 (%.0f%%, p = %s)", 100*all_perc[1], axis.test[[4]][1])
ylab=sprintf("RDA2 (%.0f%%, p = %s)", 100*all_perc[2], axis.test[[4]][2])
lim = c(-1.2,1.2)

#now draw plot
plot(rda_test, main = main, xlab = xlab, ylab = ylab)

?scores

scores(rda_test, choices = c(1,2), display = c("reg"),
       scaling = "species")

dat <- scores(rda_test, choices = c(1,2), display = c("sp","wa","cn"),
               scaling = "species")
dat <- as.data.frame(dat$sites)

ggplot(aes(x = RDA1, y=RDA2),data = dat)+
  geom_label_repel( label = rownames(dat), size = 2.5,
                   max.overlaps = 17)+
  geom_point()+
  theme_bw()+
  xlab(xlab)+ylab(ylab)+ggtitle('Gentian sites distribution over RDA axes')


##onlinelibrary.wiley.com/doi/full/10.1111/mec.13243

#check whether alcon RDA-MEM data matches connectivity data
env <- read.csv('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Vegan/environment', h = T, sep = '\t')
alc.scores <- scores(alcon.rda, display = 'sites')
rownames(alc.scores) %in% env$PopID
envalc <- env[env$PopID %in% rownames(alc.scores),] #removes missing populations from env
rownames(alc.scores) == envalc$PopID

isol.pops <- c('SW5','SW6','NW19','NE7','NE1','SE4','SE5','SE1','SE3','NW1','NW5')
envalc$RDAConn <- ifelse(envalc$PopID %in% isol.pops, 'Isolated','Connected')

test.conn <- lm(RegionalConnectivity~RDAConn, data = envalc)
summary(test.conn)
boxplot(RegionalConnectivity~RDAConn, data = envalc, xlab = ' ', 
        main ='Alcon Connectivity (p = 0.0182)')

test.conn2 <- lm(TotalConnectivity~RDAConn, data = envalc)
summary(test.conn2)
boxplot(TotalConnectivity~RDAConn, data = envalc, xlab = ' ', 
        main ='Alcon Connectivity (p = 0.0248)')
