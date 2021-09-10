#install.packages('adespatial')
#install.packages('readxl')
library(adespatial)
library(readxl)

#load data
setwd('/home/jonas/Documents/Masterthesis/manuscript/data/')
envdata <- read_excel('Cleaned SNP datasets.xlsx', sheet = 1)
alc.xy <- envdata[,6:5]
envdata$PopID
rownames(alc.xy) <- envdata$PopID
gen.xy <- read.delim('Gentian_coordinates.txt')
gen.xy <- rbind(gen.xy, alc.xy) #actual sampled population localities are at the end of the dataframe (important later)

#check which population is what number
alcdata <- read_excel('Cleaned SNP datasets.xlsx', sheet = 2)
alcpops <- alcdata[1,4:29]
colnames(alcpops) <- 1:26

gendata <- read_excel('Cleaned SNP datasets.xlsx', sheet = 4)
genpops <- gendata[1,4:29]
colnames(genpops) <- 1:26

# Compute the MEMs corresponding to all non-null eigenvalues
?dbmem
#alc.mems.all <- dbmem(alc.xy, MEM.autocor = 'all')
alc.mems <- dbmem(alc.xy, MEM.autocor = "positive", store.listw = T) #threshold = longest edge of minimum distance spanning tree
#gen.mems.all <- dbmem(gen.xy, MEM.autocor = 'all')
gen.mems <- dbmem(gen.xy, MEM.autocor = "positive", store.listw = T)

plot(gen.mems$MEM1[(495-27):495], gen.mems$MEM2[(495-27):495], xlab = 'Eigenvector 1',
     ylab = 'Eigenvector 2', main = 'dbMEM - Gentian')
?plot

#MEM figure
gen.mems.sampled <- gen.mems[(495-27):495,] #last 28 pops are those we sampled
library(ggplot2)
library(ggrepel)
ggplot(aes(x = MEM1, y = MEM2),data = gen.mems.sampled)+
  geom_label_repel(label = rownames(gen.mems.sampled),
                   max.overlaps = 17, size = 2.1)+
  geom_point()+
  ggtitle('dbMEM - Gentian')+ theme_bw()

#we observe that populations with MEM1 > 1 cluster together
rownames(gen.mems.sampled[gen.mems.sampled$MEM1>1,])

ggplot(data = alc.mems)+
  geom_point(aes(x = 1:length(alc.mems$MEM1), y = MEM1))+
  geom_label(mapping = aes(x = 21, y = 1.92), label ='SE Populations')+
  ggtitle('dbMEM - Alcon')+ theme_bw() + xlab('Index')

#plot localities
ggplot(data = alc.xy)+
  geom_point(aes(x = Longitude, y = Latitude))+
  theme_bw()+ ggtitle('Alcon Populations')+
  geom_segment(x = 1.75, y = 43.05, xend = 1.75, yend = 43.075, 
               arrow = arrow(length = unit(0.2, 'cm')))+
  geom_text(x = 1.75, y = 43.085, label = 'N')

ggplot(data = gen.xy)+
  geom_point(aes(x = Longitude, y = Latitude), shape = 20)+
  geom_point(data = alc.xy, aes(x = Longitude, y = Latitude), color = 'darkorange3')+
  theme_bw()+ ggtitle('Gentian Populations')+
  geom_segment(x = 1.9, y = 43.08, xend = 1.9, yend = 43.12, 
               arrow = arrow(length = unit(0.2, 'cm')))+
  geom_text(x = 1.9, y = 43.14, label = 'N')+
  geom_text_repel(data = alc.xy, aes(x = Longitude, y = Latitude), 
                   label = rownames(alc.xy), color = 'darkorange3', size = 2, 
                   max.overlaps = 18)

####CLUSTERING#####
#first cluster all, then check assignement of sampled localities -> problem: all in same cluster
gen.clust <- kmeans(gen.mems, 2) #cluster sampled gentian populations based on MEMs constructed from ALL known poputations
gen.clust.sampled <- gen.clust$cluster[(495-27):495] #only intrested in mems from sampled populations
gen.clust.sampled #all gentian populations in same cluster
#cluster only sampled localities (but on dbMEMs constructed based on all populations)
gen.clust.sampled2 <- kmeans(gen.mems[(495-27):495], 2)
gen.clust.sampled2$cluster #SW3 and SW4 separate

alc.clust <- kmeans(alc.mems, 2)                #cluster alcon MEMs
alc.clust.sampled <- alc.clust$cluster
alc.clust.sampled #SE popuations cluster seperatly

names(gen.clust.sampled[gen.clust.sampled==1])
names(gen.clust.sampled[gen.clust.sampled==2])

#install.packages('ggrepel')
library(ggplot2)
library(ggrepel)

#alc (sampled only)
ggplot(data = alc.xy, aes(x = Longitude, y = Latitude))+
  geom_point()+
  geom_label_repel(label = rownames(alc.xy), max.overlaps = 17, size = 3)+
  theme_bw()

#gen (all localities)
ggplot(data = gen.xy, aes(x = Longitude, y = Latitude))+
  geom_point()+
  geom_label_repel(data = gen.xy[(495-27):495,], label = rownames(gen.xy[(495-27):495,]), max.overlaps = 25, size = 3)+
  theme_bw()

 #not needed, we'll use figure