#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("LEA")
library(LEA)
setwd("/home/jonas/Documents/Masterthesis/data/Analyses/PCadapt")

##ENV VARIABLES####
#each variable needed in separate list, with a value for each individual
alc.pops <- read.csv('popsAlc', sep = '\t', h = T) #contains a row with the popID for each individual
a.indv.ppop <- table(alc.pops) #amount of individuals in each population
a.indv.ppop <- a.indv.ppop[2:length(a.indv.ppop)] #removes ker location

env <- read.csv("/home/jonas/Documents/Masterthesis/data/Analyses/Vegan/environment", h = T, sep = '\t')
env.noKer <- env[2:nrow(env),] #remove ker location

#alcon

env.alc <- env.noKer[!(env.noKer$PopID=="SE2") & !(env.noKer$PopID=="SW4"),] #no genotypic data for these pop
length(a.indv.ppop) == nrow(env.alc) #populations match
colnames(env.alc)

##Altitude
Aalt <-rep(env.alc[,7], times = as.vector(a.indv.ppop))
#last 6 indv. in .lfmm file are from ker location -> append these (at end)
Aalt <- append(Aalt, rep(env[env$PopID=="?", 7], times = 6))
write.table(Aalt, "Aalt.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

##Regional connectivity
Aconn <- rep(env.alc[,10], times = as.vector(a.indv.ppop))
#last 6 indv. in .lfmm file are from ker location -> append these (at end)
Aconn <- append(Aconn, rep(env[env$PopID=="?", 10], times = 6))
write.table(Aconn, "Aconn.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

##Patch Size
Aps <- rep(env.alc[,8], times = as.vector(a.indv.ppop))
#last 6 indv. in .lfmm file are from ker location -> append these (at end)
Aps <- append(Aps, rep(env[env$PopID=="?", 8], times = 6))
write.table(Aps, "Aps.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

##Suitability
As <- rep(env.alc[,13], times = as.vector(a.indv.ppop))
#last 6 indv. in .lfmm file are from ker location -> append these (at end)
As <- append(As, rep(env[env$PopID=="?", 13], times = 6))
write.table(As, "As.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)
#Gentian
#(see further)

####LFMMs######
###############

#Alcon

#estimate latent factors (ancestral populations)
alconGT <- lfmm2geno('Alcon.lfmm')
alc.snmf <- snmf(alconGT, K = 1:18, entropy = T, ploidy = 2, project="new")
par(mfrow=c(1,1))
plot(alc.snmf) #14
barplot(t(Q(alc.snmf, K = 13)), col = 1:13) #not working

##Altitude
Aalt.lfmm = lfmm("Alcon.lfmm", "Aalt.env", K = 14, rep = 5, project="new")
#The zscores:
zs = z.scores(Aalt.lfmm, K = 13)
zs.Aalt = apply(zs, MARGIN = 1, median)
lambda = median(zs.Aalt^2)/qchisq(0.5, df = 1) #very high GIF (12.2) -> expected because these are under expected to be under selection
adjP.Aalt = pchisq(zs.Aalt^2/lambda, df = 1, lower = FALSE)
hist(adjP.Aalt, col = "red") #not good -> GIF overly conservative
adjP.Aalt = pchisq(zs.Aalt^2/8, df = 1, lower = FALSE)
hist(adjP.Aalt, col = "lightblue", main = "Alcon-Altitude adj. P-values")
## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 122
#fdr level q
q = 0.1
w = which(sort(adjP.Aalt) < q * (1:L)/L)
cand.bh.Aalt = order(adjP.Aalt)[w] #98, 35, 79, 33 

##Regional Connectivity
Aconn.lfmm = lfmm("Alcon.lfmm", "Aconn.env", K = 14, rep = 5, project = "new")

zs = z.scores(Aconn.lfmm, K = 14)
zs.Aconn = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Aconn^2)/qchisq(0.5, df = 1) #high GIF (12.6)
adjP.Aconn = pchisq(zs.Aconn^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Aconn, col = "red") #not good
adjP.Aconn = pchisq(zs.Aconn^2/7, df = 1, lower = FALSE) #GIF overconservative
hist(adjP.Aconn, col = "lightblue", main = "Alcon-RegionalConnectivity adj. P-values") #better
#FDR control
w = which(sort(adjP.Aconn) < q * (1:L)/L)
cand.bh.Aconn = order(adjP.Aconn)[w] #98, 79, 35, 7, 108

##Patch Size
Aps.lfmm = lfmm("Alcon.lfmm", "Aps.env", K =14, rep = 5, project = "new")

zs = z.scores(Aps.lfmm, K = 14)
zs.Aps = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Aps^2)/qchisq(0.5, df = 1) #GIF = 7.66
adjP.Aps = pchisq(zs.Aps^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Aps, col = "red") #not good
adjP.Aps = pchisq(zs.Aps^2/5, df = 1, lower = FALSE)
hist(adjP.Aps, col = "lightblue", main = "Alcon-RegionalConnectivity adj. P-values")
#FDR control
w = which(sort(adjP.Aps) < q * (1:L)/L)
cand.bh.Aps = order(adjP.Aps)[w] #110 108  82  62 102  18  60  41  30  29 101  79
#-> double check environment: error at write.table-> Aalt was always written instead of Aconn, Aps or As

##Suitability
As.lfmm = lfmm("Alcon.lfmm", "As.env", K =14, rep = 5, project = "new")

zs = z.scores(As.lfmm, K = 14)
zs.As = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.As^2)/qchisq(0.5, df = 1) #GIF = 9.05
adjP.As = pchisq(zs.As^2/lambda2, df = 1, lower = FALSE)
hist(adjP.As, col = "red") #not good
adjP.As = pchisq(zs.As^2/7, df = 1, lower = FALSE)
hist(adjP.As, col = "lightblue", main = "Alcon-RegionalConnectivity adj. P-values")
#FDR control
w = which(sort(adjP.As) < q * (1:L)/L)
cand.bh.As = order(adjP.As)[w] #110 102  41  56  18 108  62  30  60  21

#########
#Gentian# -> same environmental variables but genetic data comes from different populations
#########

#estimate latent factors (ancestral populations)
gentianGT <- lfmm2geno('Gentian.lfmm')
gen.snmf <- snmf(gentianGT, K = 1:18, entropy = T, ploidy = 2, project="new")
plot(gen.snmf) #5
barplot(t(Q(alc.snmf, K = 5)), col = 1:5) #not working

#restructure env matrix based on gentian samples per population
gen.pops <- read.csv("popsGen", h = T, sep = "\t") #652 individuals (from "cleaned snp datasets", sheet "gentian SNP")
g.indv.ppop <- table(gen.pops)
setdiff(env.noKer$PopID, names(g.indv.ppop)) #NW18 must be deleted
env.gen <- env.noKer[!(env.noKer$PopID=="NW18"),]
setdiff(env.gen$PopID, names(g.indv.ppop)) #match

colnames(env.gen)

##ALTITUDE
Galt <- rep(env.gen[,7], times = as.vector(g.indv.ppop))
write.table(Galt, "Galt.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

Galt.lfmm <- lfmm("Gentian.lfmm", "Galt.env", K = 5, rep = 5, project = 'new')
#significant SNPs
zs = z.scores(Galt.lfmm, K = 5)
zs.Galt = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Galt^2)/qchisq(0.5, df = 1) #high GIF (3.0)
adjP.Galt = pchisq(zs.Galt^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Galt, col = "red") #not good
adjP.Galt = pchisq(zs.Galt^2/2.4, df = 1, lower = FALSE)
hist(adjP.Galt, col = "lightblue", main = "Gentian-Altitude adj. P-values") #OK
#FDR control
L = 105
w = which(sort(adjP.Galt) < q * (1:L)/L)
cand.bh.Galt = order(adjP.Galt)[w] #23, 39, 82

##TOTAL CONNECTIVITY
Gconn <- rep(env.gen[,9], times = as.vector(g.indv.ppop))
write.table(Gconn, "Gconn.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

Gconn.lfmm <- lfmm("Gentian.lfmm", "Gconn.env", K = 5, rep = 5, project = 'new')
#significant SNPs
zs = z.scores(Gconn.lfmm, K = 5)
zs.Gconn = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Gconn^2)/qchisq(0.5, df = 1) #GIF = 1.11
adjP.Gconn = pchisq(zs.Gconn^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Gconn, col = "red") #not good
adjP.Gconn = pchisq(zs.Gconn^2/1, df = 1, lower = FALSE) #GIF overconservative
hist(adjP.Gconn, col = "lightblue", main = "Gentian-Total Connectivity adj. P-values") #OK
#FDR control
L = 105
w = which(sort(adjP.Gconn) < q * (1:L)/L)
cand.bh.Gconn = order(adjP.Gconn)[w] #no candidates

##PATCH SIZE
Gps <- rep(env.gen[,8], times = as.vector(g.indv.ppop))
write.table(Gps, "Gps.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

Gps.lfmm <- lfmm("Gentian.lfmm", "Gps.env", K = 5, rep = 5, project = 'new')
#significant SNPs
zs = z.scores(Gps.lfmm, K = 5)
zs.Gps = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Gps^2)/qchisq(0.5, df = 1) #pretty high GIF = 1.92
adjP.Gps = pchisq(zs.Gps^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Gps, col = "red") #not good
adjP.Gps = pchisq(zs.Gps^2/1, df = 1, lower = FALSE) #GIF overconservative
hist(adjP.Gps, col = "lightblue", main = "Gentian-Total Connectivity adj. P-values") #OK
#FDR control
L = 105
w = which(sort(adjP.Gps) < q * (1:L)/L)
cand.bh.Gps = order(adjP.Gps)[w] #67, 74, 56, 61, 101, 82

##SUITABILITY
Gs <- rep(env.gen[,13], times = as.vector(g.indv.ppop))
write.table(Gs, "Gs.env", append = FALSE, sep = " ", dec = ".",
            row.names = F, col.names = F)

Gs.lfmm <- lfmm("Gentian.lfmm", "Gs.env", K = 5, rep = 5, project = 'new')
#significant SNPs
zs = z.scores(Gs.lfmm, K = 5)
zs.Gs = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Gs^2)/qchisq(0.5, df = 1) #GIF = 1.93
adjP.Gs = pchisq(zs.Gs^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Gs, col = "red") #not good
adjP.Gs = pchisq(zs.Gs^2/1, df = 1, lower = FALSE) #GIF overconservative
hist(adjP.Gs, col = "lightblue", main = "Gentian-Total Connectivity adj. P-values") #OK
#FDR control
w = which(sort(adjP.Gs) < q * (1:L)/L)
cand.bh.Gs = order(adjP.Gs)[w] #67,62


