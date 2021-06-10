setwd('/home/jonas/Documents/Masterthesis/data/Analyses/LFMM')
library(LEA)

alconGT <- lfmm2geno('AlcLFMM.lfmm')
alc.snmf <- snmf(alconGT, K = 1:18, entropy = T, ploidy = 2, project="new")
par(mfrow=c(1,1))
plot(alc.snmf) #7
barplot(t(Q(alc.snmf, K = 18)), col = 1:18) #not working

##Altitude
?lfmm
Aalt.lfmm = lfmm("AlcLFMM.lfmm", "Aalt.env", K = 7, rep = 5, project="new")
#The zscores:
zs = z.scores(Aalt.lfmm, K = 7)
zs.Aalt = apply(zs, MARGIN = 1, median)
lambda = median(zs.Aalt^2)/qchisq(0.5, df = 1) #very high GIF (9.9) -> expected because these are under expected to be under selection
adjP.Aalt = pchisq(zs.Aalt^2/lambda, df = 1, lower = FALSE)
hist(adjP.Aalt, col = "red") #not good -> GIF overly conservative
adjP.Aalt = pchisq(zs.Aalt^2/4, df = 1, lower = FALSE)
hist(adjP.Aalt, col = "lightblue", main = "Alcon-Altitude adj. P-values")
## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 127
#fdr level q
q = 0.1
w = which(sort(adjP.Aalt) < q * (1:L)/L)
cand.bh.Aalt = order(adjP.Aalt)[w] #98, 35, 79, 33 

##Regional Connectivity
Aconn.lfmm = lfmm("AlcLFMM.lfmm", "Aconn.env", K = 7, rep = 5, project = "new")

zs = z.scores(Aconn.lfmm, K = 7)
zs.Aconn = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Aconn^2)/qchisq(0.5, df = 1) #high GIF (10.6)
adjP.Aconn = pchisq(zs.Aconn^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Aconn, col = "red") #not good
adjP.Aconn = pchisq(zs.Aconn^2/5.5, df = 1, lower = FALSE) #GIF overconservative
hist(adjP.Aconn, col = "lightblue", main = "Alcon-RegionalConnectivity adj. P-values") #better
#FDR control
w = which(sort(adjP.Aconn) < q * (1:L)/L)
cand.bh.Aconn = order(adjP.Aconn)[w] #

##Suitability
As.lfmm = lfmm("AlcLFMM.lfmm", "As.env", K =7, rep = 5, project = "new")

zs = z.scores(As.lfmm, K = 7)
zs.As = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.As^2)/qchisq(0.5, df = 1) #GIF = 4.95
adjP.As = pchisq(zs.As^2/lambda2, df = 1, lower = FALSE)
hist(adjP.As, col = "red") #not good
adjP.As = pchisq(zs.As^2/3.5, df = 1, lower = FALSE)
hist(adjP.As, col = "lightblue", main = "Alcon-RegionalConnectivity adj. P-values")
#FDR control
w = which(sort(adjP.As) < q * (1:L)/L)
cand.bh.As = order(adjP.As)[w] #111 109  19  31  35 103  51  61  46  52  38 102  45  25  84 126  30

##Patch Size
Aps.lfmm = lfmm("AlcLFMM.lfmm", "Aps.env", K =7, rep = 5, project = "new")

zs = z.scores(Aps.lfmm, K = 7)
zs.Aps = apply(zs, MARGIN = 1, median)
lambda2 = median(zs.Aps^2)/qchisq(0.5, df = 1) #GIF = 3.4
adjP.Aps = pchisq(zs.Aps^2/lambda2, df = 1, lower = FALSE)
hist(adjP.Aps, col = "red") #not good
adjP.Aps = pchisq(zs.Aps^2/2.5, df = 1, lower = FALSE)
hist(adjP.Aps, col = "lightblue", main = "Alcon-RegionalConnectivity adj. P-values")
#FDR control
w = which(sort(adjP.Aps) < q * (1:L)/L)
cand.bh.Aps = order(adjP.Aps)[w] #

