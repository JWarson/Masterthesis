# 	 This file is used to plot figures for the software Bayescan in R.

#    This program, BayeScan, aims at detecting genetics markers under selection,
#	 based on allele frequency differences between population. 
#    Copyright (C) 2010  Matthieu Foll
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Arguments:
# - file is the name of your file ex: "output_fst.txt"
# - the q-value threshold corresponding to the target False Discovery Rate (FDR)
# - size is the size of the points and text labels for outliers
# - pos is the distance between the points and the labels 
# - highlight is a optional list of marker indices to display in red.
# - name_highlighted alows to write the indices of highlighted markers instead of using a point like the other markers
# - add_text adds the indices of the outlier markers

# Output:
# This function returns different paremeters in a list
# - outliers: the list of outliers
# - nb_outliers: the number of outliers

# Typical usage: 
# - load this file into R (file/source R code)
# - in R, go to the directory where "output_fst.txt" is (file/change current dir)
# - at the R prompt, type 
# > plot_bayescan("output_fst.txt",0,FDR=0.05)
# if you save the output in a variable, you can recall the different results:
# results<-plot_bayescan("output_fst.txt",0,FDR=0.05)
# results$outliers
# results$nb_outliers

#
# plotting posterior distribution is very easy in R with the output of BayeScan:
# first load the output file *.sel produced by BayeScan
# > mydata=read.table("bi.sel",colClasses="numeric")
# choose the parameter you want to plot by setting for example:
# parameter="Fst1"
# then this line will make the plot for:
# > plot(density(mydata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
# you can plot population specific Fst coefficient by setting
# parameter="Fst1"
# if you have non-codominant data you can plot posterior for Fis coefficients in each population:
# parameter="Fis1"
# if you test for selection, you can plot the posterior for alpha coefficient for selection:
# parameter="alpha1"
# you also have access to the likelihood with:
# parameter="logL"
# if you have the package "boa" installed, you can very easily obtain Highest Probability 
# Density Interval (HPDI) for your parameter of interest (example for the 95% interval):
# > boa.hpd(mydata[[parameter]],0.05)


plot_bayescan<-function(res,FDR=0.05,size=1,pos=0.35,highlight=NULL,name_highlighted=F,add_text=T)
{
  if (is.character(res))
    res=read.table(res)
  
  colfstat=5
  colq=colfstat-2
  
  highlight_rows=which(is.element(as.numeric(row.names(res)),highlight))
  non_highlight_rows=setdiff(1:nrow(res),highlight_rows)
  
  outliers=as.integer(row.names(res[res[,colq]<=FDR,]))
  
  ok_outliers=TRUE
  if (sum(res[,colq]<=FDR)==0)
    ok_outliers=FALSE;
  
  res[res[,colq]<=0.0001,colq]=0.0001
  
  # plot
  plot(log10(res[,colq]),res[,colfstat],xlim=rev(range(log10(res[,colq]))),xlab="log10(q value)",ylab=names(res[colfstat]),type="n")
  points(log10(res[non_highlight_rows,colq]),res[non_highlight_rows,colfstat],pch=19,cex=size)
  
  if (name_highlighted) {
    if (length(highlight_rows)>0) {
      text(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],row.names(res[highlight_rows,]),col="red",cex=size*1.2,font=2)
    }
  }
  else {
    points(log10(res[highlight_rows,colq]),res[highlight_rows,colfstat],col="red",pch=19,cex=size)
    # add names of loci over p and vertical line
    if (ok_outliers & add_text) {
      text(log10(res[res[,colq]<=FDR,][,colq])+pos*(round(runif(nrow(res[res[,colq]<=FDR,]),1,2))*2-3),res[res[,colq]<=FDR,][,colfstat],row.names(res[res[,colq]<=FDR,]),cex=size)
    }
  }
  lines(c(log10(FDR),log10(FDR)),c(-1,1),lwd=2)
  
  return(list("outliers"=outliers,"nb_outliers"=length(outliers)))
}
par(mfrow=c(1,1))
plot_bayescan('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/alcon/alcon_fst.txt', FDR = 0.05)#+title('Alcon')
plot_bayescan('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/gentian/gentian_fst.txt', FDR = 0.05)#+title('Gentian')

#Make own nice GGplot with clear labels
##ALCON
bayes.alc <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/alcon/alcon_fst.txt', h = T, sep = ' ')
bayes.alc <- bayes.alc[,3:7]
colnames(bayes.alc) <- c('prob','log10(PO)','qval','alpha','fst')

signSNP.bayes.alc <- bayes.alc[bayes.alc$qval <= 0.05,]
bayes.alc[bayes.alc$qval < 0.001, 3] <- 0.001 #set q-value to at least 0.001 (to show on plot)
bayes.alc$balancing <- as.factor(ifelse(test = bayes.alc$alpha<0 & bayes.alc$qval<0.05, 1, 0))


library(ggplot2)
library(ggthemes)
library(dplyr)
library(tibble)
library(ggrepel)

ggplot(data = bayes.alc, aes(x= log10(qval), y= fst, color = balancing))+
  geom_point()+
  geom_vline(xintercept = log10(0.05), color = 'tomato1')+
  scale_x_reverse()+
  xlim(0, -3) +
  theme_classic()+
  ggtitle('Alcon')+
  xlab('Log10(Q-value)')+ylab('Fst')+
  geom_label_repel(data = bayes.alc %>% filter(qval < 0.05), 
                   aes(label=rownames(signSNP.bayes.alc)), size = 3)+
  scale_color_manual(values=c('black', '#116E8A'))+ 
  theme(legend.position = "none") 


bayes.gen <- read.delim('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/gentian/gentian_fst.txt', h = T, sep = ' ')
bayes.gen <- bayes.gen[,3:7]
colnames(bayes.gen) <- c('prob','log10(PO)','qval','alpha','fst')

signSNP.bayes.gen <- bayes.gen[bayes.gen$qval <= 0.05,]
bayes.gen[bayes.gen$qval < 0.001, 3] <- 0.001 #set q-value to at least 0.001 (to show on plot)
bayes.gen$balancing <- as.factor(ifelse(test = bayes.gen$alpha<0 & bayes.gen$qval<0.05, 1, 0))

ggplot(data = bayes.gen, aes(x= log10(qval), y= fst, color = balancing))+
  geom_point()+
  geom_vline(xintercept = log10(0.05), color = 'tomato1')+
  scale_x_reverse()+
  xlim(0, -3) +
  theme_classic()+
  ggtitle('Gentian')+
  xlab('Log10(Q-value)')+ylab('Fst')+
  geom_label_repel(data = bayes.gen %>% filter(qval < 0.05), 
                   aes(label=rownames(signSNP.bayes.gen)), size = 3)+
  scale_color_manual(values=c('black', '#116E8A'))+ 
  theme(legend.position = "none") 


alcondata=read.table('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/alcon/alcon.sel',colClasses="numeric")
gentiandata=read.table('/home/jonas/Documents/Masterthesis/data/Analyses/Bayescan/BayescanOutput/gentian/gentian.sel',colClasses="numeric")
# choose the parameter you want to plot by setting for example:
parameter="Fst1"
# then this line will make the plot for:
plot(density(alcondata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
plot(density(gentiandata[[parameter]]),xlab=parameter,main=paste(parameter,"posterior distribution"))
