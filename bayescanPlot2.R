#Make own nice GGplot with clear labels
##read bayescan file
alc.all <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Bayescan/BayescanOutput/alcon/alcon_fst.txt', h = T, sep = ' ')
alc.conn <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Bayescan/BayescanOutput/alcon_conn/alcon_conn_fst.txt', h = T, sep = ' ')
alc.isol <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Bayescan/BayescanOutput/alcon_isol/alcon_isol_fst.txt', h = T, sep = ' ')


#function to clean up bayescan file
prep.fst.txt <- function(bayescanfile){
  bayescanfile <- bayescanfile[,3:7]
  colnames(bayescanfile) <- c('prob','log10(PO)','qval','alpha','fst')
  bayescanfile[bayescanfile$qval < 0.001, 3] <- 0.001 #set q-value to at least 0.001 (to show on plot)
  bayescanfile$balancing <- as.factor(ifelse(test = bayescanfile$alpha<0 & bayescanfile$qval<0.05, 1, 0))
  bayescanfile$snpid <- rownames(bayescanfile)
  return(bayescanfile)
}

alc.all <- prep.fst.txt(alc.all)
alc.conn <- prep.fst.txt(alc.conn)
alc.isol <- prep.fst.txt(alc.isol)

library(ggplot2)
library(dplyr)
library(ggrepel)

#all populations
ggplot(data = alc.all, aes(x= log10(qval), y= fst, color = balancing))+
  geom_point()+
  geom_vline(xintercept = log10(0.05), color = 'tomato1')+
  scale_x_reverse()+
  xlim(0, -3) +
  theme_classic()+
  ggtitle('Alcon - all populations')+
  xlab('Log10(Q-value)')+ylab('Fst')+
  geom_label_repel(data = alc.all %>% filter(qval < 0.05), 
                   aes(label=snpid), size = 3)+
  scale_color_manual(values=c('black', '#116E8A'))+ 
  theme(legend.position = "none") 

#connected populations only
ggplot(data = alc.conn, aes(x= log10(qval), y= fst, color = balancing))+
  geom_point()+
  geom_vline(xintercept = log10(0.05), color = 'tomato1')+
  scale_x_reverse()+
  xlim(0, -3) +
  theme_classic()+
  ggtitle('Alcon - connected populations')+
  xlab('Log10(Q-value)')+ylab('Fst')+
  geom_label_repel(data = alc.conn %>% filter(qval < 0.05), 
                   aes(label=snpid), size = 3)+
  scale_color_manual(values=c('black', '#116E8A'))+ 
  theme(legend.position = "none") 

#isolated populations only (SW1, 5 and 6 and all SE pops)
ggplot(data = alc.isol, aes(x= log10(qval), y= fst, color = balancing))+
  geom_point()+
  geom_vline(xintercept = log10(0.05), color = 'tomato1')+
  scale_x_reverse()+
  xlim(0, -3) +
  theme_classic()+
  ggtitle('Alcon - isolated populations')+
  xlab('Log10(Q-value)')+ylab('Fst')+
  geom_label_repel(data = alc.isol %>% filter(qval < 0.05), 
                   aes(label=snpid), size = 3)+
  scale_color_manual(values=c('black', '#116E8A'))+ 
  theme(legend.position = "none") 

#NOO isolated pops here based on dbMEM -> all pops are connected #Gentian
bayes.gen <- read.delim('/home/jonas/Documents/Masterthesis/thesis/data/Analyses/Bayescan/BayescanOutput/gentian/gentian_fst.txt', h = T, sep = ' ')
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

