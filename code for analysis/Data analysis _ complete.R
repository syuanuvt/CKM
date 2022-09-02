##################################################
######### code for summary and reports ###########
##################################################
library(sparcl)
library(rARPACK)
library(combinat)
library(mclust)
library(dplyr)
library(cluster)
library(RSpectra)
library(GGally)
library(CKM)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library("forcats")

##################################################################################
##### summary of simulation study 1
##################################################################################
n.cluster <- c(3,5,30)
n.noisevar <- c(5,50,250,1000)
var <- 1
n.validvar <- 50
nmu <- c(0.6,0.7,0.8,1)
con.rep <- 1:40

condition0$ckm.ce <- rep(NA, nrow(condition0))
condition0$simple.ce <- rep(NA, nrow(condition0))
condition0$skm.ce <- rep(NA, nrow(condition0))
condition0$ckm.prop <- rep(NA, nrow(condition0))
condition0$simple.prop <- rep(NA, nrow(condition0))
condition0$skm.prop <- rep(NA, nrow(condition0))

for (i in 1:nrow(condition0)){
    load(paste0(i, ".RData"))
    condition0[i,6] <- out[[1]][[1]]
    condition0[i,7] <- out[[2]][[1]]
    if(condition0$ncluster[i] != 30){
      condition0[i,8] <- out[[3]][[1]]
    }
    valid.set <- out[[1]][[3]]
    condition0[i,9] <- sum(valid.set %in% 1:condition0$nvalidvar[i]) / condition0$nvalidvar[i]
    valid.set <- out[[2]][[3]]
    condition0[i,10] <- sum(valid.set %in% 1:condition0$nvalidvar[i]) / condition0$nvalidvar[i]
    if(condition0$ncluster[i] != 30){
      valid.set <- out[[3]][[3]]
      condition0[i,11] <- sum(valid.set %in% 1:condition0$nvalidvar[i]) / condition0$nvalidvar[i]
    }
}

#################################
## plots of classification error
#################################
condition0.ce <- condition0[,1:8]
summary.condition0.ce <- tidyr::gather(condition0.ce, methods, ce, c(ckm.ce,simple.ce,skm.ce), factor_key=TRUE)
summary.condition0.ce$methods <- as.factor(summary.condition0.ce$methods)
levels(summary.condition0.ce$methods) <- c("CKM", "SAS", "SKM")
summary.condition0.ce$nnoisevar <- as.factor(summary.condition0.ce$nnoisevar)
levels(summary.condition0.ce$nnoisevar) <- c("V=5", "V=50", "V=250", "V=1000")

a <- summary.condition0.ce %>%
  filter(ncluster == 3) %>%
  ggplot(aes(x = nnoisevar, y=ce, fill = methods)) +
  geom_boxplot() +
  scale_fill_manual(values = c("gray0", "gray50", "gray80"))+
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(0,.3) +
  ylab("Classification Error")
b <- summary.condition0.ce %>%
  filter(ncluster == 5) %>%
  ggplot(aes(x = nnoisevar, y=ce, fill = methods)) +
  scale_fill_manual(values = c("gray0", "gray50", "gray80"))+
  geom_boxplot() +
  ylim(0,.3) +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Classification Error")
c <- summary.condition0.ce %>%
  filter(ncluster == 30) %>%
  ggplot(aes(x = nnoisevar, y=ce, fill = methods)) +
  scale_fill_manual(values = c("gray0", "gray50", "gray80"))+
  ylim(0,.3) +
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ) +
  ylab("Classification Error")

figure1 <- plot_grid(a, b,c, labels = "AUTO")
#################################
## summary of other results
#################################
condition0.ce %>%
  filter(ncluster == 3) %>%
  filter(nnoisevar == 1000) %>%
  summarise(ckm = mean(ckm.ce, na.rm = TRUE), ckm.median = median(ckm.ce, na.rm = TRUE))
condition0.ce %>%
  filter(ncluster == 3) %>%
  filter(nnoisevar == 1000) %>%
  summarise(sas = mean(simple.ce, na.rm = TRUE), sas.median = median(simple.ce, na.rm = TRUE))
condition0.ce %>%
  filter(ncluster == 3) %>%
  summarise(skm = mean(skm.ce, na.rm = TRUE), skm.median = median(skm.ce, na.rm = TRUE))
####################################
condition0.prop <- condition0[,c(1:5,9:11)]
summary.condition0.prop <- tidyr::gather(condition0.prop, methods, prop, c(ckm.prop,simple.prop,skm.prop), factor_key=TRUE)
summary.condition0.prop$methods <- as.factor(summary.condition0.ce$methods)
levels(summary.condition0.prop$methods) <- c("CKM", "SAS", "SKM")
condition0.prop %>%
  filter(ncluster == 3) %>%
  filter(nnoisevar == 1000) %>%
  summarise(ckm = mean(ckm.prop, na.rm = TRUE), ckm.median = median(ckm.prop, na.rm = TRUE))
condition0.prop %>%
  filter(ncluster == 3) %>%
  filter(nnoisevar == 1000) %>%
  summarise(sas = mean(simple.prop, na.rm = TRUE), sas.median = median(simple.prop
                                                                       , na.rm = TRUE))
condition0.prop %>%
  filter(ncluster == 3) %>%
  filter(nnoisevar == 1000) %>%
  summarise(skm = mean(skm.prop, na.rm = TRUE), skm.median = median(skm.prop, na.rm = TRUE))

##################################################################################
##### summary of simulation study 2
##################################################################################
n.cluster <- c(3,5,30)
n.noisevar <- c(5,50,250,1000)
var <- 1
n.validvar <- 50
nmu <- c(0.6,0.7,0.8,1)
con.rep <- 1:40

condition1 <- expand.grid(ncluster = n.cluster, nnoisevar= n.noisevar, nmu = nmu, rep = con.rep)
condition1$ckm <- rep(NA, nrow(condition1))
condition1$simple <- rep(NA, nrow(condition1))
condition1$skm <- rep(NA, nrow(condition1))
condition1$km <- rep(NA, nrow(condition1))
condition1$ckm.length <- rep(NA, nrow(condition1))
condition1$simple.length <- rep(NA, nrow(condition1))
condition1$skm.length <- rep(NA, nrow(condition1))
condition1$km.length <- rep(NA, nrow(condition1))
condition1$ckm.time <- rep(NA, nrow(condition1))
condition1$simple.time <- rep(NA, nrow(condition1))
condition1$skm.time <- rep(NA, nrow(condition1))
condition1$km.time <- rep(NA, nrow(condition1))

for (i in 1:nrow(condition1)){
    load(paste0(i, ".RData"))
    records <- results$records
    condition1[i,5:8] <- records[1,]
    condition1[i,9:12] <- records[2,]
    condition1[i,13:16] <- records[3,]
}

condition1$nnoisevar <-as.factor(condition1$nnoisevar)
levels(condition1$nnoisevar) <- c("V=5", "V=50", "V=250",
                                  "V=1000")
condition1$nmu <- as.factor(condition1$nmu)
levels(condition1$nmu) <- c(".6", ".7", ".8",
                            "1")

#################################
## plots of classification error
#################################
summary.condition1.ce <- tidyr::gather(condition1, methods, cer, c(ckm,simple,skm,km), factor_key=TRUE)
summary.condition1.ce$methods <- as.factor(summary.condition1.ce$methods)
levels(summary.condition1.ce$methods) <- c("CKM", "SAS", "SKM", "KM")
a2 <- summary.condition1.ce %>%
  filter(ncluster == 3) %>%
  ggplot(aes(x=nnoisevar, y=cer, fill=methods)) +
  scale_fill_manual(values = c("gray0", "gray40", "gray60","gray80"))+
  ylim(0,.8) +
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        #panel.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Classification Error")
b2 <- summary.condition1.ce %>%
  filter(ncluster == 5) %>%
  ggplot(aes(x=nnoisevar, y=cer, fill=methods)) +
  scale_fill_manual(values = c("gray0", "gray40", "gray60","gray80"))+
  ylim(0,.8) +
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        #panel.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Classification Error")
c2 <- summary.condition1.ce %>%
  filter(ncluster == 30) %>%
  ggplot(aes(x=nnoisevar, y=cer, fill=methods)) +
  scale_fill_manual(values = c("gray0", "gray40", "gray60","gray80"))+
  ylim(0,.8) +
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        #panel.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("Classification Error")
plot_grid(a2, b2, c2,labels = "AUTO")

#################################
## plots of numbers of selected variables
#################################
summary.condition1.length <- tidyr::gather(condition1, methods, length, c(ckm.length, simple.length, skm.length), factor_key=TRUE)
summary.condition1.length$methods <- as.factor(summary.condition1.length$methods)
levels(summary.condition1.length$methods) <- c("CKM", "SAS", "SKM")
a2 <- summary.condition1.length %>%
  filter(ncluster == 3) %>%
  ggplot(aes(x=nnoisevar, y=length, fill=methods)) +
  scale_fill_manual(values = c("gray0", "gray50", "gray80"))+
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif")) +
  ylim(30,70)+
  ylab("Number of Selected Variables")
###########################
b2 <- summary.condition1.length %>%
  filter(ncluster == 5) %>%
  ggplot(aes(x=nnoisevar, y=length, fill=methods)) +
  scale_fill_manual(values = c("gray0", "gray50", "gray80"))+
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif")) +
  ylim(30,70)+
  ylab("Number of Selected Variables")
###########################
c2 <- summary.condition1.length %>%
  filter(ncluster == 30) %>%
  ggplot(aes(x=nnoisevar, y=length, fill=methods)) +
  scale_fill_manual(values = c("gray0", "gray50", "gray80"))+
  geom_boxplot() +
  facet_grid(paste0("\u03bc=", nmu)~nnoisevar, scale="free")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size = 16, family = "serif"),
        panel.border      = element_blank(),
        #panel.background  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylim(30,70)+
  ylab("Number of Selected Variables")
plot_grid(a2, b2,c2, labels = "AUTO")

#################################
## summary of computational time
#################################
condition1 %>%
  dplyr::filter(ncluster == 5 | ncluster == 3) %>%
  summarise(ckm = mean(ckm.time), simple = mean(simple.time), skm = mean(skm.time, na.rm = TRUE), km = mean(km.time))
condition1 %>%
  dplyr::filter(ncluster == 3) %>%
  summarise(ckm = mean(ckm.time), simple = mean(simple.time), skm = mean(skm.time, na.rm = TRUE), km = mean(km.time))
condition1 %>%
  dplyr::filter(ncluster == 5) %>%
  summarise(ckm = mean(ckm.time), simple = mean(simple.time), skm = mean(skm.time, na.rm = TRUE), km = mean(km.time))

##################################################################################
##### summary of simulation study 3
##################################################################################
set.seed(970912)
n.cluster <- c(3,5,15)
n.noisevar <- c(5, 50, 250,1000)
var <- 1
n.validvar <- 50
nmu <- c(0.6,0.7, 0.8,1)
con.rep <- 1:40
sim3 <- expand.grid(ncluster = n.cluster, nnoisevar= n.noisevar, nmu = nmu, rep = con.rep)
######################
sim3.results <- as.data.frame(matrix(nrow = nrow(sim3), ncol = 12+ncol(sim3)))
sim3.results[,1:4] <- sim3
sim3.correct <- as.data.frame(matrix(nrow = nrow(sim3), ncol = 12+ncol(sim3)))
sim3.correct[,1:4] <- sim3
for (i in 1:nrow(sim3.results)){
    load(paste0(i, ".RData"))
    records <- results[[2]]
    sim3.results[i,5:ncol(sim3.results)] <- records[1:12,1]
}
for (i in 1:nrow(sim3.correct)){
  sim3.correct[i,5:ncol(sim3.correct)] <- sim3.results[i,5:ncol(sim3.results)] == sim3.results$ncluster[i]
}
colnames(sim3.results) <- c("ncluster", "nnoisevar", "mu", "rep",
                            "all_KL", "all_DIndex","all_global", "all_first",
                            "sas_global", "sas_first","sas_KL", "sas_DIndex",
                            "ckm_global", "ckm_first","ckm_KL", "ckm_DIndex")
colnames(sim3.correct) <- c("ncluster", "nnoisevar", "mu", "rep",
                            "all_KL", "all_DIndex","all_global", "all_first",
                            "sas_global", "sas_first","sas_KL", "sas_DIndex",
                            "ckm_global", "ckm_first","ckm_KL", "ckm_DIndex")

cor.results <- sim3.correct %>%
  dplyr::group_by(ncluster, nnoisevar) %>%
  dplyr::summarise_all(mean)


##################################################################################
##### Analysis of the application
##################################################################################
autism.t <- read.csv("autism_new.csv", sep = ",")
autism <- t(autism.t)
## invalid respondents
autism <- autism[-c(1,13,14,28),]
missing <- apply(autism, 2, function(x) sum(is.na(x)))
sd.missing <- apply(autism, 2, sd)
autism.clean <- autism[,-which(sd.missing == 0)]
autism.center <- scale(autism.clean, scale = FALSE)
mean.center <- apply(autism.center, 2, function(x) sqrt(sum(x^2)))
autism.st <- sweep(autism.center, 2, mean.center, FUN = "/")
#################################
## part 1: analysis on the original data sets
n.rep <- 10
set.seed(921009)
results.ckm <- CKMSelVar(autism.st, n.cluster = 2, kmeans_starts = n.rep)
set.seed(921009)
results.km <- kmeans(autism.st, centers = 2, nstart = n.rep)
set.seed(921009)
results.sas <- hill_climb_GSS(X = autism.st, k = 2, nperms=n.rep,itermax=20,threshold=0,tolerance=1, num_starts_kmeans = n.rep)
set.seed(921009)
results.skm.pre <- try(KMeansSparseCluster.permute(autism.st, K = 2, wbounds = seq(1.001,20,len= 200), nperms=n.rep))
select.wbound <- results.skm.pre$wbounds[which(results.skm.pre$nnonzerows == n.validvar)]
mean.wbound <- mean(select.wbound)
results.skm <- KMeansSparseCluster(autism.st, K = 2, wbounds = mean.wbound, nstart = 1000)
#################################
## part 2: inspect the meaning of selected genes
#### ############################
gene <- read.csv("gene_full.csv", sep = ",")
gene.clean <- dplyr::filter(gene, CONTROL_TYPE %in% c("FALSE","pos","neg"))
gene.final <- gene.clean[-which(sd.missing == 0),]
useful.gene <- gene.final$GENE_SYMBOL[stable.set]
#################################
## part 3: analysis on selected subsets
#### ############################
## create a new data set with signaling variables
autism.new <- cbind(autism.st, c(rep(1,13), rep(2,14)))
autism.new <- as.data.frame(autism.new)
names(autism.new)[43894] <- "ID"
autism.new$ID <- as.factor(autism.new$ID)
autism.effect <- rep(NA, 43893)
for(i in 1:43893){
  a <- lm(autism.new[,i] ~ autism.new[,43894])
  autism.effect[i] <- abs(a$coefficients[2])
}
sig.genes <- sort(autism.effect, index.return = TRUE, decreasing = TRUE)$ix[1:293]
rest.genes <- setdiff(1:43893, sig.genes)
set.seed(921009)
n <- 1707
rest.selected <- sample(rest.genes, n)
new.data <- autism.st[,c(sig.genes, rest.selected)]
## analyze on this new data set
n.rep <- 10
set.seed(921009)
results.km.new <- kmeans(new.data, centers = 2, nstart = n.rep)
set.seed(921009)
results.ckm.new <- CKMSelVar(new.data, 2, kmeans_starts = n.rep)
set.seed(921009)
results.sas.new <- hill_climb_GSS(X = new.data, k = 2, nperms=n.rep,itermax=20,threshold=0,tolerance=1,num_starts_kmeans = n.rep)
set.seed(921009)
results.skm.new.pre <- try(KMeansSparseCluster.permute(new.data, K = 2, wbounds = seq(1.001,20,len= 50), nperms=n.rep))
select.wbound <- results.skm.new.pre$wbounds[which(results.skm.new.pre$nnonzerows == n.validvar)]
mean.wbound <- mean(select.wbound)
results.skm <- KMeansSparseCluster(new.data, K = 2, wbounds = mean.wbound, nstart = 1000)

