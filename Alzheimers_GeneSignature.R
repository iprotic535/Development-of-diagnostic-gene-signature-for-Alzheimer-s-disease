library(dplyr)
library(ade4)
# read file

expr = read.delim("expr_cohort_B.txt", row.names="gene")

pheno = read.delim("pheno_cohort_B.txt")

# match to patient codes

z = match(colnames(expr), pheno$sample)

pheno = pheno[z,]



# Assuming `expr` is a data frame with columns corresponding to samples and 
# `pheno` is a data frame with a column "diagnosis" that corresponds to the columns in `expr`

# Control group only
expr_ctl <- expr %>%
  select(pheno %>%
           filter(diagnosis == "CTL") %>%
           pull(sample))

# Alzheimer group only
expr_ad <- expr %>%
  select(pheno %>%
           filter(diagnosis == "AD") %>%
           pull(sample))

# Mild Cognitive Impairment (MCI) group only
expr_mci <- expr %>%
  select(pheno %>%
           filter(diagnosis == "MCI") %>%
           pull(sample))

# Combined disease group (AD and MCI)
expr_ad_mci <- expr %>%
  select(pheno %>%
           filter(diagnosis %in% c("AD", "MCI")) %>%
           pull(sample))



# making p values and fold change vectors for the t.test

p_ad = rep(NA, nrow(expr_ctl)) # control x ad

p_mci = rep(NA, nrow(expr_ctl)) # control x mci

p_ad_mci = rep(NA, nrow(expr_ctl)) # control x disease

fc_ad = rep(NA, nrow(expr_ctl))

fc_mci = rep(NA, nrow(expr_ctl))

fc_ad_mci = rep(NA, nrow(expr_ctl))



# loop for t.test

for (i in 1:nrow(expr_ctl))
  
{
  
  # test for control & ad
  
  p_ad[i] = t.test(t(expr_ctl[i,]), t(expr_ad[i,]))$p.value
  
  fc_ad[i] = mean(2^t(expr_ad[i,])) / mean(2^t(expr_ctl[i,]))
  
  # test for control & mci
  
  p_mci[i] = t.test(t(expr_ctl[i,]), t(expr_mci[i,]))$p.value
  
  fc_mci[i] = mean(2^t(expr_mci[i,])) / mean(2^t(expr_ctl[i,]))
  
  # test for control & disease
  
  p_ad_mci[i] = t.test(t(expr_ctl[i,]), t(expr_disease[i,]))$p.value
  
  fc_ad_mci[i] = mean(2^t(expr_disease[i,])) / mean(2^t(expr_ctl[i,]))                    
  
}

#t-test results
result = data.frame(rownames(expr_ctl), fc_ad, p_ad, fc_mci, p_mci, fc_ad_mci, p_ad_mci) 

# adjusting p-values to obtain genes for gene signature

#CTL x AD

sig1 = result[p.adjust(result$p_ad, method="BH")<0.05 & abs(log2(result$fc_ad))>log2(2),]

#CTL x MCI

sig2 = result[p.adjust(result$p_mci, method="BH")<0.05 & abs(log2(result$fc_mci))>log2(2),]

#CTL x Disease

sig3 = result[p.adjust(result$p_ad_mci, method="BH")<0.05 & abs(log2(result$fc_ad_mci))>log2(2),]

#making gene signatures
# CTL & AD
signature1 = data.frame(sig1[,1], sign(log2(sig1$fc_ad))) 
# CTL x MCI
signature2 = data.frame(sig2[,1], sign(log2(sig2$fc_mci))) 
# CTL x Disease
signature3 = data.frame(sig3[,1], sign(log2(sig3$fc_ad_mci))) 

colnames(signature1) = c("gene", "weight")

colnames(signature2) = c("gene", "weight")

colnames(signature3) = c("gene", "weight")

#making tables with the signatures so we can use it for validation

write.table(signature1, "signature1.txt", quote=F, row.names=F, sep="\t")

write.table(signature2, "signature2.txt", quote=F, row.names=F, sep="\t")

write.table(signature3, "signature3.txt", quote=F, row.names=F, sep="\t")




# PCA analysis

expr1 = cbind(expr_ctl, expr_ad)

expr2 = cbind(expr_ctl, expr_mci)

expr3 = cbind(expr_ctl, expr_ad, expr_mci)

z = match(signature1$gene, rownames(expr1))

expr1 = expr1[z,]

pca_out1 = dudi.pca(t(expr1), scannf = FALSE, nf = 2)

z = match(signature2$gene, rownames(expr2))

expr2 = expr2[z,]

pca_out2 = dudi.pca(t(expr2), scannf = FALSE, nf = 2)

z = match(signature3$gene, rownames(expr3))

expr3 = expr3[z,]

pca_out3 = dudi.pca(t(expr3), scannf = FALSE, nf = 2)



pdf(file="discovery_cohort.pdf", width=12, height=7, useDingbats=FALSE)

split.screen(c(1, 3))

screen(1)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out1$l1[pheno$diagnosis == "CTL",]

set2 = pca_out1$l1[pheno$diagnosis == "AD",]

plot(set1[,1], set1[,2], col='blue', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Discovery Cohort B")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'AD'), pch=20, col=c('blue', 'red'), bty='n', cex=1, pt.cex=1.5)



screen(2)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out2$l1[pheno$diagnosis == "CTL",]

set2 = pca_out2$l1[pheno$diagnosis == "MCI",]

plot(set1[,1], set1[,2], col='blue', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Discovery Cohort B")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'MCI'), pch=20, col=c('blue', 'red'), bty='n', cex=1, pt.cex=1.5)



screen(3)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out3$l1[pheno$diagnosis == "CTL",]

set2 = pca_out3$l1[pheno$diagnosis == "MCI" | pheno$diagnosis == "AD" ,]

plot(set1[,1], set1[,2], col='blue', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Discovery Cohort B")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'AD & MCI'), pch=20, col=c('blue', 'red'), bty='n', cex=1, pt.cex=1.5)

close.screen(all=TRUE)

dev.off()


## Risk Score calculation
library(matrixStats)

risk.score <- function(expr, signature)
  
{
  
  z = match(rownames(expr), signature[,1])
  
  expr = expr[!is.na(z),]
  
  signature = na.omit(signature[z,])
  
  w = signature[,2]
  
  u = rowSds(expr)
  
  m = rowMeans(expr)
  
  s = c()
  
  for (i in 1:ncol(expr))
    
  {
    
    s = append(s, ((expr[,i] - m) / u) %*% w)
    
  }
  
  s
  
}

#risk score of the gene signature #1 in the discovery cohort

score1 = risk.score(as.matrix(expr1), signature1)

score2 = risk.score(as.matrix(expr2), signature2)

score3 = risk.score(as.matrix(expr3), signature3)



library(pROC)



pdf(file="risk_scores_discovery.pdf", width=10, height=12, useDingbats=FALSE)

split.screen(c(3, 2))

screen(1)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score1[pheno$diagnosis == "CTL"], score1[pheno$diagnosis == "AD"], main = "Control & AD", ylab = "score", names = c("Control", "AD"), col = c("forestgreen","magenta"))

screen(2)

cl1 = c(rep("CTL", ncol(expr_ctl)), rep("AD", ncol(expr_ad)))

roc1 = roc(cl1,score1)

plot(roc1)

screen(3)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score2[pheno$diagnosis == "CTL"], score2[pheno$diagnosis == "MCI"], main = "Control & MCI", ylab = "score", names = c("Control", "MCI"), col = c("forestgreen","magenta"))

screen(4)

expr_combined = cbind(expr_ctl, expr_mci)

cl2 = c(rep("CTL", ncol(expr_ctl)), rep("MCI", ncol(expr_mci)))

roc2 = roc(cl2,score2)

plot(roc2)

screen(5)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score3[pheno$diagnosis == "CTL"], score3[pheno$diagnosis == "AD" | pheno$diagnosis == "MCI"], main = "Control & MCI-AD", ylab = "score", names = c("Control", "MCI-AD"), col = c("forestgreen","magenta"))

screen(6)

expr_combined = cbind(expr_mci, expr_ad)

cl3 = c(rep("CTL", ncol(expr_ctl)), rep("AD", ncol(expr_ad)), rep("MCI", ncol(expr_mci)))

roc3 = roc(cl3,score3)

plot(roc3)

close.screen(all=TRUE)

dev.off()

# Reading cohort A for validation

expr_val = read.delim("expr_cohort_A.txt", row.names="gene")

pheno_val = read.delim("pheno_cohort_A.txt")



#gathering experimental groups
#control 
expr_ctl = expr_val[,pheno_val$diagnosis == "CTL"] 
# Alzheimer
expr_ad = expr_val[,pheno_val$diagnosis == "AD"]  
# MCI
expr_mci = expr_val[,pheno_val$diagnosis == "MCI"] 

expr_disease = expr_val[,pheno_val$diagnosis == "AD" | pheno_val$diagnosis == "MCI"] # combined disease group (AD and MCI)



# loading the gene signatures

signature11 = read.delim("signature1.txt") #CTL x AD signature

signature21 = read.delim("signature2.txt") #CTL x MCI signature

signature31 = read.delim("signature3.txt") #CTL x Disease signature



expr1_val = cbind(expr_ctl, expr_ad)

expr2_val = cbind(expr_ctl, expr_mci)

expr3_val = cbind(expr_ctl, expr_ad, expr_mci)



#making individual matches for each signature

z = match(signature1$gene, rownames(expr1_val))

expr1_val = expr1_val[z,]

pca_out1 = dudi.pca(t(expr1_val), scannf = FALSE, nf = 2)

z = match(signature2$gene, rownames(expr2_val))

expr2_val = expr2_val[z,]

pca_out2 = dudi.pca(t(expr2_val), scannf = FALSE, nf = 2)

z = match(signature3$gene, rownames(expr3_val))

expr3_val = expr3_val[z,]

pca_out3 = dudi.pca(t(expr3_val), scannf = FALSE, nf = 2)



#validation plot

pdf(file=" validation_a.pdf", width=16, height=7, useDingbats=FALSE)

split.screen(c(1, 3))

screen(1)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out1$l1[pheno_val$diagnosis == "CTL",]

set2 = pca_out1$l1[pheno_val$diagnosis == "AD",]

plot(set1[,1], set1[,2], col='forestgreen', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Validation of Cohort B with Cohort A")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'AD'), pch=20, col=c('forestgreen', 'red'), bty='n', cex=1, pt.cex=1.5)



screen(2)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out2$l1[pheno_val$diagnosis == "CTL",]

set2 = pca_out2$l1[pheno_val$diagnosis == "MCI",]

plot(set1[,1], set1[,2], col='forestgreen', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Validation of Cohort B with Cohort A")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'MCI'), pch=20, col=c('forestgreen', 'red'), bty='n', cex=1, pt.cex=1.5)



screen(3)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out3$l1[pheno_val$diagnosis == "CTL",]

set2 = pca_out3$l1[pheno_val$diagnosis == "MCI" | pheno_val$diagnosis == "AD" ,]

plot(set1[,1], set1[,2], col='forestgreen', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Validation of Cohort B with Cohort A")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'MCI-AD'), pch=20, col=c('forestgreen', 'red'), bty='n', cex=1, pt.cex=1.5)



close.screen(all=TRUE)

dev.off()



library(matrixStats)

risk.score <- function(expr, signature)
  
{
  
  z = match(rownames(expr), signature[,1])
  
  expr = expr[!is.na(z),]
  
  signature = na.omit(signature[z,])
  
  w = signature[,2]
  
  u = rowSds(expr)
  
  m = rowMeans(expr)
  
  s = c()
  
  for (i in 1:ncol(expr))
    
  {
    
    s = append(s, ((expr[,i] - m) / u) %*% w)
    
  }
  
  s
  
}



#calculating risk score

score11 = risk.score(as.matrix(expr1_val), signature11)

score21 = risk.score(as.matrix(expr2_val), signature21)

score31 = risk.score(as.matrix(expr3_val), signature31)



pdf(file="risk_scores_validation_a.pdf", width=10, height=12, useDingbats=FALSE)

split.screen(c(3, 2))

screen(1)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score1[pheno_val$diagnosis == "CTL"], score1[pheno_val$diagnosis == "AD"], main = "Control & AD", ylab = "score", names = c("Control", "AD"), col = c("forestgreen","magenta"))

screen(2)

cl1 = c(rep("CTL", ncol(expr_ctl)), rep("AD", ncol(expr_ad)))

roc11 = roc(cl1,score11)

plot(roc11)

screen(3)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score2[pheno_val$diagnosis == "CTL"], score1[pheno_val$diagnosis == "MCI"], main = "Control & MCI", ylab = "score", names = c("Control", "MCI"), col = c("forestgreen","magenta"))

screen(4)

cl2 = c(rep("CTL", ncol(expr_ctl)), rep("MCI", ncol(expr_mci)))

roc21 = roc(cl2,score21)

plot(roc21)

screen(5)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score3[pheno_val$diagnosis == "CTL"], score1[pheno_val$diagnosis == "AD" | pheno_val$diagnosis == "MCI"], main = "Control & MCI-AD", ylab = "score", names = c("Control", "MCI-AD"), col = c("forestgreen","magenta"))

screen(6)

cl3 = c(rep("CTL", ncol(expr_ctl)), rep("AD", ncol(expr_ad)), rep("MCI", ncol(expr_mci)))

roc31 = roc(cl3,score31)

plot(roc31)

close.screen(all=TRUE)

dev.off()

# Reading cohort C for validation

expr_val = read.delim("expr_cohort_c.txt", row.names="gene")

pheno_val = read.delim("pheno_cohort_c.txt")



expr_ctl = expr_val[,pheno_val$diagnosis == "CTL"] #control group only

expr_ad = expr_val[,pheno_val$diagnosis == "AD"]  # Alzheimer group only

expr_mci = expr_val[,pheno_val$diagnosis == "MCI"] # mild cognitive impairment group only

expr_disease = expr_val[,pheno_val$diagnosis == "AD" | pheno_val$diagnosis == "MCI"] # combined disease group (AD and MCI)



# loading the gene signatures

signature13 = read.delim("signature1.txt")

signature23 = read.delim("signature2.txt")

signature33 = read.delim("signature3.txt")



# making individual matches for each signature

expr1_val = cbind(expr_ctl, expr_ad)

expr2_val = cbind(expr_ctl, expr_mci)

expr3_val = cbind(expr_ctl, expr_ad, expr_mci)



z = match(signature1$gene, rownames(expr1_val))

expr1_val = expr1_val[z,]

pca_out1 = dudi.pca(t(expr1_val), scannf = FALSE, nf = 2)

z = match(signature2$gene, rownames(expr2_val))

expr2_val = expr2_val[z,]

pca_out2 = dudi.pca(t(expr2_val), scannf = FALSE, nf = 2)

z = match(signature3$gene, rownames(expr3_val))

expr3_val = expr3_val[z,]

pca_out3 = dudi.pca(t(expr3_val), scannf = FALSE, nf = 2)



pdf(file=" validation_c.pdf", width=16, height=7, useDingbats=FALSE)

split.screen(c(1, 3))

screen(1)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out1$l1[pheno_val$diagnosis == "CTL",]

set2 = pca_out1$l1[pheno_val$diagnosis == "AD",]

plot(set1[,1], set1[,2], col='purple', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Validation of Cohort B with Cohort C")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'AD'), pch=20, col=c('purple', 'red'), bty='n', cex=1, pt.cex=1.5)



screen(2)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out2$l1[pheno_val$diagnosis == "CTL",]

set2 = pca_out2$l1[pheno_val$diagnosis == "MCI",]

plot(set1[,1], set1[,2], col='purple', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Validation of Cohort B with Cohort C")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'MCI'), pch=20, col=c('purple', 'red'), bty='n', cex=1, pt.cex=1.5)



screen(3)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

set1 = pca_out3$l1[pheno_val$diagnosis == "CTL",]

set2 = pca_out3$l1[pheno_val$diagnosis == "MCI" | pheno_val$diagnosis == "AD" ,]

plot(set1[,1], set1[,2], col='purple', xlim=c(-3, 3), ylim=c(-3, 3), xlab='PC1', ylab='PC2', pch=20, cex=1.5, cex.lab=1, axes=FALSE, main="Validation of Cohort B with Cohort C")

points(set2[,1], set2[,2], col='red', pch=20, cex=1.2)

axis(1)

axis(2)

legend('topleft', c('Control', 'MCI-AD'), pch=20, col=c('purple', 'red'), bty='n', cex=1, pt.cex=1.5)



close.screen(all=TRUE)

dev.off()



score13 = risk.score(as.matrix(expr1_val), signature1)

score23 = risk.score(as.matrix(expr2_val), signature2)

score33 = risk.score(as.matrix(expr3_val), signature3)





pdf(file="risk_scores_validation_c.pdf", width=10, height=12, useDingbats=FALSE)

split.screen(c(3, 2))

screen(1)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score13[pheno_val$diagnosis == "CTL"], score13[pheno_val$diagnosis == "AD"], main = "Control & AD", ylab = "score", names = c("Control", "AD"), col = c("forestgreen", "magenta"))

screen(2)

cl1 = c(rep("CTL", ncol(expr_ctl)), rep("AD", ncol(expr_ad)))

roc13 = roc(cl1,score13)

plot(roc13)

screen(3)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score23[pheno_val$diagnosis == "CTL"], score23[pheno_val$diagnosis == "MCI"], main = "Control & MCI", ylab = "score", names = c("Control", "MCI"), col = c("forestgreen", "magenta"))

screen(4)

cl2 = c(rep("CTL", ncol(expr_ctl)), rep("MCI", ncol(expr_mci)))

roc23 = roc(cl2,score23)

plot(roc23)

screen(5)

par(mai=c(0.8, 0.8, 0.3, 0.3), mgp=c(2, 0.5, 0), tck=-0.03)

boxplot(score33[pheno_val$diagnosis == "CTL"], score33[pheno_val$diagnosis == "AD" | pheno_val$diagnosis == "MCI"], main = "Control & MCI-AD", ylab = "score", names = c("Control", "MCI-AD"), col = c("forestgreen", "magenta"))

screen(6)

cl3 = c(rep("CTL", ncol(expr_ctl)), rep("AD", ncol(expr_ad)), rep("MCI", ncol(expr_mci)))

roc33 = roc(cl3,score33)

plot(roc33)

close.screen(all=TRUE)

dev.off()
