# create barplot 
library(ggpubr)
par(mfrow=c(1,1))
a1=read.csv("paraloags_measure.csv")
Phenotype=a1$Phenotype
genes=a1$Number_of_Orthogroups
data=data.frame(Phenotype,genes)

genes1=a1$Total_genes


p<-ggboxplot(data, x = "Phenotype", y ="genes", ylim=c(1600,7000),
               color = "Phenotype", palette =c("red","blue"),
               xlab="Phenotype",ylab ="Number of Orthogroups", bxp.errorbar=TRUE,
          bxp.errorbar.width=0.1,width = 0.5)

p + stat_compare_means(method = "t.test", hide.ns=TRUE)
p=p+stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", label.x.npc = "center")


######
#box plot genes  number and Phenotype

p1<-ggboxplot(data, x = "Phenotype", y ="genes1",
             color = "Phenotype", palette =c("red","blue"),
             xlab="Phenotype",ylab ="Number of coding genes", bxp.errorbar=TRUE,
             bxp.errorbar.width=0.1,width = 0.5)

p + stat_compare_means(method = "t.test", hide.ns=TRUE)
p1=p1+stat_compare_means(aes(label = ..p.signif..),
                     method = "t.test", label.x.npc = "center")


ggarrange(p, p1+ rremove("x.text"),
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
#Regression

plot(Genome_size~Number_of_Orthogroups,a1, cex=1.5,,xlab="Genome size", ylab="Number of Orthogroups")
abline(lm(Number_of_Orthogroups~Genome_size,a1))
