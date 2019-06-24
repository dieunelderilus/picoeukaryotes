#author : Dieunel Derilus
This scripts could takes different datasets and plot multiple regressions in one  single plots
data=read.csv("pC_significant_R1.csv",row.names = 1)
data1=read.csv("pc_not_significant2.csv",row.names = 1)
par(mfcol=c(2,2))
#plot1
#Cysteine.and.methionine.metabolism
plot(Cysteine.and.methionine.metabolism~Proteome_size,pch=16,cex=1.5,data,xlab="Proteome size (codons)", main="A",ylim=c(0,40), ylab="Genes in pathway")
abline(lm(Cysteine.and.methionine.metabolism~Proteome_size,data),lty=1, lwd=3,col="black")
#Glycine..serine.and.threonine.metabolism
points(Glycine..serine.and.threonine.metabolism~Proteome_size,data,col="red",pch=17,cex=1.5,ylim=c(0,14))
abline(lm(Glycine..serine.and.threonine.metabolism~Proteome_size,data),lty=1, lwd=3,col="red")
#Arginine.biosynthesis
points(Arginine.biosynthesis~Proteome_size,data,col="blue",pch=18,cex=1.5)
abline(lm(Arginine.biosynthesis~Proteome_size,data),lty=1, lwd=3,col="blue")


# add legend
legend(x="topleft", legend=c("Cysteine and methionine metabolism(R = 0.61, p=0.01)","Glycine,serine and threonine metabolism(R=0.75, p=0.0007)",
                             "Arginine.biosynthesis(R=0.61,p=0.01)"), bty="n",lty=c(1,1,1),pch=c(16,17,18),text.font=2,
       col=c("black","red","blue"),cex=0.7)


########################################plot2
#Phenylalanine..tyrosine.and.tryptophan.biosynthesis..PATH.ko00400
plot(Phenylalanine..tyrosine.and.tryptophan.biosynthesis..PATH.ko00400.~Proteome_size,pch=18,cex=1.5,ylim=c(0,25),data1,xlab="Proteome size (codons)", main="C",ylab="Genes in pathway")
abline(lm(Phenylalanine..tyrosine.and.tryptophan.biosynthesis..PATH.ko00400.~Proteome_size,data1),lty=1, lwd=3,col="black")
#Alanine..aspartate.and.glutamate.metabolism 
points(Alanine..aspartate.and.glutamate.metabolism~Proteome_size,data,col="red",pch=17,cex=1.5)
abline(lm(Alanine..aspartate.and.glutamate.metabolism~Proteome_size,data),lty=1, lwd=3,col="red")
#Lysine biosynthesis  
points(Lysine.biosynthesis~Proteome_size,data,col="darkgreen",pch=16,cex=1.5,ylim=c(0,14))
abline(lm(Lysine.biosynthesis ~Proteome_size,data),lty=1, lwd=3,col="darkgreen")

# add legend
legend(x="topleft", legend=c("Phenylalanine,tyrosine & tryptophan biosynthesis(R=0.53,p=0.03)","Alanine,aspartate and glutamate metabolism(R=0.67, p=0.004)",
                             "Lysine.biosynthesis(R=0.71,p=0.002)"), bty="n",lty=c(1,1,1),pch=c(18,17,16),text.font=2,
       col=c("black","red","darkgreen"),cex=0.7)


###############################################################################plot3
#Arginine biosynthesis 
plot(Arginine.biosynthesis~Proteome_size,pch=16,cex=1.5,data,xlab="Proteome size (codons)",main="B",ylim=c(0,20), ylab="Genes in pathway")
abline(lm(Arginine.biosynthesis~Proteome_size,data),lty=1, lwd=3,col="black")
#Histidine.metabolism..PATH.ko00340.  
points(Histidine.metabolism..PATH.ko00340.~Proteome_size,data1,col="red",pch=18,cex=1.5,ylim=c(0,14))
abline(lm(Histidine.metabolism..PATH.ko00340.~Proteome_size,data1),lty=1, lwd=3,col="red")
#Tyrosine.metabolism..PATH.ko00350. 
points(Tyrosine.metabolism..PATH.ko00350.~Proteome_size,data1,col="blue",pch=16,cex=1.5,ylim=c(0,14))
abline(lm(Tyrosine.metabolism..PATH.ko00350.~Proteome_size,data1),lty=1, lwd=3,col="blue")
#Lysine degradation 
points(Lysine.degradation~Proteome_size,data,col="darkgreen",pch=18,cex=1.5)
abline(lm(Lysine.degradation~Proteome_size,data),lty=1, lwd=3,col="darkgreen")
# add legend
legend(x="topleft", legend=c("Arginine.biosynthesis(R = 0.61, p=0.01)","Histidine metabolism(R=0.56, p=0.02)",
                             "Tyrosine.metabolism(R=0.53,p=0.03)","Lysine degradation(R=0.70,p=0.002)"), bty="n",lty=c(1,1,1,1),pch=c(16,18,16,18),text.font=2,
       col=c("black","red","blue","darkgreen"),cex=0.7)




#####################################################################################plot4
#Transfer.RNA.biogenesis 
plot(Tryptophan.metabolism..PATH.ko00380.~Proteome_size,pch=16,cex=1.5,data1, main="D",xlab="Proteome size (codons)",ylim=c(0,39), ylab="Genes in pathway")
abline(lm(Tryptophan.metabolism..PATH.ko00380.~Proteome_size,data1),lty=1, lwd=3,col="black")
#Messenger RNA biogenesis  
points(Valine..leucine.and.isoleucine.biosynthesis..PATH.ko00290.~Proteome_size,data1,col="red",pch=18,cex=1.5)
abline(lm(Valine..leucine.and.isoleucine.biosynthesis..PATH.ko00290.~Proteome_size,data1),lty=1, lwd=3,col="red")

#Arginine.and.proline.metabolism
points(Arginine.and.proline.metabolism..PATH.ko00330.~Proteome_size,data1,col="darkgreen",pch=17,cex=1.5,ylim=c(0,14))
abline(lm(Arginine.and.proline.metabolism..PATH.ko00330.~Proteome_size,data1),lty=1, lwd=3,col="darkgreen")

# add legend
legend(x="topleft", legend=c("Tryptophan.metabolism(R = 0.16, p=0.53)","Valine,leucine and isoleucine biosynthesis(R=0.19, p=48)","Arginine and proline metabolism, R=0.45, p=0.08)"),
                             bty="n",lty=c(1,1,1),pch=c(16,18,17),text.font=2,
       col=c("black","red","darkgreen"),cex=0.7)

