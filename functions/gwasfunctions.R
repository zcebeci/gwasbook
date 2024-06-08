

# Kod 4.1: Listeler listesini tek listeye dönüştürme
unlist <- function (x){
   xl <- x[[1L]]
   if (is.factor(xl)){
      structure(unlist(x), class="factor", levels=levels(xl))
   } else {
      do.call(c, x)
   }
}




# Kod 6.1: Genomik şişme faktörü hesaplama
 lambdagc <- function(x){
  x <- x[!is.na(x)]
  return(median(qchisq(x, df=1, lower.tail=FALSE)) / qchisq(0.5, 1))
 }





# Kod 6.2: Normallik kontrolü için tanı grafikleri
 diagPlots <- function(x, trname="", color="orange"){
   opar <- par(mfrow=c(2,2))
   hist(x, col=color, xlab=trname, main=paste(trname, "histogramı"))
   boxplot(x, pch=19, cex=0.8, col=color, horizontal=TRUE,
     xlab=trname, main=paste(trname, "kutu grafiği"))
   qqnorm(x, pch=1, frame=F, col=color, main=paste(trname,"QQ grafiği"))
   qqline(x, col="steelblue", lwd=2)
   plot(ecdf(x), col=color, xlab=trname, main=paste(trname, "ECDF"))
   par(opar)
 }





# Kod 6.3: Ters normal dönüşüm fonksiyonu
 invnormtrans <- function(x){
   return(qnorm((rank(x, na.last="keep")-0.5)/sum(!is.na(x))))
 }





# Kod 6.4: Tukey yöntemi ile aykırı değer saptama
#
 findOutliers <- function(x, color="orange", plot=FALSE){
   bp <- boxplot(x, horizontal=TRUE, plot=plot, col=color, pch=19,
      cex=1, main="Kutu-bıyık grafiği")
   if(length(bp$out)=1){
     outliers <- which(x %in% bp$out)
   }else{
     outliers <- NULL
   }
   return(outliers)
 }





# Kod 6.5: Aykırı değer temizleme fonksiyonu
#
 removeOutliers <- function(x){
   repeat{
     bp <- boxplot(x, plot=FALSE)
     if(length(bp$out)=1){
       outliers <- which(x %in% bp$out)
       x <- x[-outliers]
     }else{
       break
     }
   }
   return(x)
 }



# Kod 6.6: Tanımlayıcı istatistiklerle genotip tamamlama fonksiyonu
#
 genImp <- function(x, method="mean", custval=1){
   .impute <- function(x){
      if(method=="custom") {
         x <- replace(x, is.na(x), custval)
      }else if(method=="mean"){
         x <- replace(x, is.na(x), mean(x, na.rm = TRUE))
      }else if(method=="median"){
         x <- replace(x, is.na(x), median(x, na.rm = TRUE))
      }else if(method=="mode"){
         x <- replace(x, is.na(x), which.max(table(unique(x)))-1)
      }else if(method=="min"){
         x <- replace(x, is.na(x), which.min(table(unique(x)))-1)
      }else if(method=="sample"){
         x <- replace(x, is.na(x), sample(x,1))
      }else if(method=="rand"){
         x <- replace(x, is.na(x), sample(c(0,1,2), 1))
      }
      return(x)
   }
   return(apply(x, 2, .impute))
 }



# Kod 6.7: Çok sınıflı problemler için başarım ölçütleri
 perfMetrics <- function(actualClasses, predictedClasses){
   cl <- factor(union(unique(actualClasses), unique(predictedClasses)))
   actualClasses <- factor(actualClasses, levels=levels(cl))
   predictedClasses <- factor(predictedClasses, levels=levels(cl))
   classesTable <- table(actualClasses, predictedClasses, 
      dnn=c("Gerçek","Tahmin"))
   confMat <- as.matrix(classesTable)
   N <- sum(confMat) 
   nClasses <- nrow(confMat) 
   cmDiags <- diag(confMat)
   rowSums <- apply(confMat, 1, sum)
   colSums <- apply(confMat, 2, sum)
   actualClassFreqs <- rowSums / N 
   predClassFreqs <- colSums / N 
   classesName <- names(cmDiags)
   if(is.null(classesName)) 
       classesName <- paste("CL",(1:nClasses), sep="")

   # Sınıflara göre ölçütleri hesapla
   pcAccuracy <- 1/N*sum(cmDiags)
   pcPrecision <- cmDiags / colSums
   pcRecall <- cmDiags / rowSums
   pcF1 <- 2 * pcPrecision * pcRecall / (pcPrecision  pcRecall)
   
   # Makro başarım ölçütleri
   macRecall <- mean(pcRecall)
   macPrecision <- mean(pcPrecision)
   macF1 <- mean(pcF1)
   
   # Sınıflara göre matrisler
   sepClassMats <- lapply(1:nClasses, function(x){
        mat <- c(confMat[x,x],
        rowSums[x] - confMat[x,x],
        colSums[x] - confMat[x,x],
        N-rowSums[x] - colSums[x]  confMat[x,x])
        return(matrix(mat, nrow=2, byrow=TRUE))}
   )
   # Toplamlar matrisi
   combMat <- matrix(0, nrow=2, ncol=2)
   for(i in 1:nClasses){
     combMat <- combMat  sepClassMats[[i]]
   }
   
   # Ortalama doğruluk ve mikro başarım ölçütleri
   avgAccuracy <- sum(diag(combMat))/sum(combMat)
   micRecall <- (diag(combMat) / apply(combMat,1, sum))[1]
   micPrecision <- micF1 <- micRecall
   
   # Çoğunluk sınıfı ölçütleri
   majorClassIdx <- which(rowSums==max(rowSums))[1]
   majorClassAccuracy <- as.numeric(actualClassFreqs[majorClassIdx]) 
   majorClassRecall <- 0*actualClassFreqs
   majorClassRecall[majorClassIdx] <- 1
   majorClassPrecision <- 0*actualClassFreqs
   majorClassPrecision[majorClassIdx] <- actualClassFreqs[majorClassIdx]
   majorClassF1 <- 0*actualClassFreqs 
   majorClassF1[majorClassIdx] = 2*majorClassPrecision[majorClassIdx] /
     (majorClassPrecision[majorClassIdx]1)
   
   # Kappa istatistiği 
   expectedAccuracy <- sum(actualClassFreqs * predClassFreqs)
   kappaStat <- (pcAccuracy - expectedAccuracy) / (1 - expectedAccuracy)
   
   classMetrics <- data.frame(
     Accuracy = pcAccuracy,
     Precision = pcPrecision,
     Recall = pcRecall,
     F1 = pcF1)
   macroAverage <- data.frame(
     Precision = macPrecision,
     Recall = macRecall,
     F1 = macF1)
   microAverage <- data.frame(
     Accuracy = avgAccuracy,
     Precision = micPrecision,
     Recall = micRecall,
     F1 = micF1)
   majorClassMetrics <- data.frame(
     Class=majorClassIdx,
     Accuracy = majorClassAccuracy,
     Precision = majorClassPrecision[majorClassIdx],
     Recall = majorClassRecall[majorClassIdx],
     F1 = majorClassF1[majorClassIdx])

   results <- list(cm=confMat, classmetrics=classMetrics,
     macrometrics=macroAverage, micrometrics=microAverage,
     majoritymetrics=majorClassMetrics, Kappa=kappaStat)
   return(results)
 }




# Kod 6.8: t-SNE tereddüt (perplexity) puanları denemesi
 tsnePerplexityTuning <- function(x, d=2, pp=10, 
   itermax=500, eta=200, labs=NULL, colvec=NULL, seed=NULL){
   if(!is.null(seed)) set.seed(seed) 
   resTSNE <- Rtsne(x, dims=d, perplexity=pp, 
      max_iter=itermax, eta=eta,
      verbose=FALSE)
   expName <- paste("perplexity=",pp, ", max_iter = ", itermax, 
     ", eta = ", eta)
   opar <- par(ask=TRUE)
   if(!is.null(labs)){
     plot(resTSNE$Y, t='n', "cex.main"=1, "cex.lab"=1.5,
       xlab="Boyut 1", ylab="Boyut 2", main=expName)
      text(resTSNE$Y, labels=labs, col=colvec[as.factor(labs)])
   }else{
     plot(resTSNE$Y, pch=19, "cex.main"=1, "cex.lab"=1.5,
       xlab="Boyut 1", ylab="Boyut 2", main=expName, col="skyblue")
   }
   par(opar)
 }



# Kod 6.9: Genotiplerin normalleştirilmesi
 normalizegeno <- function(geno){
  m <- nrow(geno)
  p <- colMeans(geno)
  normGeno <- geno - rep(p, each=m)
  normGeno <- normGeno/sqrt(2*rep(p/2*(1-p/2), each=m))
  return(normGeno)
 }



# Kod 6.10: Forni yöntemleriyle G matrisi
 forniG <- function(geno, option=1){
   m <- ncol(geno)
   if(option==1){ 
     p <- 0.5
     W <- scale(x=geno, center=rep(2*p, m), scale=FALSE)
     Dnom <- 2*p*(1-p)*m 
     G <- tcrossprod(W)/Dnom
   }else if(option==2){ 
     p <- apply(X=geno, 2, FUN=function(x){p=mean(x)/2})
     p <- ifelse(p0.5, (1-p), p)
     p <- mean(p)
     W <- scale(x=geno, center=rep(2*p, m), scale=FALSE)
     Dnom <- 2*p*(1-p)*m 
     G <- tcrossprod(W)/Dnom
   }else if(option==3){
     W <- scale(x=geno, center=TRUE, scale=FALSE)
     WW <- tcrossprod(W)
     Dnom <- sum(diag(WW))/m
     G <- G/Dnom
   }
   return(G)
 }




# Kod 8.1: GWAS sonuç simülasyonu
 simGwasResult <- function(nchr=23, nsnp=1000, chrX=NULL, ntraits=1,
   posrange=c(1e6, 2e6, 3e6, 4e6, 5e6), maxpos=NULL, sortpos=TRUE, 
   nsigp=1, minsplev=1e-16, maxsplev=1e-03, seed=NULL){
   if(!missing(seed)) set.seed(seed)
   CHR <- c(); POS <- c(); SNP <- c()
   for(i in 1:nchr) {
     CHR <- c(CHR, rep(i, nsnp))
     spos <- ifelse(is.null(maxpos), sample(posrange, 1), maxpos)
     chrPOS <- sample(1:spos, nsnp, replace=FALSE)
     if(sortpos) chrPOS <- sort(chrPOS)
     POS <- c(POS, chrPOS)
     SNP <- c(SNP, paste0(i, "_", 1:nsnp))
   }
   if(!is.null(chrX)) {
     chr <- ifelse(chrX=="n", nchr1, "X")
     CHR <- c(CHR, rep(chr, nsnp))
     spos <- ifelse(is.null(maxpos), sample(posrange, 1), maxpos)
     chrPOS <- sample(1:spos, nsnp, replace=FALSE)
     chrPOS <- if(sortpos) chrPOS <- sort(chrPOS)
     POS <- c(POS, chrPOS)
     SNP <- c(SNP, paste0(chr, "_", 1:nsnp))
   }
   n <- length(POS)
   for(i in 1:ntraits){
     pval <- runif(n)
     sigpidx <- sample(1:n, nsigp, replace=FALSE)
     pval[sigpidx] <- runif(nsigp, minsplev, maxsplev)
     if(i==1)
       P <- as.data.frame(pval)
     else
       P <- cbind(P, pval)
   }
   colnames(P) <- paste0("P", 1:ntraits)
   gwasresdf <- data.frame(SNP, CHR, POS, P)
   return(gwasresdf)
 }




# Kod 8.2: QQ grafiği fonksiyonu
 gwas_qq <-function(pvals){
    n <- length(pvals)
    exppvals <- 1:n
    sobspvals <- sort(pvals, decreasing=FALSE)
    neglogobspvals <- -log(sobspvals, base=10)
    neglogexppvals <- -log(exppvals / (n1), base=10)
    maxobspval <- max(neglogobspvals)
    maxexppval <- max(neglogexppvals)
    pvdf <- data.frame(Gozlenen=neglogobspvals, Beklenen=neglogexppvals)
    plot(c(0, maxexppval), c(0, maxobspval), 
      col="red", lwd=2, type="l", 
      xlab=expression(Beklenen:  -log[10] (P)),
      ylab=expression(Gözlenen:  -log[10] (P)),
      xlim=c(0, maxexppval), ylim=c(0, maxobspval), 
      xaxs="i", yaxs="i", bty="l", las=1)
    points(neglogexppvals, neglogobspvals, pch=19, cex=0.8, bg="gray") 
    return(pvdf)
 }



# Kod 8.3: Manhattan grafiği
 gwas_manhattan <- function(gwasresult, palette, bcline=TRUE){
   completes <- !is.na(gwasresult$P)
   chr <- gwasresult$CHR[completes]
   p <- gwasresult$P[completes]
   neglogp <- -log(p, base=10)
   ylim <- c(0, max(neglogp)*1.1)
   if(missing(palette))
     palette = c("dodgerblue","orange")
   plot(neglogp, pch=20, 
      col=ifelse(as.integer(chr)%%2==0, palette[1], palette[2]),
      ylim=ylim, axes=FALSE,
      xlab="Kromozom", ylab=expression(-log[10] (P))) 
   axis(2)
   if(bcline){
     bonfcor <- -log(0.05/length(p), base=10)
     abline(h=bonfcor, lty=2, col=1, lwd=2)
     legend("topright","Bonferroni eşiği", lty=2, col=1, lwd=2)
   }
 }







# Kod 9.1: Genotip tamamlama
 imputegeno <- function(geno){
   cat("Imputation started..","\n")
   impidx <- c()
   for(i in 1:nrow(geno)){
     idx <- which(is.na(geno[i,]))
     if(length(idx)0){
     for(j in idx)
       geno[i,j] <- which.max(table(geno[,j])) - 1
       impdidx <- c(impidx, idx)
     }
     cat(i, ".örneklem kontrol ediliyor...\n")
   }
   cat("Toplam", length(impidx),"Genotip tamamlandı\n")
   return(geno)
 }



# Kod 9.2: Fenotip simülasyonu
 simulatepheno <- function(geno, h2, nq){
   n <- nrow(geno)
   nsnp <- ncol(geno)
   nt <- length(h2)
   y <- matrix(NA, nrow=n, ncol=nt)
   colnames(y) <- paste0("Trait",1:nt)
   if(nq0){
      QTL <- matrix(NA, nrow=nq, ncol=nt)
      colnames(QTL) <- paste0("Trait",1:nt)
   }
   geno <- t(geno)
   for(i in 1:nt){
     u <- rep(0, nsnp) 
     if(nq0){
       QTL[,i] <- sort(sample(1:nsnp, nq, replace=FALSE)) 
       u[QTL[,i]] <- 1  
     }
     g <- as.vector(crossprod(geno, u))    
     y[,i] <- g  rnorm(n, mean=0, 
       sd=sqrt((1-h2[i])/h2[i]*var(g, na.rm=TRUE))) 
   }
   if(nq0)
     pheno <- list(y=y, QTL=QTL)
   else
     pheno <- list(y=y)
   return(pheno)
 }





