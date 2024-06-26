# SIMER ile sim�lasyon
#
if(!require(simer)){
  if(!require(devtools))
    install.packages("devtools", repos="https://cloud.r-project.org")
  devtools::install_github("xiaolei-lab/SIMER")
}
library(simer)

# Dosyalar�n tan�mlanmas�
dataFolder <- "D:/gwasbook/data/simdata"
if (!file.exists(dataFolder))
   dir.create(dataFolder)
dataFile <- "cattle"
dataFormat <- "plink"

# Senaryo parametreleri

nInd <- 1000 			# Birey say�s�
nSnp <- 5000			# SNP say�s�
nChr <- 29			# Kromozom say�s�
nChrLen <- 1.5e+08		# Kromozom uzunlu�u
nQtl <- c(600,500,600,400)	# Karakterlere g�re sim�le edilecek QTL say�s�
dQtl <- c("norm","norm","norm","norm") # Karakterlere g�re sim�le edilecek QTL da��l��lar� (norm', 'geom', 'gamma', 'beta)
mQtl <- c("A","A","A", "A")   # Karakterlere g�re sim�le edilecek QTL model t�r� (A eklemeli, D dominans)
mutRate <- c(1e-8, 1e-8)	# Mutasyon oranlar� (SNP ve QTL i�in)
speciesName <- "cattle"		# T�rler: "arabidopsis", "cattle", "chicken", "dog", "horse", "human", "maize", "mice", "pig",  "rice".
h2a <- c(0.1, 0.3, 0.5, 0.7)    # Dar anlamda kal�t�m dereceleri
h2d <- c(0.01, 0.05, 0.15, 0.25)    # Dominans kal�t�m dereceleri
pVar <- c(28130, 2000, 0.137, 50)   # Fenotipik varyanslar

# Genel parametreler (atla)
SP <- NULL
SP <- param.global(
  seed.sim = 123,		# Rastlant�sal say� �ekirde�i
  out = dataFile,		# ��k�� dosyalar�n�n ad�
  outpath = dataFolder,		# ��k�� kla�r� ad�
  out.format = dataFormat,	# ��k�� format� ('numeric' veya 'plink')
  #pop.gen = 5,		 	# Sim�le edilen populasyon say�s�
  #out.geno.gen = 1:5,		# Genotip verisinin generasyon say�s� 
  #out.pheno.gen = 1:5,	# Fenotip verisinin generasyon say�s� 
  #useAllGeno = TRUE,		# Fenotipleri sim�le etmek i�in t�m genotipleri kullanma durumu
  #ncpus=NULL,			# Kullan�lacak CPU �ekirdek say�s�
  verbose = TRUE		# Gevezelik durumu
)


# Genotipik map parametreleri
popMap <- generate.map(
  pop.marker = nSnp, 
  num.chr = nChr,
  len.chr = nChrLen)

# Annotation parametreleri
SP <- param.annot(
   SP = SP,
   pop.map = popMap,  		# A�a��daki species kullan�lmadan jenerik map kullan�lacaksa devreye al
   #species = speciesName, 	# pop.map tan�ml� ise gerek yok
   #pop.marker = nSnp,
   num.chr = nChr,
   len.chr = nChrLen,
   qtn.var = list(tr1=1, tr2=1, tr3=1, tr4=1),
   qtn.num = list(tr1=nQtl[1], tr2=nQtl[2], tr3=nQtl[3], tr4=nQtl[4]),
   qtn.dist = list(tr1=dQtl[1], tr2=dQtl[2], tr3=dQtl[3], tr4=dQtl[4]),
   qtn.model = list(tr1=mQtl[1], tr2=mQtl[2], tr3=mQtl[3], tr4=mQtl[4])
)

# Genotipik parametreler
SP <- param.geno(
   SP = SP, 
   pop.marker = nSnp, 
   pop.ind = nInd,
   incols = 1,
   rate.mut=mutRate) 

# �evre fakt�rleri
popEnvironment <- list(
  F1 = list( # Sabit etki 1
    level = c("L1", "L2"),
    effect = list(tr1 = c(0.1, 0.2), tr2=c(0.2, 0.3), tr3=c(0.1, 0.2), tr4=c(0.05, 0.15))
  ), 
  F2 = list( # Sabit etki 2
    level = c("S1", "S2", "S3"),
    effect = list(tr1 = c(0.1, 0.2, 0.3), tr2=c(0.1, 0.2, 0.3), tr3=c(0.1, 0.2, 0.3), tr4=c(0.1, 0.2, 0.3))
  ),
  C1 = list( # Kovaryet 1
    level = c(0.70, 0.75, 0.90),
    slope = list(tr1 = 0.2, tr2=0.4, tr3=0.2, tr4=0.5)
  ),
  R1 = list( # Rastlant�sal etki 1
    level = c("R1", "R2", "R3"),
    ratio = list(tr1 = 0.1, tr2=0.1, tr3=0.1 , tr4=0.1)
  )
)

# Fenotipik parametreler
SP <- param.pheno(
   SP = SP, 
   #pop = pop,
   pop.env = popEnvironment,
   phe.type = list(tr1 = list(case = 0.2, control = 0.8), tr2 = "continuous", tr3 = "continuous", tr4 = "continuous"),
   phe.model = list(
     tr1 = "Trait1 = A + F1 + F2 + C1 + R1 + E",  # �evresel etkiler dahil edilmeyecekse yaln�z A+E
     tr2 = "Trait2 = A + F1 + F2 + C1 + R1 + E",  
     tr3 = "Trait3 = A + F1 + F2 + C1 + R1 + E",  
     tr4 = "Trait4 = A + F1 + F2 + C1 + R1 + E"),
    # phe.h2D = list(tr1=h2d[1], tr2=h2d[2], tr3=h2d[3], tr4=h2d[4]),  #Dominans kal�t�m dereceleri (varsa)
    # phe.var = list(tr1=pVar[1], tr2=pVar[2], tr3=pVar[3], tr4=pVar[4]),
    phe.h2A = list(tr1=h2a[1], tr2=h2a[2], tr3=h2a[3], tr4=h2a[4]),  #Kal�t�m dereceleri
    phe.corA = matrix(c( # Eklemeli genetik korelasyonlar
    1, 0, 0, 0,
    0, 1, 0, 0,
    0, 0, 1, 0,
    0, 0, 0, 1),
   nrow=4, ncol=4) # Korelasyon yok, ili�kisiz
)

# Sim�lasyon uygula
seed <- 28
set.seed(seed)
SP <- annotation(SP, verbose = TRUE)
set.seed(seed)
SP <- genotype(SP, verbose = TRUE)
set.seed(seed)
SP <- phenotype(SP, verbose = TRUE)

# QTL ve Map inceleme
# Map  yap�s� ve i�eri�i 
map <- SP$map$pop.map
dim(map)
map[1:10, ]

# Simer'in rastlant�sal olarak atad��� QTN etki de�erleri (a)
map[!is.na(map$QTN1_A),]
map[!is.na(map$QTN2_A),]
map[!is.na(map$QTN3_A),]
map[!is.na(map$QTN4_A),]

# Simer'in rastlant�sal olarak belirledi�i QTL indisleri
qtnIdx1 <- SP$map$qtn.index$tr1
qtnIdx2 <- SP$map$qtn.index$tr2
qtnIdx3 <- SP$map$qtn.index$tr3
qtnIdx4 <- SP$map$qtn.index$tr4

# Etkileri de�i�tir
a <- 1 #QTN de�eri, NULL atan�rsa Simer'in atad��� etkiler kullan�l�r
a <- NULL
if(!is.null(a)){
  avec1 <- c(-2*a,-a,-0.5*a, 0.5*a, a, 2*a) #Senaryo 1
  avec2 <- c(a, a, a, a, a, a) #Senaryo 2
  avec3 <- c(-a, -a, -a, -a, -a, -a) #Senaryo 3
  avec4 <- c(2*a, 2*a, 2*a, 2*a, 2*a, 2*a) #Senaryo 4
  pvec <- rep(1/6, 6)
  
  # Senaryo uygula
  avec <- avec2
  map$QTN1_A[qtnIdx1] <- sample(avec, nQtl[1], replace=TRUE, prob=pvec)
  map$QTN2_A[qtnIdx2] <- sample(avec, nQtl[2], replace=TRUE, prob=pvec)
  map$QTN3_A[qtnIdx3] <- sample(avec, nQtl[3], replace=TRUE, prob=pvec)
  map$QTN4_A[qtnIdx4] <- sample(avec, nQtl[4], replace=TRUE, prob=pvec)
}

# QTL indisleri tablosu
qtlIdx <- cbind(qtnIdx1, qtnIdx2, qtnIdx3, qtnIdx4)
colnames(qtlIdx) <- c("Trait1", "Trait2", "Trait3", "Trait4")
qtlIdx
qtlVal <- cbind(map$QTN1_A[qtnIdx1], map$QTN2_A[qtnIdx2], map$QTN3_A[qtnIdx3], map$QTN4_A[qtnIdx4])
qtlVal

# Fenotip matrisi
pheno <- SP$pheno$pop$gen1  
dim(pheno)
str(pheno)
pheno[1:10, 1:18]
pheno1 <- pheno[,c(1:15)]
head(pheno1)

# Genotip (SNP) matrisi
geno <- SP$geno$pop.geno$gen1
options(bigmemory.allow.dimnames=TRUE)
dim(geno)
geno[1:10, 1:10]

# Genotipi 0,1,2 kodlu matrise d�n��t�r
geno1 <- simer::geno.cvt1(geno)  # Convert genotype matrix from (0, 1) to (0, 1, 2).
# Genotipi 0,1,2 kodundan 0,1 kodlu matrise d�n��t�r
geno1_2 <- simer::geno.cvt2(geno1)   # Convert genotype matrix from  (0, 1, 2) to (0, 1)

#rownames(geno1) <- paste0("snp", 1:nrow(geno))
colnames(geno1) <- pheno_1$index
rownames(geno1) <- map$SNP

geno1[1:10, 1:10]
geno1_2[1:10, 1:10]

#------------------------------------------------------------------------------------------
# Kurucu populasyon verilerini kaydet
setwd(dataFolder)
write.table(geno1, file="cattle_sim_geno1.dat", sep="\t", row.names=T, col.names=T, quote=F)
write.table(map, file="cattle_sim.map", sep="\t", row.names=F, col.names=T, quote=F)
write.table(pheno1, file="cattle_sim_pheno1.dat", sep="\t", row.names=F, col.names=T, quote=F)
write.table(qtlIdx, file="cattle_sim_qtlidx1.dat", sep="\t", row.names=F, col.names=T, quote=F)
write.table(qtlVal, file="cattle_sim_qtlval1.dat", sep="\t", row.names=F, col.names=T, quote=F)

# �kinci generasyon
# Tek �zellik i�in seleksiyon parametresi (bireysel seleksiyon)
# SP <- param.sel(
#   SP = SP, 
#   sel.single = "ind",  # 'ind', 'fam', 'infam' ve 'comb' olabilir
#   sel.crit = "TBV", #seleksiyon �l��t� (TBV, TGV ve pheno olabilir)
#   ps=c(0.1, 0.8),   # Erkek ve di�i seleksiyon oran�
# )

# Birden fazla karakter i�in seleksiyon parametresi (�ndeks seleksiyonu)
SP <- param.sel(
   SP = SP, 
   sel.multi = "tdm",  # 'index', 'indcul' ve 'tdm' olabilir
   sel.crit = "TBV",
   index.wt = c(1, 1, 1, 1), # Indeks seleksiyonu karakterlerin ekonomik a��rl�klar�
   index.tdm = c(3, 1, 2, 4), #Tandem seleksiyon i�in �zelliklerin indisleri - otomatik b�rak�lmal�
   ps=c(0.1, 0.8), #Cinsiyetlere g�re se�ilenlerin y�zdesi
   #goal.perc = 0.1, # Ortalamas�ndan daha fazla olmas� ama�lanan birey y�zdesi
   #pass.perc = 0.9 # Beklenen m�kemmel bireylerin y�zdesi
   decr = TRUE  # Seleksiyon �l��t�n� b�y�kten k����e s�ralayarak se�im
)

# Seleksiyon uygula
SP <- selects(SP)

# �reme parametreleri ayarla
SP <- param.reprod(
  SP = SP, 
  pop.gen = 5,  #Sim�lasyon yap�lacak generasyonlar�n say�s�
  reprod.way = "assort",  #randmate, 'clone', 'dh', 'selfpol', 'randmate', 'randexself', 'assort', 'disassort', '2waycro', '3waycro', '4waycro', 'backcro', and 'userped'.
  sex.rate = 0.5,   	  # Erkek oran�
  prog = 1.07 		  # �iftle�me ba��na d�l say�s� (%7 ikizlik oldu�unu varsay�ld�)
)

# �iftle�tir
SP <- reproduces(SP)

# Generasyonlara ait verilerin al�nmas�

# Fenotip 2. generasyon
pheno <- SP$pheno$pop$gen2  
pheno2 <- pheno[,c(1:15)]
head(pheno2)
dim(pheno2)

# Fenotip 3. generasyon
pheno <- SP$pheno$pop$gen3  
pheno3 <- pheno[,c(1:15)]
head(pheno3)
dim(pheno3)

# Fenotip 4. generasyon
pheno <- SP$pheno$pop$gen4  
pheno4 <- pheno[,c(1:15)]
head(pheno4)
dim(pheno4)

# Fenotip 5. generasyon
pheno <- SP$pheno$pop$gen5  
pheno5 <- pheno[,c(1:15)]
head(pheno5)
dim(pheno5)

# Genotip 2. generasyon
geno <- SP$geno$pop.geno$gen2
geno2 <- simer::geno.cvt1(geno)  
colnames(geno2) <- pheno2$index
rownames(geno2) <- map$SNP
dim(geno2)
geno2[1:10, 1:10]

# Genotip 3. generasyon
geno <- SP$geno$pop.geno$gen3
geno3 <- simer::geno.cvt1(geno)  
colnames(geno3) <- pheno3$index
rownames(geno3) <- map$SNP
dim(geno3)
geno3[1:10, 1:10]

# Genotip 4. generasyon
geno <- SP$geno$pop.geno$gen4
geno4 <- simer::geno.cvt1(geno)  
colnames(geno4) <- pheno4$index
rownames(geno4) <- map$SNP
dim(geno4)
geno4[1:10, 1:10]

# Genotip 5. generasyon
geno <- SP$geno$pop.geno$gen5
geno5 <- simer::geno.cvt1(geno)  
colnames(geno5) <- pheno5$index
rownames(geno5) <- map$SNP
dim(geno5)
geno5[1:10, 1:10]

allgeno <- rbind(t(geno1), t(geno2), t(geno3), t(geno4), t(geno5))
allgeno <- t(allgeno)
head(allgeno)

allpheno <- rbind(pheno1,pheno2,pheno3,pheno4,pheno5)
head(allpheno)

setwd(dataFolder)
write.table(geno2, file="cattle_sim_geno2.dat", sep="\t", row.names=F, col.names=F, quote=F)
write.table(pheno2, file="cattle_sim_pheno2.dat", sep="\t", row.names=F, col.names=T, quote=F)
write.table(geno3, file="cattle_sim_geno3.dat", sep="\t", row.names=F, col.names=F, quote=F)
write.table(pheno3, file="cattle_sim_pheno3.dat", sep="\t", row.names=F, col.names=T, quote=F)
write.table(geno4, file="cattle_sim_geno4.dat", sep="\t", row.names=F, col.names=F, quote=F)
write.table(pheno4, file="cattle_sim_pheno4.dat", sep="\t", row.names=F, col.names=T, quote=F)
write.table(geno5, file="cattle_sim_geno5.dat", sep="\t", row.names=F, col.names=F, quote=F)
write.table(pheno5, file="cattle_sim_pheno5.dat", sep="\t", row.names=F, col.names=T, quote=F)
write.table(allgeno, file="cattle_sim_geno_all.dat", sep="\t", row.names=F, col.names=F, quote=F)
write.table(allpheno, file="cattle_sim_pheno_all.dat", sep="\t", row.names=F, col.names=T, quote=F)

#--------------------------------------------------------------------------------

# Kurucu pop�lasyonu olu�turma �rne�i
# �lk %10 hayvan erkek, di�erleri di�i
pop <- generate.pop(pop.ind = nInd, from = 1, ratio = 0.1)
#table(pop$sex)