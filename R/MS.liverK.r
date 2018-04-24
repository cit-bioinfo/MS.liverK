# Rd
# description >> internal
# argument
# item >> d >> ...
# item >> pbxgene >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
hypoxiaScore <- function(d,pbxgene=NULL){
  nam <-  names(d)
  d <- d - rowMeans(d)
  if(!is.null(pbxgene)){
    pbxgene <- pbxgene[which(pbxgene[,2] %in% c("CCNG2","EGLN3","ERO1L","WDR45L","FGF21","MAT1A","RCL1")),]
    d <- cit.dfAggregate(d[pbxgene[,1],],pbxgene[,2])
  }
  
  d1 <- colMeans(d[intersect(rownames(d),c("CCNG2","EGLN3","ERO1L","WDR45L")),])
  d2 <- colMeans(d[intersect(rownames(d),c("FGF21","MAT1A","RCL1")),])
  res <- as.numeric(d1-d2)
  resD <- c("good-Prognosis","bad-Prognosis")[1+(res>.35)]
  tmp <- factoall(as.data.frame(cbind("hypoxia.score"=res,"hypoxia.group"=resD)))
  rownames(tmp) <- nam
  attr(tmp,"available genes") <- intersect(rownames(d),c("CCNG2","EGLN3","ERO1L","WDR45L","FGF21","MAT1A","RCL1"))
  tmp
}


# Rd
# description >> internal
# argument
# item >> d >> ...
# item >> dist.meth >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
leePredictor <- function(d.,dist.meth="pearson"){
  data("prediction_Lee",  envir=sys.frame(sys.nframe()))
  gok <- setdiff(intersect(prediction_Lee$Genes,rownames(prediction_Lee$GEP)),"")
  gok <- intersect(rownames(d.),gok)
  dok <- prediction_Lee$GEP[gok,]
  aok <- prediction_Lee$ClinicalData
  
  centro <- suppressWarnings(cit.centroids( dok,aok[,2],rowCentering =median.na,dist.meth=dist.meth))
  
  xx <- suppressWarnings(cit.distToCentroids(d.[rownames(centro[[1]]$mean),],centro[[1]] , dist.meth = dist.meth, d.isPretreated = FALSE))
  
  res <- xx$pred
  
  attr(res,"available genes") <- gok
  
  res
}


# Rd
# description >> internal
# argument
# item >> d >> ...
# item >> dist.meth >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
roesslerPredictor <- function(d.,dist.meth="pearson"){    # dlda
  data("prediction_Roessler",  envir=sys.frame(sys.nframe()))
  
  gok <- intersect(rownames(d.),setdiff(intersect(prediction_Roessler$Genes,rownames(prediction_Roessler$GEP)),""))
  dok <- prediction_Roessler$GEP[gok,]
  aok <-  prediction_Roessler$ClinicalData
  
  centro <- suppressWarnings(cit.centroids( dok,aok[,2],rowCentering =median.na,dist.meth=dist.meth))
  
  xx <- suppressWarnings(cit.distToCentroids(d.[rownames(centro[[1]]$mean),],centro[[1]] , dist.meth = dist.meth, d.isPretreated = FALSE))
  
  res <- xx$pred
  
  attr(res,"available genes") <- gok
  
  res
}


# Rd
# description >> internal
# argument
# item >> d. >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
citPredictor <- function(d.){
  tmp <- list(data.frame(matrix(ncol=6,nrow=ncol(d.))))
  try(tmp <- citLiverTumorPrognosis(data=d.,data.type =c("RTQPCR.DDCTvalue","microarray.log2intensity")[2],trainingSet="boyault",more=T),silent=T)
  
  tmp
  
}


# Rd
# description >> internal
# argument
# item >> dag >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
newBoyaultPredictor <- function(dag, allHCC){
  
  data("prediction_Boyault",  envir=sys.frame(sys.nframe()))
  
  
  
  funG <- function(g){
    tmp <-  strsplit(g," /// ")
    names(tmp) <- g
    tmp <- cit.unsplit(tmp)
    colnames(tmp) <- c("single","multi")
    tmp
  }
  
  genes <- prediction_Boyault$version.2017$genes
  
  tmp=funG(genes)
  
  
  if(sum(rownames(dag) %in% genes) <  sum(rownames(dag) %in% tmp[,1])) {
    
    w = which(tmp[,1] %in% rownames(dag))
    
    dag <- cit.dfAggregate(dag[tmp[w,1],],tmp[w,2])
    
  }
  
  genes <- intersect(genes,rownames(dag))
  
  TMP <- NULL
  for(meth in c("pearson","manhattan","dlda")){
    
    suppressWarnings(cc <- cit.centroids(prediction_Boyault$version.2017$GEP[genes,], prediction_Boyault$version.2017$classes, rowCentering = median.na, dist.meth = meth))
    
    suppressWarnings(d2cc <- cit.distToCentroids(dag[genes,], cc[[1]], dist.meth =meth, maxDist = 0.5, d.isPretreated = FALSE))
    
    if(is.null(TMP)){
      
      TMP <- d2cc$pred
      
    }else{
      
      TMP <- cbind(TMP,d2cc$pred)
      
    }
    
  }
  
  p= apply(TMP,1,function(z){tt<-table(z);if(length(tt)==3){z[1]}else{names(which.max(tt))}})
  
  if(allHCC){
    
    p = c("G1" = "G1", "G2" = "G2", "G3" = "G3", "G4" = "G4", "G5" = "G5", "G6" = "G6", "NTorHCA" = "G4")[p]
    
  }
  
  names(p) <- colnames(dag)
  
  p
  
}

# Rd
# description >> internal
# argument
# item >> dag >> ...
# item >> dist.meth >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
chiangPredictor <- function(dag,dist.meth=c("dqda","cosine","pearson","dlda","euclidian")){
  data("prediction_Chiang",  envir=sys.frame(sys.nframe()))
  
  d  <- prediction_Chiang$GEP
  cl <- prediction_Chiang$classes
  gok <-prediction_Chiang$genes
  gok <- intersect(gok,rownames(dag))
  
  predChiang <- suppressWarnings(cit.centroids(d[gok,],cl,dist.meth=dist.meth[1]))
  p <- suppressWarnings(cit.distToCentroids(dag[gok,],predChiang[[1]] ,dist.meth=dist.meth[1]))
  p$genes <- gok
  
  if(length(dist.meth)>1){
    for(dm in dist.meth[-1]){
      predChiang <- suppressWarnings(cit.centroids(d[gok,],cl,dist.meth=dm))
      p. <- suppressWarnings(cit.distToCentroids(dag[gok,],predChiang[[1]] ,dist.meth=dm))
      p$pred <- cbind(p$pred,p.$pred)
    }
    p$pred <- apply(p$pred,1,function(z){tt <- table(z);names(tt)[which.max(tt)]})
    names(p$pred) <- names(p.$pred)
  }
  p
}

# Rd
# description >> internal
# argument
# item >> dag >> ...
# item >> dist.meth >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
hoshidaPredictor <- function(dag,dist.meth="cosine"){
  data("prediction_Hoshida",  envir=sys.frame(sys.nframe()))
  
  d  <- prediction_Hoshida$GEP
  cl <- prediction_Hoshida$classes
  gok <-prediction_Hoshida$genes
  gok <- intersect(gok,rownames(dag))
  
  predHoshida <- suppressWarnings(cit.centroids(d[gok,],cl,dist.meth=dist.meth))
  p <- suppressWarnings(cit.distToCentroids(dag[gok,],predHoshida[[1]] ,dist.meth=dist.meth))
  p$genes <- gok
  p
}


# Rd
# description >> internal
# argument
# item >> v >> ...
# item >> qlim >> ...
# item >> php >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
my.discr <- function(v,qlim=c(.3,.7),php=.1){
  tut <- cit.peaks(v,percentHighestPeak =php)
  if(is.null(tut[[2]])){
    gp <- c("1"="0","2"=NA,"3"="1")[cit.discretize(v,qlim,quant=T)]
  } else{
    gp <- as.numeric(v > tut[[2]][1])
  }
  if(min(table(gp))<.1*length(v)) gp <- c("1"="0","2"=NA,"3"="1")[cit.discretize(v,qlim,quant=T)]
  gp
}


# Rd
# description >> internal
# argument
# item >> dag >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
epcamAfpPredictor <- function(dag){
  EPCAM <-  my.discr(unlist(dag["EPCAM",]))
  AFP   <-  my.discr(unlist(dag["AFP",]))
  if(!"EPCAM" %in% rownames(dag) | !"AFP" %in% rownames(dag))return(NA)
  gp <- c("0 0"="-","0 1"=NA,"1 0"=NA,"1 1"="+")[paste(EPCAM,AFP)]
  names(gp) <- names(dag)
  gp
}


# Rd
# description >> internal
# argument
# item >> dag >> ...
# value >> ...
# author >> A de Reynies, F Petitprez
# keyword >> internal
# end
sigPredictor <- function(dag, scoreChoice = c("Down","-Down","Up","Up-Down")[4]){
  data("supervised_signatures",  envir=sys.frame(sys.nframe()))
  
  Lsig <- supervised_signatures
  g <- intersect(rownames(dag),unlist(Lsig))
  dag <- dag[g,]-rowMeans(dag[g,],na.rm=T)
  
  res<-lapply(names(Lsig),function(namsig){ 
    sig <- lapply(Lsig[[namsig]],intersect,g)
    if(length(sig)==1){
      tmp <- list(colMeans(dag[unlist(sig),],na.rm=T))
      names(tmp) <- namsig
      tmp
    }else{
      if(scoreChoice == "Down"){
        tmp <- colMeans(dag[unlist(sig[[1]]),],na.rm=T)
      }
      if(scoreChoice == "-Down"){
        tmp <- -colMeans(dag[unlist(sig[[1]]),],na.rm=T)
      }
      if(scoreChoice == "Up"){
        tmp <- colMeans(dag[unlist(sig[[2]]),],na.rm=T)
      }
      if(scoreChoice == "Up-Down"){
        tmp <- colMeans(dag[unlist(sig[[2]]),],na.rm=T)-colMeans(dag[unlist(sig[[1]]),],na.rm=T)
      }
    }
  })
  res <- as.data.frame(res)
  colnames(res) = names(Lsig)
  res      
}



# Rd
# description >> internal
# argument
# item >> pbxgene >> ...
# item >> sep >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
addTypeKeep <- function(pbxgene,sep=" /// "){
  
  pbxgene$type <- "single"
  
  pbxgene$type[grep(sep,pbxgene[,2],fixed=T)] <- "multiple"
  
  gsingle <- unique(pbxgene[which(pbxgene$type=="single"),2])
  
  gmultiple <-  unique(unlist(strsplit(pbxgene[which(pbxgene$type=="multiple"),2],sep,fixed=T)))
  
  gmultipleonly <-  setdiff( unique(unlist(strsplit(pbxgene[which(pbxgene$type=="multiple"),2],sep,fixed=T))) , gsingle)
  
  l<- strsplit(pbxgene[,2],sep,fixed=T)
  
  names(l) <- pbxgene[,1]
  
  tmp <- cit.unsplit(l)
  
  tmp <- cbind(tmp,"multipleOnly"=as.numeric(tmp[,1] %in% gmultipleonly))
  
  pbsm <- unique(tmp[which(tmp[,3]==1),2])
  
  pbxgene$"keep" <- 0
  
  pbxgene[which(pbxgene[,1] %in% pbsm),"keep"] <- 1
  
  pbxgene
  
}




# Rd
# description >> internal
# argument
# item >> probesData >> ...
# item >> probesSymbols >> ...
# value >> ...
# author >> F Petitprez
# keyword >> internal
# end
liverCancerSubtypes.testDataFormat <- function(probesData,probesSymbols){
  if(!(is.data.frame(probesData)&&is.data.frame(probesSymbols))){
    stop("probesData and probesSymbols must be provided as data frames. See help for further information")
  }
  if(length(rownames(probesData))==0){
    stop("probes names must be specified as row names of probesData")
  }
  if(!all(rownames(probesData) %in% probesSymbols[,1])){
    stop("one or more probe has no matching symbol")
  }
}


signatureCoverages <- function(genesList){
  data("geneUse",  envir=sys.frame(sys.nframe()))
  coverage <- unlist(lapply(colnames(geneUse), function(x){
    sigGenes <- rownames(geneUse)[which(geneUse[,x]==1)]
    return(length(which(sigGenes %in% genesList))/length(sigGenes))
  }))
  names(coverage) <- colnames(geneUse)
  zeroCov <- names(coverage[which(coverage==0)])
  warningCov <- names(coverage[intersect(which(coverage<.9),which(coverage>0))])
  if(length(zeroCov)>0){
    warning(paste("The following signatures are not represented in the genes measured. Their scores will not be calculated: ",paste(zeroCov, sep="", collapse=", "),sep=""))
  }
  if(length(warningCov)>0){
    warning(paste("The following signatures have less than 90% of their genes in the provided data. Their scores should be viewed cautiously: ",paste(warningCov, sep="", collapse=", "),sep=""))
  }
}


# Rd
# description >> gives the molecular subtype of the samples according to various classifications
# argument
# item >> probesData >> data frame containing the expression data. Probes must be in lines and samples in columns. Probes names must be specified as rownames.
# item >> probesSymbols >> data frame containing the correspondance between probes (first column) and gene symbols (second column).
# value >> ...
# author >> A de Reynies, F Petitprez
# keyword >> ...
# end
MS.liverK.subtypes <- function(probesData,probesSymbols = NULL, multipleProbeSep = " /// ", signatureChoice = c("Down","-Down","Up","Up-Down")[4], allHCC = TRUE, PrognosticSignatures = FALSE){
  
  data("prediction_Lee",  envir=sys.frame(sys.nframe()))
  data("prediction_Roessler",  envir=sys.frame(sys.nframe()))
  genesHypoxia <- c("CCNG2","EGLN3","ERO1L","WDR45L","FGF21","MAT1A","RCL1")
  genesLee <- prediction_Lee$Genes
  genesRoessler <- prediction_Roessler$Genes
  
  
  if(!is.null(probesSymbols)){
    liverCancerSubtypes.testDataFormat(probesData,probesSymbols)
    
    signatureCoverages(probesSymbols[,2])
    
    d <- probesData
    pbxgene <- probesSymbols
    
    S <- names(d)
    
    genesCIT <- c(TAF9 = "202168_at", RAMP3 = "205326_at",HN1 = "217755_at", KRT19 = "201650_at", RAN = "200750_s_at")
    
    pbxgene = addTypeKeep(pbxgene,sep = multipleProbeSep)
    
    pbxgeneok <- pbxgene
    
    try(pbxgeneok <- pbxgene[-which(pbxgene$type=="multiple" & pbxgene$"keep"==0),],silent = T)
    
    d. <- cit.dfAggregate(d[pbxgeneok[,1],],pbxgeneok[,2])
    wnona <- which(apply(d.,1,nb.na)==0)
    d. <- d.[wnona,]
    
    
    
    if(all(genesCIT %in% rownames(d))){
      dcit <- d[genesCIT,]
      rownames(dcit) <- names(genesCIT)
    }else{
      dcit <- d.
    }
  }
  else{
    warning("Expression data is assumed to be aggregated by gene symbol, as no correspondance to gene symbols was provided.")
    d. <- probesData
    dcit <- d.
    S <- names(d.)
  }
  
  
  data("prediction_Boyault",  envir=sys.frame(sys.nframe()))
  
  chiang <- hoshida <- epcam <- met <- tgfb <- boyault.new <- lee <- roe <- cit <- rep(NA,ncol(d.))
  hyp <- list(data.frame(matrix(ncol=2,nrow=ncol(d.))))
  cit <- NULL
  try( hyp       <- hypoxiaScore(d.)  , silent=T)
  try( lee       <- leePredictor(d.)  , silent=T)
  try( roe       <- roesslerPredictor(d.)  , silent=T)
  try( cit       <- citPredictor(dcit) , silent=T)
  if(is.null(cit)){
    cit <- as.data.frame(array(NA,dim=c(ncol(d.),6)))
    colnames(cit) <- c("cox.continuous","cox.discrete","pearson.medianRowCentering","pearson.noRowCentering","dqda.medianRowCentering","dqda.noRowCentering") 
    rownames(cit) <- S
    cit <- list(cit)
  }
  try( chiang    <- chiangPredictor(d.) , silent=T)
  try( hoshida   <- hoshidaPredictor(d.) , silent=T)
  try( epcam     <- epcamAfpPredictor(d.), silent=T)
  try( boyault.new <- newBoyaultPredictor(d.,allHCC), silent=T)
  
  supsig <- NULL
  try( supsig    <- sigPredictor(d., scoreChoice = signatureChoice), silent=T)
  
  if(!is.null(lee)){
    lee <- c("subgroup A"="A","subgroup B"="B")[lee]
  }
  
  prediction <- as.data.frame(cbind("Boyault" = boyault.new,
                                    "Chiang" = chiang$pred,
                                    "Hoshida" = hoshida$pred,
                                    "Lee" = lee,
                                    "Roessler" = roe,
                                    "EPCAM.AFP" = epcam))
  rownames(prediction) <- S
  
  signatures <- NULL
  if(!is.null(supsig)){
    signatures <- as.data.frame(cbind("Andersen_Cholangio_subtype1" = supsig$ANDERSEN_CHOLANGIO_SUBTYPE1,
                                      "Andersen_KRT19_AvB" = supsig$ANDERSEN_KRT19_AvB,
                                      "Boyault_G1" = supsig$`Boyault - up in G1`,
                                      "Boyault_G2" = supsig$`Boyault - up in G2`,
                                      "Boyault_G3" = supsig$`Boyault - up in G3`,
                                      "Boyault_G4" = supsig$`Boyault - up in G4`,
                                      "Boyault_G5" = supsig$`Boyault - up in G5`,
                                      "Boyault_G6" = supsig$`Boyault - up in G6`,
                                      "Boyault_G123" = supsig$`Boyault - up in G123`,
                                      "Boyault_G456" = supsig$`Boyault - up in G456`,
                                      "Boyault_G12" = supsig$`Boyault - up in G12`,
                                      "Boyault_G13" = supsig$`Boyault - up in G13`,
                                      "Boyault_G23" = supsig$`Boyault - up in G23`,
                                      "Boyault_G56" = supsig$`Boyault - up in G56`,
                                      "Cairo_hepatoblastoma_classes" = supsig$CAIRO_HEPATOBLASTOMA_CLASSES,
                                      "Chiang_CTNNB1" = supsig$`Chiang - CTNNB1`,
                                      "Chaing_prolif" = supsig$`Chiang - Prolif`,
                                      "Chiang_Inflam" = supsig$`Chiang - Inflam`,
                                      "Chiang_poly7" = supsig$`Chiang - poly7`,
                                      "Chiang_unannot" = supsig$`Chiang - unannot`,
                                      "Hoshida_S1" = supsig$HOSHIDA_LIVER_CANCER_SUBCLASS_S1,
                                      "Hoshida_S2" = supsig$HOSHIDA_LIVER_CANCER_SUBCLASS_S2,
                                      "Hoshida_S3" = supsig$HOSHIDA_LIVER_CANCER_SUBCLASS_S3,
                                      "Lee_A" = supsig$`Lee - up in A`,
                                      "Lee_B" = supsig$`Lee - up in B`))
    rownames(signatures) <- S
    wnona <- which(apply(signatures,2,nb.na)==0)
    signatures <- signatures[,wnona]
  }
  
  molecularSubtype <- list(prediction = prediction, signatures = signatures)
  
  prediction <- as.data.frame(cbind("Nault" = cit[[1]][3]))
  colnames(prediction) <- c("Nault")
  rownames(prediction) <- S
  
  signatures <- NULL
  if(PrognosticSignatures){
    if(!is.null(supsig)){
      signatures <- as.data.frame(cbind("Budhu_liverK_metastasis" = supsig$BUDHU_LIVER_CANCER_METASTASIS,
                                        "Hao_survival" = -supsig$HAO_SURVIVAL,
                                        "Hoshida_liverK_late_recurrence" = supsig$HOSHIDA_LIVER_CANCER_LATE_RECURRENCE,
                                        "Hoshida_liverK_survival" = supsig$HOSHIDA_LIVER_CANCER_SURVIVAL,
                                        "Iizuka_liverK_early_recurrence" = supsig$IIZUKA_LIVER_CANCER_EARLY_RECURRENCE,
                                        "Kim_poor_survival" = supsig$KIM_POOR_SURVIVAL,
                                        "Kurokawa_liverK_early_recurrence" = -supsig$KUROKAWA_LIVER_CANCER_EARLY_RECURRENCE,
                                        "Lee_liverK_survival" = -supsig$LEE_LIVER_CANCER_SURVIVAL,
                                        "Minguez_vascular_invasion" = supsig$MINGUEZ_VASCULAR_INVASION,
                                        "Okamoto_liverK_multicentric_occurrence" = supsig$OKAMOTO_LIVER_CANCER_MULTICENTRIC_OCCURRENCE,
                                        "Roessler_liverK_metastasis" = supsig$ROESSLER_LIVER_CANCER_METASTASIS,
                                        "Wang_recurrent_liverK" = supsig$WANG_RECURRENT_LIVER_CANCER,
                                        "Woo_Cholangio-like" = supsig$WOO_CHOLANGIO.LIKE,
                                        "Woo_liverK_recurrence" = supsig$WOO_LIVER_CANCER_RECURRENCE,
                                        "Ye_metastatic_liverK" = supsig$YE_METASTATIC_LIVER_CANCER,
                                        "Yoshioka_liverK_early_recurrence" = supsig$YOSHIOKA_LIVER_CANCER_EARLY_RECURRENCE))
      rownames(signatures) <- S
      wnona <- which(apply(signatures,2,nb.na)==0)
      signatures <- signatures[,wnona]
    }
  }
  
  
  prognosis <- list(prediction = prediction, signatures = signatures)
  
  biologicalPathwaysSignatures <- NULL
  if(!is.null(supsig)){
    biologicalPathwaysSignatures <- as.data.frame(cbind("Coulouarn_temporal_TGFB1_signature" = supsig$COULOUARN_TEMPORAL_TGFB1_SIGNATURE,
                                                        "Kaposi_liverK_met" = supsig$KAPOSI_LIVER_CANCER_MET,
                                                        "Oishi_Cholangio_stem_cell_like" = supsig$OISHI_CHOLANGIO_STEM_CELL_LIKE,
                                                        "Yamashita_liverK_stem_cell" = supsig$YAMASHITA_LIVER_CANCER_STEM_CELL,
                                                        "Yamashita_liverK_with_EPCAM" = supsig$YAMASHITA_LIVER_CANCER_WITH_EPCAM))
    rownames(biologicalPathwaysSignatures) <- S
    wnona <- which(apply(biologicalPathwaysSignatures,2,nb.na)==0)
    biologicalPathwaysSignatures <- biologicalPathwaysSignatures[,wnona]
  }
  
  biologicalPathwaysSignatures = as.data.frame(cbind(biologicalPathwaysSignatures,
                                                     "Van_Malenstein_hypoxia" = hyp[[1]]))
  rownames(biologicalPathwaysSignatures) <- S
  
  result = list(molecularSubtype = molecularSubtype, prognosis = prognosis, biologicalPathwaysSignatures = biologicalPathwaysSignatures)
  
  result
}

# Rd
# description >> draws a segment plot representing the classes from liverCancerSubtypes()
# argument
# item >> subtypes >> result from liverCancerSbtypes()
# item >> col >> named vector containing the colors to use for each class. Names must contain each of "A","B","S1","S2","S3","CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated","G1","G2","G3","G4","G5","G6".
# value >> ...
# author >> F Petitprez
# keyword >> ...
# end
MS.liverK.plot <- function(subtypes, col=NULL){
  if(is.null(col)){
    
    allHCC = length(unique(subtypes$molecularSubtype$prediction$Boyault))<=6
    if(allHCC){
      col <- c("orange","yellow","orange","yellow","red","pink","yellow","pink","orange","red","yellow","green","darkgreen","green","lightblue","purple","orange","yellow","orange","yellow")
      names(col) <- c("A","B","subgroup A","subgroup B","S1","S2","S3","G1","G2","G3","G4","G5","G6","CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated","+","-")
    }else{
      col <- c("orange","yellow","orange","yellow","red","pink","yellow","pink","orange","red","yellow","green","darkgreen","grey","green","lightblue","purple","orange","yellow","orange","yellow")
      names(col) <- c("A","B","subgroup A","subgroup B","S1","S2","S3","G1","G2","G3","G4","G5","G6","NTorHCA","CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated","+","-")
    }
  }
  #reorder
  df <- subtypes$molecularSubtype$prediction[c("Lee","Hoshida","Boyault","Chiang","Roessler","EPCAM.AFP")]
  ref <- 1
  expr <- ("df <- df[order(")
  expr <- paste(expr, "df[,", ref, "],", sep = "")
  if (length(setdiff(1:(ncol(df) - 1), ref)) > 0) {
    for (c in setdiff(1:(ncol(df) - 1), ref)) expr <- paste(expr, 
                                                            "df[,", c, "],", sep = "")
    expr <- paste(expr, "df[,", c + 1, "]),]", sep = "")
  }
  else {
    c <- setdiff(1:(ncol(df)), ref)
    expr <- paste(expr, "df[,", c, "]),]", sep = "")
  }
  eval(parse(text = expr))
  samplesOrder = rownames(df)
  
  #plot
  if(!is.null(subtypes[[2]]$signatures)){
    layout(matrix(c(1,2,3,4),4,1),heights = c(25,35,26,14))
    cit.dfSegmentplot(df, labelscolors = col,reorder = F, main = "Molecular subtypes",lwd=1)
    x1 = (95*par("mai")[2])/(100*par("pin")[1])
    x2 = (70*par("mai")[2])/(100*par("pin")[1])
    x3 = (60*par("mai")[2])/(100*par("pin")[1])
    y1 = 0
    y2 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(3.5/20)/(par("pin")[2])
    y3 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(7/20)/(par("pin")[2])
    y4 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(11.5/20)/(par("pin")[2])
    y5 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(1/20)/(par("pin")[2])
    y6 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(7.5/20)/(par("pin")[2])
    legend("topleft", legend = c("A","B"), fill = col[c("A","B")], title = "Lee", bty="n",xpd=T, inset = c(-x1,y1), cex = .8)
    legend("topleft", legend = c("A", "B"), fill = col[c("subgroup A","subgroup B")], title = "Roessler", bty="n",xpd=T, inset = c(-x1,y2), cex = .8)
    legend("topleft", legend = c("S1", "S2","S3"), fill = col[c("S1", "S2","S3")], title = "Hoshida", bty="n",xpd=T, inset = c(-x1,y3), cex = .8)
    legend("topleft", legend = c("+","-"), fill = col[c("+","-")], title = "EPCAM.AFP", bty="n",xpd=T, inset = c(-x1,y4), cex = .8)
    legend("topleft", legend = c("CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated"), fill = col[c("CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated")], title = "Chiang", bty="n",xpd=T, inset = c(-x2,y5), cex = .8)
    if(allHCC){
      legend("topleft", legend = c("G1","G2","G3","G4","G5","G6"), fill = col[c("G1","G2","G3","G4","G5","G6")], title = "Boyault", bty="n",xpd=T, inset = c(-x3,y6), cex = .8)
    }else{
      legend("topleft", legend = c("G1","G2","G3","G4","G5","G6","NT or HCA"), fill = col[c("G1","G2","G3","G4","G5","G6","NTorHCA")], title = "Boyault", bty="n",xpd=T, inset = c(-x3,y6), cex = .8)
    }
    cit.image(subtypes$molecularSubtype$signatures[samplesOrder,ncol(subtypes$molecularSubtype$signatures):1],labColPos = NA, main = "Molecular subtypes signatures")
    cit.image(cbind(subtypes$prognosis$signatures[,ncol(subtypes$prognosis$signatures):1],subtypes$prognosis$prediction)[samplesOrder,],labColPos = NA, main = "Prognosis", colpos = "orangered", colneg = "royalblue")
    cit.image(subtypes$biologicalPathwaysSignatures[samplesOrder,ncol(subtypes$biologicalPathwaysSignatures):1],labColPos = NA, main = "Biological pathways signatures")
  } else{
    layout(matrix(c(1,2,3,4),4,1),heights = c(25,35,8,14))
    #layout(matrix(c(1,2,3),3,1),heights = c(25,35,14))
    cit.dfSegmentplot(df, labelscolors = col,reorder = F, main = "Molecular subtypes",lwd=1)
    x1 = (95*par("mai")[2])/(100*par("pin")[1])
    x2 = (70*par("mai")[2])/(100*par("pin")[1])
    x3 = (60*par("mai")[2])/(100*par("pin")[1])
    y1 = 0
    y2 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(3.5/20)/(par("pin")[2])
    y3 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(7/20)/(par("pin")[2])
    y4 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(11.5/20)/(par("pin")[2])
    y5 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(1/20)/(par("pin")[2])
    y6 = (par("pin")[2]+par("mai")[1]+par("mai")[3])*(7.5/20)/(par("pin")[2])
    legend("topleft", legend = c("A","B"), fill = col[c("A","B")], title = "Lee", bty="n",xpd=T, inset = c(-x1,y1), cex = .8)
    legend("topleft", legend = c("A", "B"), fill = col[c("subgroup A","subgroup B")], title = "Roessler", bty="n",xpd=T, inset = c(-x1,y2), cex = .8)
    legend("topleft", legend = c("S1", "S2","S3"), fill = col[c("S1", "S2","S3")], title = "Hoshida", bty="n",xpd=T, inset = c(-x1,y3), cex = .8)
    legend("topleft", legend = c("+","-"), fill = col[c("+","-")], title = "EPCAM.AFP", bty="n",xpd=T, inset = c(-x1,y4), cex = .8)
    legend("topleft", legend = c("CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated"), fill = col[c("CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated")], title = "Chiang", bty="n",xpd=T, inset = c(-x2,y5), cex = .8)
    if(allHCC){
      legend("topleft", legend = c("G1","G2","G3","G4","G5","G6"), fill = col[c("G1","G2","G3","G4","G5","G6")], title = "Boyault", bty="n",xpd=T, inset = c(-x3,y6), cex = .8)
    }else{
      legend("topleft", legend = c("G1","G2","G3","G4","G5","G6","NT or HCA"), fill = col[c("G1","G2","G3","G4","G5","G6","NTorHCA")], title = "Boyault", bty="n",xpd=T, inset = c(-x3,y6), cex = .8)
    }
    cit.image(subtypes$molecularSubtype$signatures[samplesOrder,ncol(subtypes$molecularSubtype$signatures):1],labColPos = NA, main = "Molecular subtypes signatures")
    cit.image(cbind("Nault" = subtypes$prognosis$prediction[samplesOrder,]," " = subtypes$prognosis$prediction[samplesOrder,]),labColPos = NA, colpos = "orangered", colneg = "royalblue", main = "Prognosis")
    #cit.image(matrix(subtypes$prognosis$prediction$Nault[samplesOrder],ncol=1,dimnames = list(NULL,"Nault")),labColPos = NA, colpos = "orangered", colneg = "royalblue", main = "Prognosis")
    cit.image(subtypes$biologicalPathwaysSignatures[samplesOrder,ncol(subtypes$biologicalPathwaysSignatures):1],labColPos = NA, main = "Biological pathways signatures")
  }
  return(NULL)
}




MS.liverK.discretize <- function(subtypes){
  S <- rownames(subtypes$molecularSubtype$prediction)
  discSubtypes <- as.data.frame(apply(subtypes$molecularSubtype$signatures,2,cit.discretize,"lim"=c(.2,.4,.6,.8),"quant"=T))
  rownames(discSubtypes) <- S
  discProg <- NULL
  if(!is.null(subtypes$prognosis$signatures)){
    discProg <- as.data.frame(apply(subtypes$prognosis$signatures,2,cit.discretize,"lim"=c(.2,.4,.6,.8),"quant"=T))
    rownames(discProg) <- S
  }
  discPath <- as.data.frame(apply(subtypes$biologicalPathwaysSignatures,2,cit.discretize,"lim"=c(.2,.4,.6,.8),"quant"=T))
  rownames(discPath) <- S
  molecularSubtype_discretized = list(prediction = subtypes$molecularSubtype$prediction, signatures = discSubtypes)
  prognosis_discretized = list(prediction = subtypes$prognosis$prediction, signatures = discProg)
  biologicalPathwaysSignatures_discretized = discPath
  res = list(molecularSubtype_discretized=molecularSubtype_discretized,prognosis_discretized=prognosis_discretized,biologicalPathwaysSignatures_discretized=biologicalPathwaysSignatures_discretized)
  return(res)
}


MS.liverK.convert <- function(userData){
  
  if(ncol(userData)<50) warning("This function is designed to convert datasets of at least 50 samples. Below this limit, the results are subject to caution.")
  
  data("GSE20238",  envir=sys.frame(sys.nframe()))
  
  expGSE20238 = cit.dfAggregate(GSE20238$GEP,GSE20238$AnnotProbeset[rownames(GSE20238$GEP),"Gene.Symbol"])
  
  commonGenes = sort(intersect(rownames(userData),rownames(expGSE20238)))
  
  expGSE20238 = expGSE20238[commonGenes,]
  toConvert = userData[commonGenes,]
  
  sdGSE20238 = apply(expGSE20238,1,sd,na.rm=T)
  sdUser = apply(toConvert,1,sd,na.rm=T)
  
  genesUsed = intersect(commonGenes[which(sdGSE20238 > quantile(sdGSE20238,.1))],
                        commonGenes[which(sdUser > quantile(sdUser,.1))])
  
  
  expGSE20238 = expGSE20238[genesUsed,]
  toConvert = userData[genesUsed,]
  
  arrayConv = data.frame(matrix(nrow = length(genesUsed),ncol = 2))
  rownames(arrayConv) = genesUsed
  colnames(arrayConv) = c("x","intercept")
  
  for(g in genesUsed){
    y=quantile(expGSE20238[g,],probs = seq(0,1,.10),na.rm = T)
    x=quantile(toConvert[g,],probs = seq(0,1,.10),na.rm = T)
    model = lm(y~x)
    arrayConv[g,] = c(model$coefficients[2], model$coefficients[1])
  }
  
  converted = data.frame(t(apply(cbind(arrayConv,toConvert),1,function(line){
    x = line[3:length(line)]
    return(c(rep(NA,2),line[1]*x+line[2]))
  })))[,3:(ncol(toConvert)+2)]
  
  return(converted)
  
}



