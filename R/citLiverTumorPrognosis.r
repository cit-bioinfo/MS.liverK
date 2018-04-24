# Rd
# description >> predicting prognosis of hepatocellular carcinomas based on 5 transcripts RTQPCR measurements 
# argument
# item >> data >>    dataframe, with rows=transcripts measurements and columns=(test set) samples  
# item >> trainingSet >>  character string, designating the training set to be used, either  "nault" (default) or "boyault" , respectively corresponding to RTQPCR ddCt  and Affymetrix HG U133A RMA log2 data
# item >> more >> boolean, if TRUE additional predictors are calculated and returned (default : FALSE => only the Cox model derived predictor is returned)
# item >> genes >>  named vector of characters, with elements corresponding to Affymetrix HG-U133A probe set ids and names  to related HUGO Gene Symbols. This parameter is given for internal use and should not be modified by the user.
# value >> dataframe , rows=samples, columns=prognostic predicted group according to five methods
# author >> A de Reynies
# keyword >> methods                                           
# end
citLiverTumorPrognosis <- function( data=NULL,
                                    data.type = c("RTQPCR.DDCTvalue","microarray.log2intensity")[1],
                                    trainingSet="boyault",
                                    more=FALSE,
                                    genes=c("TAF9"="202168_at","RAMP3"="205326_at","HN1"="217755_at","KRT19"="201650_at","RAN"="200750_s_at")
                                   ){
       
       LPRED <- list()                            
                                   
       if(!all(names(genes) %in% rownames(data))) {
         stop("Error - function citLiverTumorPrognosis : predictive genes are not found in data rownames")
       }else{
         extdata <- data[names(genes),]    
       } 

       data("prognosis_Nault",  envir=sys.frame(sys.nframe()))
       traindata  <- prognosis_Nault$GEP
       rownames(traindata) <- names(genes)
       trainannot <- prognosis_Nault$ClinicalData
       
       if((data.type=="RTQPCR.DDCTvalue" & trainingSet!="nault") | 
          (data.type=="microarray.log2intensity" & trainingSet!="boyault")){
          #  cat("m")
            data <- -data
       }    

       
       trainda<-as.data.frame(cbind(t(traindata),trainannot[,c("OSS","OSS.delay")]))
       model <- mycoxph( datannot=trainda,
                          vars=names(trainda)[1:length(genes)],
                          col.event="OSS",
                          col.delai="OSS.delay",
                          cutDelai=NULL,
                          silent=TRUE,
                          meth="breslow")
       TEMP <-   predict.coxph.adr(model$coxModel,  t(extdata) )
       cox.continuous <- TEMP$pred
       LPRED <- c(LPRED,list("COX"=TEMP$formul))
       
       if(abs(median(cox.continuous,na.rm=TRUE))<1){
           coxtmp <- cit.discretize( cox.continuous,0,quant=FALSE)
       }else{
           coxtmp <- cit.discretize( cox.continuous,0.5,quant=TRUE) 
       }
       if(more){
             pred <- cit.centroids(traindata,trainannot[,"OSS"],dist.meth="pearson",rowCentering =median.na)
             LPRED <- c(LPRED,list("CENTROIDS.medCen"=pred$centroids[1:4]))
             pearson.medianRowCentering <- 1+as.numeric(cit.distToCentroids(extdata,pred[[1]],d.isPretreated=F,dist.meth="pearson")$pred)
             
             pred <- cit.centroids(traindata,trainannot[,"OSS"],dist.meth="pearson",rowCentering =NA)
             LPRED <- c(LPRED,list("CENTROIDS.uncen"=pred$centroids[1:4]))
             pearson.noRowCentering <- 1+as.numeric(cit.distToCentroids(extdata,pred[[1]],d.isPretreated=F,dist.meth="pearson")$pred)
      
             pred <- cit.centroids(traindata,trainannot[,"OSS"],dist.meth="dqda",rowCentering =median.na)
             dqda.medianRowCentering <- 1+as.numeric(cit.distToCentroids(extdata,pred[[1]],d.isPretreated=F,dist.meth="dqda")$pred)
             
             pred <- cit.centroids(traindata,trainannot[,"OSS"],dist.meth="dqda",rowCentering =NA)
             dqda.noRowCentering <- 1+as.numeric(cit.distToCentroids(extdata,pred[[1]],d.isPretreated=F,dist.meth="dqda")$pred)
       
             res <- as.data.frame(cbind(pearson.medianRowCentering,pearson.noRowCentering,dqda.medianRowCentering,
                                        cox.continuous=cox.continuous[,1],cox.discrete=coxtmp,dqda.noRowCentering))
             rownames(res) <- names(extdata)
      
             res <- res[,c(4,5,1:3,6)]
      }else{
             res <- as.data.frame(cbind(cox.continuous=cox.continuous[,1],cox.discrete=coxtmp))                                        
             rownames(res) <- names(extdata)
      }
      PRED <-   as.data.frame(LPRED)
      if(more){
        PRED <- PRED[,1:10]
        w <- grep("var",names(PRED))
        names(PRED)[w] <- gsub(".medCen","",names(PRED)[w])   
        PRED <- PRED[,c(1:4,9:10,5:6)]
      }
      list("results"=res,"predictors"=PRED)
      
}
