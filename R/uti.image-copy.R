# Rd
# description >> graphical function, to plot a dataframe in a heatmap style, but with distinct color codes for numerical and nominal variables
# argument
# item >> df >> dataframe, containing the values to be plotted   
# item >> orderBy >> vector, containing ordered columns index/names, used to sort the rows of df. If NULL, the original order will be used. (AIM to modify the order of the rows)
# item >> colList >> list, containing one item per variable (ie column of df), giving the colors used to represent the values of this variable. NB: for any nominal variable, one can specify the correspondence between colors and modalities, simply by naming the vector of colors with the modalities
# item >> columnsOrder >> vector, containing ordered columns index/names, used to modify the columns order in the plot. If NULL, the original order will be used.  (AIM to modify the order of the columns)
# item >> transpose >> boolean, if TRUE the plotted matrix is transposed (rows become columns and vice-versa)
# item >> refColumnForPvalues >> index/name of column, corresponding to a NOMINAL variable, and used to perform tests with all other variables (chi2 tests with other nominal variables, and anova with numerical variables) (optional)
# item >> strat >>  ...
# item >> line.ind >>  ...
# item >> labColPos >> character string, indicating the position of the columns labels (names) in the plot (above, below, NA). If NA the columns  labels won't be plotted.
# item >> labRowPos >> character string, indicating the position of the row labels (names) in the plot (above, below, NA). If NA the row  labels won't be plotted.
# item >> pvalPos >>  character string, indicating the position of the pvalues calculated from the association tests between all the variables and  the variable refColumnForPvalues 
# item >> col.pval >> character string, color used to plot the pvalues
# item >> colna >>  character string, color used for NA values
# item >> colpos >> character string, color used for positive values (NB: this parameter is used only if colList is NULL)
# item >> colneg >>  character string, color used for negative values (NB: this parameter is used only if colList is NULL)
# item >> colfun >>  function, to be used to generate colors (NB: this parameter is used only if colList is NULL)
# item >> palette >> character vector, three colors used to get a gradient of colors using rampPalette function (for numerical variables) (NB: this parameter is used only if colList is NULL)
# value >> a list of 3 items, item 1 is a list containing one item per variable (ie column of df) giving the colors used to represent the values of this variable, item 2 is the matrix of colors derived from df, item 3 is the ordered vector of distinct colors used in item 2 
# author >> Aurelien de Reynies
# keyword >> methods
# end  
cit.image <- function(df,orderBy=NULL,colList=NULL,columnsOrder=NULL,transpose=FALSE,refColumnForPvalues=NA,strat=NULL,line.ind=NULL,labColPos=c("above","below",NA)[1],labRowPos=c("left","right",NA)[1],pvalPos=c("left","right")[2],col.pval="red",colna = "grey", colpos =  "black", colneg =  "white", colfun = cit.rainbow,palette=c( "royalblue", "white","orangered"),...) {
  df <- as.data.frame(df)
  
  if(!is.null(orderBy)){
        o <- cit.orderDf(df, orderBy)
        df <- df[o,]
  }
  
  if(!is.null(columnsOrder)){
        df <- df[,columnsOrder]
        if(!is.null(colList)) colList <- colList[columnsOrder]
  }
  
  tmp <- cit.dfToColor (df,colList=colList, colna = colna, colpos = colpos, colneg =colneg, colfun = colfun,palette=palette)  
 
  matcol <- tmp[[1]]
  colList <- tmp[[2]]

  if(is.null(columnsOrder)) columnsOrder <- colnames(matcol)
  
  pvalues <- NULL
  fun <- function(dfcol){
       colors <- names(table(unlist(dfcol)))
       res <- as.data.frame(matrix(as.numeric(factor(unlist(dfcol),labels=colors)),ncol=ncol(dfcol)))
       dimnames(res) <- dimnames(dfcol)
       list(res ,colors)
  }
  
  temp <- fun(matcol)
  csc <- as.matrix(temp[[1]][,columnsOrder])
  csc.colors <- temp[[2]]
  

  if(transpose){
      image(t(csc[,columnsOrder]), col = as.vector(csc.colors), axes = FALSE, ...)
  }else{
      image(csc[,columnsOrder], col = as.vector(csc.colors), axes = FALSE, ...)
  }
  
  if(!is.na(labColPos)){axis(c("below"=1,"left"=2,"above"=3,"right"=4)[labColPos],at=seq(0,1,1/(nrow(matcol)-1)),labels=dimnames(matcol)[[1]],las=2,cex.axis=.7)}
  
  if(!is.na(labRowPos) & length(colnames(matcol)) > 0){
    axis(c("below"=1,"left"=2,"above"=3,"right"=4)[labRowPos], 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),columnsOrder,las = 2,tick = FALSE)
  }
  
  if(!is.na(refColumnForPvalues)){
      axis(c("below"=1,"left"=2,"above"=3,"right"=4)[pvalPos], 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),pvalues[columnsOrder] ,las=2,cex.axis=.7,col.axis=col.pval,tick=F)
  }
  
  box()
  list(colList=colList,colMatrix=csc[,columnsOrder],col=as.vector(csc.colors))
}



# Rd
# description >> ordering a dataframe according to one or several of its columns
# argument
# item >> df >> dataframe to be ordered
# item >> orderCols >>  vector, containing ordered columns index/names, used to sort the rows of df. (AIM to modify the order of the rows)
# value >> a numerical vector  , containing the order index of the  rows of df
# author >> Aurelien de Reynies
# keyword >> methods
# end  
cit.orderDf <- function(df, orderCols){
  o <- NULL
  orderCols <- getNums(orderCols,dimnames(df)[[2]])
  expr <- ("o <- order(")
  for(i in orderCols) expr <- paste(expr, "df[,", i, "],", sep = "")
  k <- nchar(expr)
  expr <- paste(substr(expr,1,k-1),")")
  eval(parse(text = expr))
  o
}


# Rd
# description >> getting a dataframe of colors from a dataframe of numerical values
# argument
# item >> numDf >> dataframe, with numerical variables
# item >> colList >> list, containing one item per variable (ie column of numDf), giving the colors used to represent the values of this variable   (optional)
# item >> palette >> ...
# item >> N >> ...
# item >> colna >> ...
# value >> list
# author >> Aurelien de Reynies
# keyword >> methods
# end  
cit.numDfToColor <- function(numDf,colList=NULL,palette=c("royalblue", "white", "orangered"),N=50,colna = "grey"){
    if(is.null(colList)){
         L <- lapply(names(numDf),function(k)cit.numValToColor(numDf[,k],palette=palette,N=N,colna = colna))
    }else{
         L <- lapply(names(numDf),function(k)cit.numValToColor(numDf[,k],col=colList[[k]],palette=palette,N=N,colna = colna))
    }          
    colList <- lapply(L,function(z)z[[2]])
    names(colList) <- names(numDf)
    dfCol <- as.data.frame( lapply(L,function(z)z[[1]]) )
    dimnames(dfCol)  <- dimnames(numDf)    

    list(dfCol=dfCol,colList=colList) 
}

# Rd
# description >> getting a vector of colors from a vector of numerical values
# argument
# item >> val >> vector of numerical values
# item >> col >> colors used to represent the numerical values 
# item >> palette >> ...
# item >> N >> ...
# item >> colna >> ...
# value >> list of 2 items, item 1 = vector of colors derived from the vector of numbers, item 2 = ordered vector of distincts colors used in item 1
# author >> Aurelien de Reynies
# keyword >> methods
# end  
cit.numValToColor <- function(val,col=NULL,palette=c("royalblue", "white", "orangered"),N=50,colna = "grey"){
  if(is.null(col)) col <- colorRampPalette(palette)(N)
  col <- unique(col)
  N <- length(col)
  v <- sort(setdiff(val,NA))
  m <- v[1]
  M <- v[length(v)]
  delta <- (M-m)/(N-1)
  lim <- seq(m+(delta/2),M-(delta/2),by=delta)
  valDiscr <- 1+sapply(val,function(z)sum(z>lim))
  valToCol <- col[valDiscr]
  valToCol[which(is.na(valDiscr))] <- colna
  list(valToCol=valToCol,col=col)  
}


# Rd
# description >> getting a dataframe of colors from a dataframe of nominal values
# argument
# item >> nomDf >> dataframe, with nominal variables
# item >> colList >> list, containing one item per variable (ie column of numDf), giving the colors used to represent the values of this variable   (optional)
# item >> colna >> ...
# item >> colpos >> ...
# item >> colneg >> ...
# item >> colfun >> ...
# value >> list of 2 items, item 1 = dataframe of colors derived from the dataframe of nominal variables, item 2 = list of (per variable) color codes
# author >> Aurelien de Reynies
# keyword >> methods
# end 
cit.nomDfToColor<- function(nomDf,colList=NULL,colna = "grey", colpos = "black", colneg = "white",colfun = cit.rainbow){
  if(is.null(colList)){
         L <- lapply(names(nomDf),function(k)cit.nomValToColor(nomDf[,k],colna =colna, colpos=colpos,colneg =colneg,colfun=colfun))
    }else{
         L <- lapply(names(nomDf),function(k)cit.nomValToColor(nomDf[,k],col=colList[[k]],colna =colna, colpos=colpos,colneg =colneg,colfun=colfun))
    }          
    colList <- lapply(L,function(z)z[[2]])
    names(colList) <- names(nomDf)
    dfCol <- as.data.frame( lapply(L,function(z)z[[1]]) )
    dimnames(dfCol)  <- dimnames(nomDf)    

    list(dfCol=dfCol,colList=colList) 
}


# Rd
# description >> getting a vector of colors from a vector of nominal values
# argument
# item >> val >> vector of nominal values
# item >> col >> vector of colors used to represent the nominal values 
# item >> colna >> ...
# item >> colpos >> ...
# item >> colneg >> ...
# item >> colfun >> ...
# value >> list of 2 items, item 1 = vector of colors derived from the vector of nominal variables, item 2 = ordered vector of distinct colors used in item 1
# author >> Aurelien de Reynies
# keyword >> methods
# end 
cit.nomValToColor <- function(val,col=NULL,colna = "grey", colpos = "black", colneg = "white",colfun = cit.rainbow){ 
  tval <- table(val)
  vals <- names(tval)
  n <- length(tval)
  if(is.null(col)){
     if(n==2){
        col <- c(colneg,colpos)
        neg <- intersect(vals,c("WT","Wt","wt","mss","MSS","-","0","N","n","NEG","Neg","No","NO","NON","Non","non"))
        if(length(neg)>0){
           names(col) <- c(neg,setdiff(vals,neg))
        }        
     }else{
        col <- colfun(n) 
     }
  }
  
  if(!is.null(names(col))) {
       valToCol <- col[as.character(val)]      
  }else{
       valToCol <- col[as.numeric(as.factor(val))]
  }
  valToCol[which(is.na(valToCol))] <- colna
  list(valToCol=valToCol,col=col)  
}

    
    


# Rd
# description >>  getting a dataframe of colors from a dataframe
# argument
# item >> df >> dataframe
# item >> colList >> list, containing one item per variable (ie column of df), giving the colors used to represent the values of this variable   (optional)
# item >> colna >> ...
# item >> colpos >> ...
# item >> colneg >> ...
# item >> colfun >> ...
# item >> palette >> ...
# value >> list of 2 items, item 1 = dataframe of colors derived from the dataframe, item 2 = list of (per variable, ie column of df) color codes
# author >> Aurelien de Reynies
# keyword >> methods
# end 
cit.dfToColor <- function (df, colList = NULL, colna = "grey", colpos = "black", 
    colneg = "white", colfun = cit.rainbow, palette = c("royalblue", 
        "white", "orangered")){ 

    vars <- rep("other",ncol(df))
    names(vars) <- colnames(df)
    
    wNumerical <- which(sapply(1:ncol(df), function(i) is.numeric(df[,i]))) 
    vars[wNumerical] <- "numerical"
    
    wBinary <- which(sapply(1:ncol(df), function(i) length(setdiff(unique(df[,i]),NA))==2)) 
    vars[wBinary] <- "binary"    
    
    wOther <- which(vars=="other")
    wOB <- c(wOther,wBinary)
    
    wNumerical <- setdiff(wNumerical, wOB)

    a <- b <- NULL
    if (length(wOB) > 0) 
        a <- cit.nomDfToColor(df[, wOB, drop = F], colList = colList, 
            colna = colna, colpos = colpos, colneg = colneg, 
            colfun = colfun)
    if (length(wNumerical) > 0) 
        b <- cit.numDfToColor(df[, wNumerical, drop = F], 
            colList = colList, palette = palette, colna = "grey")
    if (is.null(a)) {
        m <- b[[1]]
        colList <- b[[2]]
    }    else {
        if (is.null(b)) {
            m <- a[[1]]
            colList <- a[[2]]
        }
        else {
            m <- as.data.frame(cbind(a[[1]], b[[1]]))
            names(m) <- colnames(df)[c(wOB, wNumerical)]
            colList <- c(a[[2]], b[[2]])[names(m)]
        }
    }
    for (i in names(m)) m[, i] <- as.character(m[, i])
   
    list(dfToCol = m[,colnames(df)], colList = colList[colnames(df)])
}

