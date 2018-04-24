DIST.METHODS <- c( "euclidian", "maximum" ,  "manhattan", "canberra", "binary",
                  "minkowski", "pearson" ,  "spearman", "cosine")
COR.USE <- c( "pairwise.complete.obs" ,"all.obs"    , "complete.obs")



# Rd
# description >> as \link[base]{load} but allow to change the name of the object to load. 
# argument
# item >> filename >> the name of the \file{.RData} to load
# value >> required
# author >> A de Reynies, M Guedj
# keyword >> utilities
# seealso >> \link[cit.utils]{cit.data}, \link[base]{load}
# end
cit.load <- function(filename){  
# description : comme load pour pouvoir utiliser l'objet sans conna?tre son nom
    if(file.exists(filename)) return(eval(parse(text=load(filename))))
    cat(paste("error - function cit.load : file ",filename," doesn't exist\n"))
    NULL
}

# Rd
# description >> internal
# argument
# item >> obj >> ...
# item >> fun >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
is.ok <- function (obj, fun = all) {
    if (is.null(obj)) return(FALSE)

    if (length(obj) == 0)  return(FALSE)

    fun(!is.na(obj) & !is.nan(obj) & !is.infinite(obj))
}

# Rd
# description >> internal
# argument
# item >> x >> ...
# item >> nbdec >> ...
# item >> nmax >> ...
# item >> pow >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.sciFormat  <-function (x, nbdec = 3, nmax = 4, pow = c("e", "10")) {

    pow = match.arg(pow)
    w <- which(!is.na(x))
    if (length(w) == 0)
        return(x)
    L = strsplit(as.character(format(signif(x[w], nbdec), scientific = TRUE)),
        split = "e")
    tmp <- NULL
    for (i in 1:length(L)) {
        if (abs(as.numeric(L[[i]][2])) < nmax) {
            tmp <- c(tmp, format(signif(x[w][i], nbdec), scientific = FALSE))
        }
        else {
            tmp <- c(tmp, format(signif(x[w][i], nbdec), scientific = TRUE))
        }
    }
    if (pow == "10")
        tmp <- gsub("e", " 10", tmp)
    x[w] <- tmp
    x
}

# Rd
# description >> internal
# argument
# item >> x >> ...
# item >> lim >> ...
# item >> quant >> ...
# item >> addlevels >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.discretize <- function (x, lim, quant = FALSE, addlevels = FALSE) {

    lim <- sort(lim)
    if (quant & (any(lim > 1) | any(lim < 0)))
        stop("lim must be [0;1] as quant=T\n")
    res <- rep(NA, length(x))
    if (quant)
        lim <- quantile(x, probs = lim, na.rm = TRUE)
    n <- length(lim)
    for (i in n:1) res[which(x < lim[i])] <- i
    res[which(x >= lim[n])] <- n + 1
    if (addlevels) {
        res <- as.factor(res)
        if (quant)
            lim <- gsub(" ", "", prettyNum(lim, format = "g",
                digits = 1))
        if (length(lim) == 1)
            lev <- c(paste("<", lim[1], sep = ""), paste(">=",
                lim[1], sep = ""))
        else lev <- c(paste("<", lim[1], sep = ""), paste("[",
            lim[-length(lim)], ";", lim[-1], "[", sep = ""),
            paste(">=", lim[length(lim)], sep = ""))
        levels(res) <- lev[as.numeric(levels(res))]
    }
    res
}

# Rd
# description >> internal
# argument
# item >> data >> ...
# item >> partition >>  ...
# item >> MARGIN >> ...
# item >> fAggreg >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.dfAggregate <- function (data, partition, MARGIN = 1, fAggreg = mean.na) {
    cMARGIN <- setdiff(c(1, 2), MARGIN)
    n <- length(partition)
    N <- dim(data)[MARGIN]
    p <- dim(data)[cMARGIN]
    if (n != N)
        stop("ERROR - cit.dfAggregate : size of partition doesn't correspond to data dimension")
    l <- split(1:N, partition)
    d <- data
    if (MARGIN == 2)d <- t(data)

    d <- matrix(sapply(l, function(i) if (length(i) == 1) {
                                         unlist(d[i, ])
                                      }else {
                                         apply(d[i, ], 2, fAggreg)
                                      }), ncol = p, byrow = TRUE)



    d <- as.data.frame(d)
    rownames(d) <- names(l)
    names(d) <- dimnames(data)[[cMARGIN]]
    if (MARGIN == 2)
        d <- as.data.frame(t(d))
    d
}

# Rd
# description >> internal
# argument
# item >> x >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
mean.na <- function (x){
    mean(x, na.rm = TRUE)
}

# Rd
# description >> internal
# argument
# item >> x >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
median.na <-function (x){ 
    median(x, na.rm = TRUE)
}


# Rd
# description >> internal
# argument
# item >> x >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
var.na <- function (x) {
    var(x, na.rm = TRUE)
}


# Rd
# description >> internal
# argument
# item >> x >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
sum.na <- function (x) {
    sum(x, na.rm = TRUE)
}


# Rd
# description >> internal
# argument
# item >> x >> ...
# item >> meth >> ...
# item >> use >> ...
# item >> diag >> ...
# item >> upper >> ...
# item >> p >>  ...
# item >> replaceNA >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.dist <- function (x, meth = DIST.METHODS[7], use = COR.USE[1], diag = FALSE,
                      
                      upper = FALSE, p = 2, replaceNA = TRUE)
  
{
  
  res <- NULL
  
  if (is.numeric(meth))
    
    meth <- DIST.METHODS[meth]
  
  if (is.numeric(use))
    
    use <- COR.USE[use]
  
  if (meth == "euclidean") {
    
    meth = "euclidian"
    
  }
  
  if (!is.na(meth) & !is.null(meth) & meth %in% DIST.METHODS) {
    
    if (meth == "propNonEqual") {
      
      nc <- ncol(x)
      
      nr <- nrow(x)
      
      m <- apply(x, 1, function(xi) rowSums(t(t(x) != xi),
                                            
                                            na.rm = T) + rowSums(is.na(t(t(x) != xi))))/nc
      
      res <- as.dist(m, diag = diag, upper = upper)
      
      maxdist <- ceiling(max(res, na.rm = TRUE))
      
    }
    
    if (meth == "cosine") {
      
      tmp <- x %*% t(x)
      
      norm <- sqrt(diag(tmp)) %*% t(sqrt(diag(tmp)))
      
      res <- as.dist(1 - (tmp/norm), diag = diag, upper = upper)
      
      maxdist <- ceiling(max(res, na.rm = TRUE))
      
    }
    
    if (meth %in% c("euclidian", "maximum", "manhattan",
                    
                    "canberra", "binary", "minkowski")) {
      
      res <- dist(x, method = meth, diag = diag, upper = upper,
                  
                  p = p)
      
      maxdist <- ceiling(max(res, na.rm = TRUE))
      
    }
    
    if (meth %in% c("pearson", "spearman")) {
      
      maxdist <- 2
      
      if (!is.ok(use)) {
        
        use <- "all.obs"
        
        
        
      }
      
      res <- as.dist(1 - cor(t(x), use = use, method = meth),
                     
                     diag = diag, upper = upper)
      
    }
    
  }
  
  wNA <- which(is.na(res))
  
  if (length(wNA) > 0) {
    
    if (replaceNA) {
      
      res[wNA] <- maxdist
      
      
    }
    
  }
  
  if (is.null(res))
    
    stop("ERROR - function cit.dist : uncorrect value for parameter 'meth' and/or 'use'")
  
  return(res)
  
}

# Rd
# description >> internal
# argument
# item >> d >> ...
# item >> margin >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.rowMeansVars <- function (d, margin = 1){

    if (margin == 2)
        d <- t(d)
    m <- rowMeans(d, na.rm = T)
    dif <- d - m
    ssd <- rowSums(dif^2, na.rm = T)
    list(means = m, sumsquaredif = ssd, vars = ssd/(ncol(d) - 1), `centered rows` = dif)

}



# Rd
# description >> internal
# argument
# item >> v >> ...
# item >> index >> ...
# item >> bipartition >> ...
# item >> nomModalite >> ...
# item >> clasNA >> ...
# value >> ...
# author >> A de Reynies, M Guedj
# keyword >> internal
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
cit.recode <- function(v=NULL, index=NULL, bipartition=TRUE,nomModalite=TRUE,clasNA=NULL){
	if(is.null(index))stop("error - recode : index=NULL")
	if(is.null(v))stop("error - recode :v=NULL")

	v.num <- as.numeric(as.factor(v))
	modalites <- levels(as.factor(v))
	p<-length(modalites)
	recod <- as.data.frame(rep(NA,p))
	if(bipartition){
		recod[index,1] <- paste(index,collapse="-")
		recod[-index,1] <- paste(c(1:p)[-index],collapse="-")
	}else{
		n <- length(levels(as.factor(index)))
		for(i in 1:n){
			x<-which(index==i)
			recod[x,1] <- paste(x,collapse="-")
		}
		if(!is.null(clasNA)) recod[which(index==clasNA),1]<-NA
	}
	if(nomModalite){
		for(i in 1:dim(recod)[1]){
			if(!is.na(recod[i,1]))
				recod[i,1]<- paste(modalites[as.numeric(unlist(strsplit(recod[i,1],"-")))],collapse="-") # modified by MG (01.12.08)
        # recod[i,1]<- paste(modalites[getNumeric(recod[i,1])],collapse="-")
		}
	}
	recod[v.num,1]
}

# Rd
# description >> internal
# argument
# item >> u >> ...
# item >> ascol >> ...
# value >> ...
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
cit.specialDim <- function(u,ascol=TRUE){
  o <- c(1,2)
  if(!ascol) o <- c(2,1)
	r <- dim(u)
	if(is.null(r)) r<-c(length(u),1)[o]
	r
}

# Rd
# description >> internal
# argument
# item >> L >> ...
# item >> x >> ...
# value >> ...
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
addToList <- function(L, x){
 c(L,list(x))
}

# Rd
# description >> internal
# argument
# item >> file >> a file name
# item >> width >> the width of the pdf
# item >> height >> the height of the pdf
# value >> no value
# author >> E Thomas, M Guedj
# keyword >> internal
# end
pdf2 <- function(file=NULL,width=16.53,height=11.69){
    pdf(file=file,width=width,height=height)
}

# Rd
# description >> internal
# argument
# item >> co >>  value  to be checked
# item >> nam >>  vector of names of the dataframe
# value >> a character string  
# author >> A de Reynies
# keyword >> internal
# end
getName <- function(co,nam){
    if(is.null(co)) return(NULL)
    if(is.numeric(co)){
         setEmptyToNULL(nam[intersect(co,1:length(nam))])
    }else{
         setEmptyToNULL(intersect(co,nam))
    }
}

# Rd
# description >> internal
# argument
# item >> vco >>  vector  to be checked
# item >> nam >>  vector of names of the dataframe
# item >> silent >> ...
# value >> a character vector  
# author >> A de Reynies
# keyword >> internal
# end
getNames <- function(vco,nam,silent=TRUE){
    if(is.numeric(vco))return(nam[vco])
    temp <- NULL
    for(co in vco) temp <- c(temp,getName(co,nam))
    if(!silent & length(temp)<length(vco))cat("Warning - funtion getNames : some elements are not in the reference set\n")
    temp
}

# Rd
# description >> internal  
# argument
# item >> x >>  a possibly empty vector
# value >> a vector or NULL  
# author >> A de Reynies
# keyword >> internal
# end
setEmptyToNULL <- function(x){
    if(length(x)>0) return(x)
    NULL
}

# Rd
# description >> internal
# argument
# item >> p1 >> required
# item >> p2 >> required
# value >> 
# author >> A de Reynies
# keyword >> internal
# end
cit.disym <- function(p1, p2){      
# description : distance de la difference symetrique entre les partitions p1 et p2 normalisee ? un

	n <- length(p1) 
	nij <- table(p1, p2)
	ni <- apply(nij, 1, sum)
	nj <- apply(nij, 2, sum)
	(sum(ni * ni) + sum(nj * nj) - 2 * sum(nij * nij))/(n*(n-1))
}

# Rd
# description >> internal
# argument
# item >> z >> required
# value >> 
# author >> A de Reynies
# keyword >> internal
# end
nb.na <- function(z){
   sum(is.na(z))
}

# Rd
# description >> internal
# argument
# item >> n >> required
# item >> k >> required
# value >> 
# author >> A de Reynies
# keyword >> internal
# end
Cnk  <- function (n, k) {
    if (n == k) 
        return(1)
    a <- 1
    b <- 1
    c <- n
    d <- 1
    for (i in 1:min(n - k, k)) {
        a <- a * i
        b <- b * c
        c <- c - 1
    }
    b/a
}



# Rd
# description >> internal
# argument
# item >> u >> a first vector of any data type (numeric, character, factor...)
# item >> v >> a second vector of any data type (numeric, character, factor...)
# value >> a \code{(length(u)*length(v), 2)} data.frame with each row corresponding to one combination of the elements of \code{u} and \code{v} 
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
cit.combine <- function(u,v){
	# code modified by MG (01.08.08) in order to be simpler and to take any type of data in u and v
  U = NULL
  for (iii in 1:length(u))
    U = c(U, rep(u[iii], length(v)))
  V = rep(v, length(u))
  return(data.frame(u = U, v = V))
}


# Rd
# description >> internal
# argument
# item >> u >> required
# item >> v >> required
# item >> colla >> required
# item >> inv >> required
# value >> required
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
cit.paste <- function(u,v,colla=",",inv=FALSE){ 
# description :  N x M paste
   temp <- cit.combine(u,v)
   if(inv) temp <- temp[,c(2,1)]
   apply(temp,1,function(x)paste(x,collapse=colla))
}

# Rd
# description >> internal
# argument
# item >> n >> the number of colors (>= 1) to be in the palette
# item >> default2Col >> default colors if \code{n = 2}. If \code{n = 1}, the function returns \code{default2Col[1]}
# item >> default6Col >>  default colors if \code{n > 2 & n <= 6}
# item >> colorMode >> color mode: \code{"rainbow"} by default or \code{"gray"/"grey"} to obtain a gradient of gray
# value >>  A character vector of colors of length \code{n}
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
cit.rainbow <- function(n, default2Col = c("black","red"), default6Col = c("red", "blue", "green", "orange", "purple", "gray" ), colorMode = c("rainbow", "gray", "grey")){
  colorMode = match.arg(colorMode)
  if (colorMode == "grey")
    colorMode = "gray"
  if (n == 1)
    return(default2Col[1])
  if (n == 2)
    return(default2Col )
  if (n > 6 & colorMode == "rainbow")
    return(rainbow(n))
  if (n < 7 & n > 2 & colorMode == "rainbow")
    return (default6Col[1:n])
  if (n > 2 & colorMode == "gray")
    return(gray (seq (1,0, length = n))) #modif JL 140410, previous : seq(0.8,0.2,...)
}



# Rd
# description >> internal
# argument
# item >> annot >> required
# item >> default2Col >> default colors if \code{n = 2} in a \link[cit.utils]{cit.rainbow} call
# item >> colorMode >> the color mode (\code{rainbow}, \code{gray} or \code{grey}) in a \link[cit.utils]{cit.rainbow} call
# item >> defaultBgTxtCol >> required
# item >> ignor.col >> required
# value >> required
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
cit.createListColannot <- function(annot,default2Col=c("black","white"),colorMode= c("rainbow","gray")[1],defaultBgTxtCol=c("white","black"),ignor.col=FALSE){
        L <- NULL
        if(is.null(dim(annot))){
            n <- length(setdiff(unique(annot),c(NA)))
            res <- cit.rainbow(n,default2Col=default2Col,colorMode=colorMode) # modified by MG (01.08.08)
            res2 <- NULL
            for(i in 1:n) res2 <- c(res2,cit.getColors(res,i,defaultBgTxtCol=defaultBgTxtCol,ignor.col=ignor.col,colorMode=colorMode))   # modified by MG (01.08.08)
            L <- list(res2)
        }else{
            for(j in 1:(dim(annot)[2])){
                 n <- length(setdiff(unique(annot[,j]),c(NA)))
                 res <- cit.rainbow(n,default2Col=default2Col,colorMode=colorMode) # modified by MG (01.08.08)
                 res2 <- NULL
                 for(i in 1:n) res2 <- c(res2,cit.getColors(res,i,defaultBgTxtCol=defaultBgTxtCol,ignor.col=ignor.col,colorMode=colorMode))   # modified by MG (01.08.08)
                 L <- addToList(L,res2)
            }
        }
        L
}

# Rd
# description >> internal
# argument
# item >> colannot >> required
# item >> indexModal >> required
# item >> defaultBgTxtCol >> required
# item >> ignor.col >> required
# value >> required
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
cit.getColors <- function ( colannot,
                            indexModal,
                            defaultBgTxtCol = c("white","black"),
                            ignor.col = FALSE,
                            colorMode=c("rainbow","gray")[1]) {
# description : return a vector (background.color, text.color)
    n <- length(colannot)
    if(n==1)return( defaultBgTxtCol)
    
    bgcol  <- defaultBgTxtCol[1]
    texcol <- defaultBgTxtCol[2]

    if(colorMode=="gray"){
       if(n==2){
            bgcol <- colannot[ifelse(indexModal > n/2,1,n)]
            texcol <- colannot[indexModal]
        }
        if (n > 2 & !ignor.col) {
            bgcol <- colannot[indexModal]
            texcol <- colannot[ifelse(indexModal > n/2,1,n)]
        }
    }else{
        if (n >= 2 & !ignor.col) {
            if (n == 2) {
                bgcol <- colannot[ifelse(indexModal > n/2,1,n)]
                texcol <- colannot[indexModal]
            }
            else {
                bgcol <- colannot[indexModal]
                
            }
        }
    }
    return(c(bgcol, texcol))
}


# Rd
# description >> internal
# argument
# item >> d >> a data.frame
# item >> ncmax >> number of character 'max' for an integer 
# value >> the resulting data.frame
# author >> A de Reynies, M Guedj
# keyword >> internal
# end
factoall <- function(d,ncmax = 10){
    n <- ncol(d)
    for(i in 1:n){
        if(is.factor(d[,i])){
            d[,i] <- as.character(d[,i])
            na <- which(is.na(d[,i]))
            num <- suppressWarnings(as.numeric(d[,i]))
            nanum <- which(is.na(num))
            if(length(nanum) == length(na)){
                 int <- suppressWarnings(as.integer(d[,i]))
                 naint <- which(is.na(int))
                 nc <- nchar(num)
                 if(length(naint) == length(nanum) & all(nc < ncmax)) {
                      d[,i] <- int
                 }else{
                      d[,i] <- num
                 }
            }
        }
    }
    d
}

# Rd
# description >> identify rows in \code{annot} verifying the condition \code{variable} \code{op} \code{values} and return corresponding rownames  
# argument
# item >> annot >> a \code{data.frame}
# item >> variable >> a column in \code{annot}
# item >> values >> a vector of values that will be compared to \code{annot[,variable]} 
# item >> op >> an operator used to compare \code{annot[,variable]} and \code{values} , NB : g is for grep
# item >> colId >> can be \code{NULL}, \code{NA} or a column of \code{annot}  
# item >> ignore.case >> for grep operator only
# item >> uniqueVal >> boolean : apply \code{unique()} or not to the result
# value >> a vector : rownames or index or values in \code{annot[,colId]} for rows verifying the given condition
# author >> A de Reynies
# details >>  if \code{colId = NULL} returns rownames, if \code{colId = NA} returns index , else returns corresponding values in \code{annot[,colId]}
# keyword >> utilities
# end
 cit.which  <- function (annot, variable, values, op = c("==", "!=", "in", "!in",
    "<", "<=", ">", ">=", "g","highest","lowest")[1], colId = NULL, ignore.case = FALSE,
    uniqueVal = TRUE)
{
    if (length(values) == 0)
        stop("ERROR - function cit.which  : values empty or null")
    OP <- NA
    try(OP <- match.arg(op, c("in", "!in", "==", "!=", "<", "<=",
        ">", ">=", "g","highest","lowest")), silent = TRUE)
    if (is.na(OP))
        stop("ERROR - function cit.which  : op undefined")

    

    
    if(op %in% c("==","!=","in", "!in")){
        wna <- NULL
        if (any(is.na(values))) wna <- which(is.na(annot[, variable]))

        w <- c(wna, which(annot[, variable] %in% setdiff(values,NA)))
        
        if(op %in% c("!=","!in")) w <- setdiff(1:nrow(annot),w)

    }else{
        if (OP == "g")   w <- mgrep(values, annot[, variable],ignore.case = ignore.case)
        if (OP == "<")   w <- which(annot[, variable] < values[1])
        if (OP == "<=")  w <- which(annot[, variable] <= values[1])
        if (OP == ">")   w <- which(annot[, variable] > values[1])
        if (OP == ">=")  w <- which(annot[, variable] >= values[1])
        if (OP == "highest")  w  <-  order(annot[, variable],decreasing=TRUE)[1:values[1]]
        if (OP == "lowest")  w  <-  order(annot[, variable])[1:values[1]]
    }
    if (length(w) == 0) return(NULL)

    if (is.null(colId))
        return(rownames(annot)[w])
    if (is.na(colId))
        return(w)
    else if (uniqueVal)
        return(unique(annot[w, colId]))
    else return(annot[w, colId])
}

# Rd
# description >> internal
# argument
# item >> patterns >> required
# item >> x >> required
# item >> exact >> required
# item >> aslist >> required
# item >> ignore.case >> required
# item >> value >> required
# value >> required
# author >> A de Reynies
# keyword >> internal
# end
mgrep <- function(patterns,x,exact=FALSE,aslist=FALSE,ignore.case=FALSE,value=FALSE){
    adEl <- function(x,y){c(x,y)}
    z<-NULL
    if(aslist){
       adEl <- function(x,y){c(x,list(y))}
       z <- list()
    }

    if(exact) for(p in patterns) z<-adEl(z,which(x==p))
    else for(p in patterns) z<-adEl(z,grep(p,x,ignore.case=ignore.case))
    
    if(aslist) names(z) <- patterns
    
    index <- unique(z[!is.na(z)])
    if(value) return(x[index])
    index
}

# Rd
# description >> internal 
# argument
# item >> v >> ...
# item >> newmin >>  ...
# item >> newmax >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.changeRange <-function (v, newmin = 1, newmax = 10) 
{
    oldmin <- min(v, na.rm = TRUE)
    oldmax <- max(v, na.rm = TRUE)
    newmin + ((newmax - newmin) * (v - oldmin)/(oldmax - oldmin))
}

# Rd
# description >> internal 
# argument
# item >> d >> ...
# item >> cols >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
suppress.NA <- function (d, cols = 1:2) 
{
    w <- NULL
    for (i in cols) w <- unique(c(w, which(is.na(d[, i]))))
    if (!is.null(w)) 
        if (length(w) > 0) 
            d <- d[-w, ]
    d
}

# Rd
# description >> internal 
# argument
# item >> data >>  ...
# item >> annot >> ...
# item >> var.annot >> ...   
# item >> var.data >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
get.datannot <- function (data = NULL, annot = NULL, var.annot = NULL, var.data = NULL) 
{
    if (is.null(data) | is.null(annot)) 
        stop("ERROR - function get.datannot : missing parameter data or annot")
    if (is.null(var.data)) {
        var.data <- 1:nrow(data)
    }
    else {
        var.data <- getNums(var.data, rownames(data))
    }
    vd <- rownames(data)[var.data]
    if (is.null(var.annot)) {
        var.annot <- 1:ncol(annot)
    }
    else {
        var.annot <- getNums(var.annot, names(annot))
    }
    va <- names(annot)[var.annot]
    S <- unique(c(names(data), rownames(annot)))
    datannot <- as.data.frame(cbind(annot[S, var.annot], t(data[var.data, 
        S])))
    rownames(datannot) <- S
    names(datannot) <- c(va, vd)
    datannot
}

# Rd
# description >> internal 
# argument
# item >> vco >> ...
# item >> nam >>  ...
# item >> silent >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
getNums  <- function (vco, nam, silent = TRUE) 
{
    temp <- NULL
    for (co in vco) temp <- c(temp, getNum(co, nam))
    if (!silent & length(temp) < length(vco)) 
        cat("Warning - funtion getNums : some elements are not in the reference set\n")
    temp
}

# Rd
# description >> internal 
# argument
# item >> co >> ...
# item >> nam >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
getNum  <-function (co, nam) 
{
    if (is.null(co)) 
        return(NULL)
    if (is.numeric(co)) {
        setEmptyToNULL(intersect(co, 1:length(nam)))
    }
    else {
        setEmptyToNULL(which(nam == co))
    }
}

# Rd
# description >> internal 
# argument
# item >> x >> ...
# item >> percentHighestPeak >>  ...
# item >> maxNbPeaks >>  ...
# item >> minDeltaBetweenPeaks >>  ...
# item >> deltaApproach >>  ...
# item >> doplot >>  ...
# item >> bw >>  ...
# item >> ... >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.peaks <-function (x, percentHighestPeak = 0.2, maxNbPeaks = NULL, minDeltaBetweenPeaks = 0.03, 
                      deltaApproach = 1, doplot = FALSE, bw = "nrd0", ...) 
{
  if (is.na(minDeltaBetweenPeaks)) 
    minDeltaBetweenPeaks <- NULL
  v <- cit.density(x, pc = percentHighestPeak, bw = bw)
  m <- v$y[which(v$x %in% v$top)]
  w <- which(m/max(m) > percentHighestPeak)
  xw <- v$top[w]
  which.eq <- function(z) {
    order(abs(v$x - z))[1]
  }
  if (deltaApproach == 1 & !is.null(minDeltaBetweenPeaks) & 
      length(xw) > 1) {
    for (i in 1:(length(xw) - 1)) {
      if (xw[i + 1] - xw[i] < minDeltaBetweenPeaks) {
        xw[i:(i + 1)] <- xw[c(i:(i + 1))][which.max(m[i:(i + 
                                                           1)])]
      }
    }
    xw <- unique(xw)
  }
  if (!is.null(maxNbPeaks)) {
    if (length(xw) > maxNbPeaks) 
      for (i in 1:(length(xw) - maxNbPeaks)) {
        oneamong <- xw[c(0, 1) + which.min(diff(xw))]
        w1 <- which.eq(oneamong[1])
        w2 <- which.eq(oneamong[2])
        out <- oneamong[which.min(v$y[c(w1, w2)])]
        xw <- setdiff(xw, out)
      }
  }
  umin <- NULL
  if (length(xw) > 1) {
    for (i in 1:(length(xw) - 1)) {
      possiblemin <- v$down[which(v$down >= xw[i] & v$down <= 
                                    xw[i + 1])]
      if (length(possiblemin) > 1) {
        w <- sapply(possiblemin, which.eq)
        possiblemin <- possiblemin[which.min(v$y[w])]
      }
      umin <- c(umin, possiblemin)
    }
  }
  L <- list(`x abciss big peaks` = xw, `x abciss inter-peaks` = umin, 
            `nb big peaks` = length(xw))
  peaks <- L[[1]]
  interpeaks <- L[[2]]
  temp <- as.data.frame(t(sapply(peaks, function(pic) {
    wleft <- which(interpeaks < pic)
    if (length(wleft) > 0) {
      wleft <- max(interpeaks[wleft])
    }
    else {
      wleft <- min(x, na.rm = TRUE)
    }
    wright <- which(interpeaks > pic)
    if (length(wright) > 0) {
      wright <- min(interpeaks[wright])
    }
    else {
      wright <- max(x, na.rm = TRUE)
    }
    c(pic, wleft, wright, length(which(x >= wleft & x <= 
                                         wright)))
  })))
  names(temp) <- c("peak", "left born", "right born", "size")
  L$"peaks x size" <- temp
  if (deltaApproach == 2 & !is.null(minDeltaBetweenPeaks) & 
      L$"nb big peaks" > 1) {
    Lini <- L
    names(Lini) <- paste("initial", names(L))
    while (any(diff(L$"x abciss big peaks") < minDeltaBetweenPeaks) & 
           L$"nb big peaks" > 1) {
      w <- which.min(abs(diff(L$"x abciss big peaks")))
      lb <- min(L$"peaks x size"[w:(w + 1), "left born"])
      rb <- max(L$"peaks x size"[w:(w + 1), "right born"])
      si <- sum(L$"peaks x size"[w:(w + 1), "size"])
      pe <- median(x[which(x >= lb & x <= rb)])
      L$"peaks x size"[w, ] <- c(pe, lb, rb, si)
      L$"peaks x size" <- L$"peaks x size"[-(w + 1), ]
      L$"nb big peaks" <- L$"nb big peaks" - 1
      L$"x abciss big peaks"[w] <- pe
      L$"x abciss big peaks" <- L$"x abciss big peaks"[-(w + 
                                                           1)]
      L$"x abciss inter-peaks" <- L$"x abciss inter-peaks"[-(w + 
                                                               1)]
    }
    L <- c(Lini, L)
  }
  if (doplot) {
    plot(v, ...)
    abline(v = L$"x abciss big peaks", col = "red")
    if (length(L$"x abciss inter-peaks") > 0) 
      abline(v = L$"x abciss inter-peaks", col = "green", 
             lty = 3)
  }
  L
}

# Rd
# description >> internal 
# argument
# item >> x >> ...
# item >> doplot >>  ...
# item >> pc >>  ...
# item >> ... >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.density <-function (x, doplot = FALSE, pc = 0.05, ...) 
{
  dx <- density(x, na.rm = TRUE, ...)
  ymax <- diff(range(dx$y))
  n <- length(dx$y)
  croissance <- as.numeric((dx$y[-1] - dx$y[-n]) > 0)
  wB <- 1 + which(diff(croissance) == 1)
  wH <- 1 + which(diff(croissance) == -1)
  pointsBas <- dx$x[wB]
  pointsHauts <- dx$x[intersect(wH, which(dx$y > pc * ymax))]
  if (length(pointsHauts) > 0) 
    pointsBas <- sapply(split(pointsBas, cit.discretize(pointsBas, 
                                                        pointsHauts)), median)
  if (doplot) {
    plot(dx, ...)
    abline(v = pointsBas, lty = 3, col = "blue")
    abline(v = pointsHauts, lty = 3, col = "red")
  }
  L <- c(dx, list(down = pointsBas, top = pointsHauts))
  attr(L, "class") <- "density"
  L
}

# Rd
# description >> internal 
# argument
# item >> f >> ...
# item >> r >>  ...
# item >> sep >>  ...
# item >> na.strings >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
read2 <-function (f, r = 1, sep = "\t", na.strings = "NA") 
{
  read.csv(f, header = TRUE, sep = sep, row.names = r, na.strings = na.strings)
}


# Rd
# description >> internal 
# argument
# item >> d >> ...
# value >> ...
# author >> E Thomas, M Guedj
# keyword >> internal
# end
factorize <- function (d) 
{
  for (i in 1:ncol(d)) d[, i] <- as.factor(d[, i])
  d
}



# Rd
# description >> internal 
# argument
# item >> df >> ...
# item >> reorder >>  ...
# item >> reorder.method >>  ...
# item >> ref >>  ...
# item >> labelscolors >>  ...
# item >> default2cols >>  ...
# item >> defaultcolfun >>  ...
# item >> addrefdelim >>  ...
# item >> addrowdelim >>  ...
# item >> addsamplenames >>  ...
# item >> addPval >>  ...
# item >> addPval.function >>  ...
# item >> cex.annot >>  ...
# item >> lwd >>  ...
# value >> ...
# author >> L Marisa
# keyword >> internal
# end
cit.dfSegmentplot = function (df, reorder = TRUE, reorder.method = c("all", "paired")[1], 
                              ref = 1, labelscolors = NULL, default2cols = c("black", "lightgrey"), 
                              defaultcolfun = cit.rainbow, addrefdelim = FALSE, 
                              addrowdelim = TRUE, addsamplenames = NULL, addPval = FALSE, 
                              addPval.function = chisq.test, cex.annot = 1, lwd = 2, 
                              ...) 
{
  df <- factorize(df)
  df <- as.data.frame(df)
  if (!is.numeric(ref)) 
    ref <- which(names(df) == ref)
  if (addPval) {
    pvs <- NULL
    for (i in 1:ncol(df)) {
      t <- table(df[, c(ref, i)])
      t <- t[apply(t, 1, sum) > 0, ]
      t <- t[, apply(t, 2, sum) > 0]
      p <- NA
      suppressWarnings(try(p <- cit.sciFormat(addPval.function(t)$p.value)))
      pvs <- c(pvs, p)
    }
    names(df) <- paste(names(df), "   p=", pvs, sep = "")
  }
  if (reorder) {
    if (reorder.method == "all") {
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
    }
    else {
      expr <- paste("df. <- data.frame(", sep = "")
      if (length(setdiff(1:(ncol(df) - 1), ref)) > 0) {
        for (c in 1:(ncol(df) - 1)) expr <- paste(expr, 
                                                  "df[order(df[,", ref, "], df[,", c, "]),", 
                                                  c, "], ", sep = "")
        expr <- paste(expr, "df[order(df[,", ref, "], df[,", 
                      c + 1, "]),", c + 1, "])", sep = "")
      }
      else {
        c <- setdiff(1:(ncol(df)), ref)
        expr <- paste(expr, "df[order(df[,", ref, "], df[,", 
                      c, "]),", c, "])", sep = "")
      }
      eval(parse(text = expr))
      names(df.) <- names(df)
      rownames(df.) <- rownames(df)
      df <- df.
    }
  }
  mat <- matrix(NA, ncol = ncol(df), nrow = nrow(df))
  lcols <- NULL
  for (c in 1:ncol(df)) {
    if (c == 1) {
      offset <- 0
    }
    else {
      offset <- min(mat[, (c - 1)], na.rm = T) + length(levels(df[, 
                                                                  (c - 1)])) - 1
    }
    mat[, c] <- as.numeric(df[, c]) + offset
    cols <- levels(df[, c])
    if ((max(mat[, c], na.rm = T) - offset) != length(cols)) 
      cols <- cols[-length(cols)]
    names(cols) <- cols
    if (is.null(labelscolors)) {
      ncols <- cols
      if (length(cols) < 3) {
        cols <- default2cols
      }
      else {
        cols <- defaultcolfun(length(cols))
      }
      names(cols) <- ncols
      lcols <- c(lcols, cols)
    }
    else {
      if (all(cols %in% names(labelscolors))) {
        cols <- labelscolors[cols]
        lcols <- c(lcols, cols)
      }
      else {
        if (length(cols) < 3) {
          cols <- default2cols
        }
        else {
          cols <- defaultcolfun(length(cols))
        }
        lcols <- c(lcols, cols)
      }
    }
  }
  namessize <- max(strwidth(names(df), units = "inches")) * 
    (cex.annot+.5)
  par(mai = c(0.39375, namessize + 1.6, 0.39375, 0.39375), 
      oma = rep(0, 4))
  image(mat[, ncol(mat):1], col = lcols, axes = F, ...)
  boxl <- par("usr")
  sepl <- seq(boxl[3], boxl[4], ((boxl[4] - boxl[3])/ncol(mat)))
  sepl <- sepl[-c(1, length(sepl))]
  if (addrowdelim) 
    abline(h = sepl, col = "black", lwd = lwd)
  if (addrefdelim) {
    sepg <- ((cumsum(table(mat[, ref])))/(nrow(mat))) * (boxl[2] - 
                                                           boxl[1]) + boxl[1]
    abline(v = sepg, col = "black", lty = 3, lwd = 1)
  }
  box(lwd = lwd)
  if (reorder.method == "all") {
    if (!is.null(addsamplenames)) 
      if (addsamplenames == "top") 
        axis(3, at = seq(0, 1, (1/(nrow(mat) - 1))), 
             labels = rownames(df), tick = FALSE, las = 2, 
             cex.axis = cex.annot)
    else axis(1, at = seq(0, 1, (1/(nrow(mat) - 1))), 
              labels = rownames(df), tick = FALSE, las = 2, 
              cex.axis = cex.annot)
  }
  axis(2, at = seq(0, 1, (1/(ncol(mat) - 1))), labels = rev(names(df)), 
       cex.axis = cex.annot, las = 2, tick = FALSE)
  if (is.null(labelscolors)) 
    lcols
}


# internal
cit.rainbow <- function (n, default2Col = c("black", "red"), default6Col = c("red", 
                                                              "blue", "green", "orange", "purple", "gray"), colorMode = c("rainbow", 
                                                                                                                          "gray", "grey")) 
{
  colorMode = match.arg(colorMode)
  if (colorMode == "grey") 
    colorMode = "gray"
  if (n == 1) 
    return(default2Col[1])
  if (n == 2) 
    return(default2Col)
  if (n > 6 & colorMode == "rainbow") 
    return(rainbow(n))
  if (n < 7 & n > 2 & colorMode == "rainbow") 
    return(default6Col[1:n])
  if (n > 2 & colorMode == "gray") 
    return(gray(seq(1, 0, length = n)))
}



cit.unsplit <- function (x) 
{
  nL <- unlist(lapply(x, length))
  part <- NULL
  for (i in 1:length(x)) part <- c(part, rep(names(x)[i], nL[i]))
  tmp <- cbind(unlist(x), part)
  attr(tmp, "dimnames") <- NULL
  return(tmp)
}
