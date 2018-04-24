# Rd
# description >> internal
# argument
# item >> DATA >> ...
# item >> cl >>  ...
# item >> alternative >> ...
# item >> equalVar >> ...
# item >> conf.level >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.Ttest  <- function (DATA, cl, alternative = c("two.sided", "greater", "less")[1],
    equalVar = F, conf.level = 0.95) {

    if (length(table(cl)) != 2)
        stop("Error -  my.t.test : cl must define exactly 2 (non empty) classes.")
    cl <- as.factor(as.character(cl))
    C1 <- which(as.numeric(cl) == 1)
    C2 <- which(as.numeric(cl) == 2)
    nb.Genes <- nrow(DATA)
    nb.Samples <- ncol(DATA)
    C1.size <- length(C1)
    C2.size <- length(C2)
    stat.C1 <- cit.rowMeansVars(DATA[, C1])
    stat.C2 <- cit.rowMeansVars(DATA[, C2])
    diffmean.C1C2 <- stat.C1$means - stat.C2$means
    if (equalVar) {
        info <- "EQUAL variances hypothesis"
        df <- nb.Samples - 2
        pooledSqrtVar.C1C2 <- sqrt((1/C1.size + 1/C2.size) *
            (stat.C1$sumsquaredif + stat.C2$sumsquaredif)/df)
        tstat <- diffmean.C1C2/pooledSqrtVar.C1C2
    }
    else {
        info <- "UNEQUAL variances hypothesis"
        a1 <- (stat.C1$vars/C1.size)
        a2 <- (stat.C2$vars/C2.size)
        a <- sqrt(a1 + a2)
        df <- a^4/(a1^2/(C1.size - 1) + a2^2/(C2.size - 1))
        tstat <- diffmean.C1C2/a
    }
    if (alternative == "less") {
        pval <- pt(tstat, df)
        cint <- cbind(rep(-Inf, nb.Genes), tstat + qt(conf.level,
            df))
    }
    else if (alternative == "greater") {
        pval <- pt(tstat, df, lower.tail = FALSE)
        cint <- cbind(tstat - qt(conf.level, df), rep(Inf, nb.Genes))
    }
    else {
        pval <- 2 * pt(-abs(tstat), df)
        alpha <- 1 - conf.level
        cint <- qt(1 - alpha/2, df)
        cint <- cbind(tstat - cint, tstat + cint)
    }

        
    res <- as.data.frame(cbind(pval, tstat,cint))
    names(res) <- paste("T-test (uneq. var.)", c("p.values","statistics", "C.I. inf", "C.I. sup"))
    res <- cbind(res, cit.GMFC(DATA, cl))

    list(res,df = df,  info = info)


}

# Rd
# description >> internal
# argument
# item >> DATA >> ...
# item >> cl >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.GMFC <- function (DATA, cl)
{
    groups <- split(dimnames(DATA)[[2]], as.character(cl))
    nag <- names(groups)
    res <- NULL
    for (g in groups) {
        if (length(g) > 1) {
            gm <- rowMeans(DATA[, g], na.rm = T)
        }
        else {
            gm <- DATA[, g]
        }
        if (is.null(res)) {
            res <- gm
        }
        else {
            res <- cbind(res, gm)
        }
    }
    res <- as.data.frame(res)
    nam <- paste("GM", names(groups), sep = ".")
    n <- length(groups)
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            res <- cbind(res, res[, i] - res[, j])
            nam <- c(nam, paste("FC", paste(nag[c(i, j)], collapse = "vs"),
                sep = "."))
        }
    }
    names(res) <- nam
    res
}




# Rd
# description >> internal 
# argument
# item >> d >> ...
# item >> classes >>  ...
# item >> rowCentering >> ...
# item >> rowClassesForAggregation >> ...
# item >> rowClassesToKeep >>  ...
# item >> dist.meth >> ...
# item >> maxDist >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.centroids <- function (d, classes, rowCentering = c(NA, mean.na, median.na)[[3]], 
    rowClassesForAggregation = NULL, rowClassesToKeep = NULL, 
    dist.meth = c("spearman", "euclidian", "maximum", "manhattan", 
        "canberra", "binary", "minkowski", "pearson", "dlda", 
        "dqda","cosine"), maxDist = 0.5, ...){ 

    n <- length(classes)
    N <- dim(d)[2]
    if (n != N) 
        stop("ERROR - cit.centroids: size of classes doesn't correspond to number of columns in data")
    dist.meth <- match.arg(dist.meth)
    if (!is.null(rowClassesForAggregation)) {
        d <- cit.dfAggregate(d, rowClassesForAggregation)
        if (!is.null(rowClassesToKeep)) {
            rowClassesToKeep <- intersect(rowClassesToKeep, rownames(d))
            if (length(rowClassesToKeep) < 2) 
                stop("ERROR - cit.centroids  : less than 2 rowClassesToKeep are in d")
            d <- d[rowClassesToKeep, ]
        }
    }
    if (!suppressWarnings(is.na(rowCentering))) 
        d <- sweep(d, 1, apply(d, 1, rowCentering))
    w <- which(!is.na(classes))
    if (length(w) == 0) 
        stop("ERROR - cit.centroids: classes are undefined (NA)")
    nc <- table(classes)
    m <- cit.dfAggregate(d[, w], classes[w], MARGIN = 2, fAggreg = mean.na)
    v <- cit.dfAggregate(d[, w], classes[w], MARGIN = 2, fAggreg = var.na)
    va <- apply(v, 1, function(vg) sum((nc - 1) * vg)/(length(w) - 
        length(nc)))
    vg <- apply(d[, w], 1, var.na)
    L <- list(mean = m, var = v, `aggregated var` = va, `global var` = vg, 
        rowCentering = rowCentering, rowClassesForAggregation = rowClassesForAggregation, 
        rowClassesToKeep = rowClassesToKeep, centroidsdata = list(samples = data.frame(samplename = colnames(d[, 
            w]), class = classes[w], stringsAsFactors = F), data = d[, 
            w]))
    list(centroids = L, distToCentroids = cit.distToCentroids(d, 
        L, dist.meth = dist.meth, maxDist = maxDist, d.isPretreated = TRUE, 
        ...))
}

# Rd
# description >> internal 
# argument
# item >> d >> ...
# item >> centroids >>  ...
# item >> dist.meth >> ...
# item >> maxDist >> ...
# item >> d.isPretreated >> ...
# item >> sdifftop >> ...
# item >> sdisttocent >>  ...
# item >> verbose >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cit.distToCentroids  <-function (d, centroids, dist.meth = c("spearman", "euclidian", 
    "maximum", "manhattan", "canberra", "binary", "minkowski", 
    "pearson", "dlda", "dqda","cosine"), maxDist = 0.5, d.isPretreated = FALSE, 
    sdifftop = NULL, sdisttocent = NULL, verbose = FALSE) {

    dist.meth <- match.arg(dist.meth)
    srcs <- colnames(d)
    if (!d.isPretreated) {
        if (!is.null(centroids$rowClassesForAggregation)) {
            d <- cit.dfAggregate(d, centroids$rowClassesForAggregation)
            if (!is.null(centroids$rowClassesToKeep)) {
                centroids$rowClassesToKeep <- intersect(centroids$rowClassesToKeep, 
                  rownames(d))
                if (length(centroids$rowClassesToKeep) < 2) 
                  stop("ERROR - cit.distToCentroids: less than 2 rowClassesToKeep are in d")
                d <- d[centroids$rowClassesToKeep, ]
            }
        }
        if (!suppressWarnings(is.na(centroids$rowCentering))) 
            d <- sweep(d, 1, apply(d, 1, centroids$rowCentering))
    }
    if (!is.null(centroids$centroidsdata)) {
        dc <- centroids$centroidsdata$data
        colnames(dc) <- paste("centroid.", colnames(dc), sep = "")
        d <- cbind(d, dc)
    }
    N <- ncol(d)
    n <- ncol(centroids$mean)
    if (dist.meth %in% c("dlda", "dqda")) {
        sumlogvar <- apply(log(centroids$var), 2, sum.na)
        if (dist.meth == "dlda") {
            tmp <- apply(d, 2, function(z) apply((z - centroids$mean)^2/centroids$"aggregated var", 
                2, sum.na))
            tdist <- as.data.frame(t(tmp))
        }
        if (dist.meth == "dqda") {
            tmp <- apply(d, 2, function(z) apply((z - centroids$mean)^2/centroids$var, 
                2, sum.na) + sumlogvar)
            tdist <- as.data.frame(t(tmp))
        }
    }
    else {
        d2 <- t(cbind(d, centroids$mean))
        tdist <- as.matrix(cit.dist(d2, meth = dist.meth, diag = TRUE))
        tdist <- as.data.frame(tdist[1:N, (N + 1):(N + n)])
    }
    rownames(tdist) <- names(d)
    names(tdist) <- names(centroids$mean)
    pred <- apply(tdist, 1, function(z) names(centroids$mean)[which.min(z)])
    mind <- apply(tdist, 1, min)
    difftofirst <- function(x) {
        m <- min(x)
        x - m
    }
    difftop <- apply(tdist, 1, function(x) {
        diff(sort(x)[1:2])
    })
    if (is.null(sdifftop)) {
        if (length(grep("centroid.", names(difftop), value = T)) > 
            0) 
            sdifftop <- round(quantile(difftop[grep("centroid.", 
                names(difftop), value = T)], 0.01, na.rm = T), 
                3)
        else sdifftop <- round(quantile(difftop, 0.01, na.rm = T), 
            3)
        if (verbose) 
            cat("Estimated difftop cut-off  =", cit.sciFormat(sdifftop), 
                "\n")
    }
    predf <- sapply(1:nrow(tdist), function(i) {
        if (difftop[i] <= sdifftop) {
            w <- which(difftofirst(tdist[i, ]) <= sdifftop)
            paste(sort(colnames(tdist[i, ])[w]), collapse = "")
        }
        else {
            colnames(tdist[i, ])[which.min(tdist[i, ])]
        }
    })
    names(predf) <- names(pred)
    if (length(grep("centroid.", names(difftop), value = T)) > 
        0) {
        sam <- grep("centroid.", names(difftop), value = T)
        wsure <- sam[centroids$centroidsdata$samples[match(sub("centroid.", 
            "", sam), centroids$centroidsdata$samples$samplename), 
            "class"] == predf[sam]]
        if (length(wsure) == 0) 
            stop("No concordance between pred and predf on centroids data. sdifftop is too high.\n")
        coresettab <- data.frame(row.names = wsure)
        coresettab$groups <- predf[wsure]
        coresettab <- cbind(coresettab, as.matrix(tdist[wsure, 
            ]))
        coresettab$disttocent <- mind[wsure]
        coresettab$difftop <- apply(tdist[wsure, ], 1, function(x) {
            diff(sort(x)[1:2])
        })
        refcoreset <- NULL
        for (g in names(tdist)) {
            L <- split(coresettab[, g], coresettab[, "groups"] == 
                g)["TRUE"]
            inf <- lapply(L, function(x) c(median(x), max(x), 
                mad(x)))
            refcoreset <- cbind(refcoreset, matrix(unlist(inf), 
                ncol = 1, dimnames = list(c("med", "max", "mad"), 
                  g)))
        }
        if (is.null(sdisttocent)) {
            sdisttocent <- max(round((refcoreset[2, ] - refcoreset[1, 
                ])/refcoreset[3, ]))
            if(verbose){
              cat("Estimated disttocentmad cut-off  =", cit.sciFormat(sdisttocent), 
                  "\n")
            }
            sdisttocent <- refcoreset["med", pred] + sdisttocent * 
                refcoreset["mad", pred]
        }
        else {
            if (length(sdisttocent) == 1) 
                sdisttocent <- refcoreset["med", pred] + sdisttocent * 
                  refcoreset["mad", pred]
        }
        if (verbose) {
            print(refcoreset)
            print(unique(cbind(pred, round(sdisttocent, 3))))
        }
        scoregroup <- c(`TRUE` = "CORE.OUTLIER", `FALSE` = "CORE")[as.character(mind > 
            sdisttocent)]
        scoregroup[which(difftop <= sdifftop)] <- "MIXTE"
        names(scoregroup) <- names(difftop)
        pred2 <- pred
        pred2[!scoregroup %in% "CORE"] <- NA
        if (!is.null(maxDist) & dist.meth %in% c("pearson", "spearman")) {
            pred2[which(mind > maxDist)] <- NA
            scoregroup[which(mind > maxDist)] <- "OUTLIER"
        }
        else {
            maxDist <- NULL
        }
    }
    else {
        scoregroup <- rep("ND", length(pred))
        scoregroup[which(difftop <= sdifftop)] <- "MIXTE"
        pred2 <- pred
        pred2[scoregroup %in% "MIXTE"] <- NA
        if (!is.null(maxDist) & dist.meth %in% c("pearson", "spearman")) {
            pred2[which(mind > maxDist)] <- NA
            scoregroup[which(mind > maxDist)] <- "OUTLIER"
        }
        else {
            pred2[which(mind > quantile(mind, 0.95))] <- NA
            scoregroup[which(mind > quantile(mind, 0.95))] <- "OUTLIER"
            maxDist <- quantile(mind, 0.95)
        }
    }
    tdist <- tdist[srcs, ]
    difftop <- difftop[srcs]
    mind <- mind[srcs]
    pred <- pred[srcs]
    predf <- predf[srcs]
    pred2 <- pred2[srcs]
    scoregroup <- scoregroup[srcs]
    list(dist.scores = tdist, distToNearestCentroid = mind, diffDistTopCentroids = difftop, 
        pred = pred, predwmixed = predf, predCore = pred2, group.confidence = scoregroup, 
        dist.meth = dist.meth, cutoffdiffdist = sdisttocent, 
        cutoffdisttocent = sdifftop, cutoffdistmax = maxDist)
}

# Rd
# description >> internal 
# argument
# item >> object >> ...
# item >> newdata >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
predict.coxph.adr <- function (object, newdata) 
{
    if (!is.null(dim(newdata))) 
        if (min(dim(newdata)) == 1) 
            newdata <- unlist(newdata)
    if (is.vector(newdata)) {
        x <- newdata - object$means
        pred <- x * object$coef
    }
    else {
        x <- sweep(as.matrix(newdata), 2, object$means)
        pred <- x %*% object$coef
    }
    list(pred = pred, formul = cbind(means = object$means, coefs = object$coef))
}

# Rd
# description >> internal 
# argument
# item >> d >> ...
# item >> colEvent >>  ...
# item >> colDelai >> ...
# item >> cutDelai >> ...
# item >> suppressNA >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
cut.delay  <- function (d, colEvent = 1, colDelai = 2, cutDelai = 60, suppressNA = TRUE) 
{
    if (suppressNA) 
        d <- suppress.NA(d, cols = c(colEvent, colDelai))
    w <- which(d[, colDelai] > cutDelai)
    if (length(w) > 0) {
        d[w, colDelai] <- cutDelai
        d[w, colEvent] <- 0
    }
    d
}

# Rd
# description >> internal 
# argument
# item >> modele >> ...
# item >> conf.int >>  ...
# item >> digits >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
getSummary.Cox <-function (modele, conf.int = 0.95, digits = 3) 
{
    logtest <- -2 * (modele$loglik[1] - modele$loglik[2])
    beta <- modele$coefficients
    nabeta <- !(is.na(beta))
    beta2 <- beta[nabeta]
    df <- length(beta2)
    Rsquare <- 1 - exp(-logtest/modele$n)
    maxRsquare <- 1 - exp(2 * modele$loglik[1]/modele$n)
    Likelihood.ratio.test <- logtest
    pvalue.Likelihood.ratio.test <- 1 - pchisq(logtest, df)
    Wald.test <- modele$wald.test
    pvalue.Wald.test <- 1 - pchisq(modele$wald.test, df)
    Logrank.test <- modele$score
    pvalue.Logrank.test <- 1 - pchisq(modele$score, df)
    se <- sqrt(diag(modele$var))
    z <- qnorm((1 + conf.int)/2, 0, 1)
    tmp <- cbind(beta, exp(beta), se, beta/se, signif(1 - pchisq((beta/se)^2, 
        1), digits), exp(-beta), exp(beta - z * se), exp(beta + 
        z * se))
    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
        "se(coef)", "z", "p", "exp(-coef)", paste("exp(coef) lower .", 
            round(100 * conf.int, 2), sep = ""), paste("exp(coef) upper .", 
            round(100 * conf.int, 2), sep = "")))
    return(list(coefs = tmp, df = df, Rsquare = Rsquare, maxRsquare = maxRsquare, 
        Likelihood.ratio.test = Likelihood.ratio.test, pvalue.Likelihood.ratio.test = pvalue.Likelihood.ratio.test, 
        Wald.test = Wald.test, pvalue.Wald.test = pvalue.Wald.test, 
        Logrank.test = Logrank.test, pvalue.Logrank.test = pvalue.Logrank.test, 
        method = modele$method))
}

# Rd
# description >> internal 
# argument
# item >> datannot >> ...
# item >> data >>  ...
# item >> annot >> ...
# item >> vars >> ...
# item >> col.event >> ...
# item >> col.delai >> ...
# item >> cutDelai >> ...
# item >> strat >> ...
# item >> stratVal >> ...
# item >> silent >> ...
# item >> method >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
mycoxph  <- function (datannot = NULL, data = NULL, annot = NULL, vars = NULL, 
    col.event = 1, col.delai = 2, cutDelai = NULL, strat = NULL, 
    stratVal = NULL, silent = TRUE, meth = "breslow") 
{
    
    if (is.null(datannot)) {
        datannot <- get.datannot(data = data, annot = annot, 
            var.annot = c(getNum(col.event, names(annot)), getNum(col.delai, 
                names(annot)), getNum(strat, names(annot))), 
            var.data = vars)
    }
    if (!is.null(cutDelai)) {
        datannot <- cut.delay(datannot, col.event, col.delai, 
            cutDelai)
    }
    if (!is.null(strat)) 
        if (is.na(strat)) 
            strat <- NULL
    if (!is.null(stratVal)) {
        datannot <- cbind(datannot, XstrataX = stratVal)
        strat = "XstrataX"
    }
    if (max(datannot[, col.event], na.rm = TRUE) > 1) {
        if (!silent) {
            cat("Warning - Cox.getModel2 : col.event should be logical")
            cat(" with 1=T=event,0=F=no event ... translate 1->1, 2->0\n")
        }
        if (length(table(datannot[, col.event])) > 2) 
            stop("Error - Cox.getModel2 :  col.event has more than 2 possible values !!")
        datannot[, col.event] <- datannot[, col.event]%%2
    }
    ivars <- getNums(vars, names(datannot))
    icol.event <- getNum(col.event, names(datannot))
    icol.delai <- getNum(col.delai, names(datannot))
    istrata <- getNum(strat, names(datannot))
    names(datannot) <- make.names(names(datannot))
    nvars <- getNames(vars, names(datannot))
    ncol.event <- getName(col.event, names(datannot))
    ncol.delai <- getName(col.delai, names(datannot))
    nstrata <- getName(strat, names(datannot))
    L <- list(time = datannot[, icol.delai], status = datannot[, 
        icol.event])
    for (iv in ivars) {
        L <- c(L, list(unlist(datannot[, iv])))
    }
    names(L)[-(1:2)] <- nvars
    if (!is.null(istrata)) {
        L <- c(L, list(unlist(datannot[, istrata])))
        names(L)[length(L)] <- nstrata
        xstrata <- paste("+strata(", nstrata, ")", sep = "")
    }
    else {
        xstrata <- ""
    }
    x <- paste("Surv(time, status) ~", paste(nvars, collapse = "+"), 
        xstrata, collapse = "")
    modele <- NULL
    tests.Maxpvalue <- NULL
    try(modele <- coxph(eval(as.formula(x)), L, method = meth), 
        silent = silent)
    if (!is.null(modele)) {
        tests.Maxpvalue <- getMaxPvalue(modele)
    }
    return(list(vars = vars, coxModel = modele, scoresHistory = tests.Maxpvalue))
}

# Rd
# description >> internal 
# argument
# item >> modele >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
getMaxPvalue <- function (modele) 
{
    modele.summary <- getSummary.Cox(modele)
    tests.Maxpvalue <- max(modele.summary$pvalue.Likelihood.ratio.test, 
        modele.summary$pvalue.Wald.test, modele.summary$pvalue.Logrank.test)
    return(tests.Maxpvalue)
}


# Rd
# description >> internal 
# argument
# item >> delta >> ...
# item >> threshold >>  ...
# value >> ...
# author >> authors of pamr package
# keyword >> internal
# end
soft.shrink <-function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- sign(delta) * dif * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}

# Rd
# description >> internal 
# argument
# item >> prediction >> ...
# item >> truth >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
successRate <- function (prediction, truth) 
{
    sum(as.numeric(prediction == truth), na.rm = T)/length(truth)
}


# Rd
# description >> internal 
# argument
# item >> score >> ...
# item >> scoreRef >>  ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
isGreater <- function (score, scoreRef) 
{
    score > scoreRef
}

# Rd
# description >> internal 
# argument
# item >> x1 >> ...
# item >> y1 >>  ...
# item >> x2 >> ...
# item >> y2 >> ...
# item >> criteria >> ...
# item >> isBest >> ...
# value >> ...
# author >> A de Reynies
# keyword >> internal
# end
pam.adr <- function(x1,y1,x2=NULL,y2=NULL,criteria= successRate, isBest = isGreater){
      
      # learn the model
      model  <- pamr.train(data = list(x=t(x1),y=as.factor(y1)))  #
      
      # get parametersof the model
      scores <- apply(model$yhat,2,function(pred)  criteria(pred,y1))
      threshold <- model$threshold[which.max(scores)]
      prior <- model$prior
      threshold.scale <- model$threshold.scale
      sd <- model$sd
      centroid.overall <- model$centroid.overall
      centroids <- model$centroids
      se.scale <- model$se.scale      
      delta <- scale((centroids - centroid.overall)/sd, FALSE,threshold.scale * se.scale)
      delta.shrunk <- scale(soft.shrink(delta, threshold), FALSE, 1/(threshold.scale * se.scale))  #
      posid <- drop(abs(delta.shrunk) %*% rep(1, length(prior))) >  0
      centroids <- delta.shrunk[posid,  , drop = FALSE] 
      dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
      names(dd0) <- NULL    
      
      
      # apply model to the complete series
      x12 <- x1
      y12 <- y1
      if(!is.null(x2)){
         S <- setdiff(dimnames(x2)[[1]],dimnames(x1)[[1]])
         if(length(S)>0){
              w <- which(dimnames(x2)[[1]] %in% S)
              x12 <- rbind(x1,x2[w,,drop=FALSE])
              if(is.null(y2)) y2 <- rep(NA,dim(x2)[1])
              y12 <- c(as.character(y1),as.character(y2[w]))              
         }
      }  
      x <- (t(x12) - centroid.overall)/sd
      x <- x[posid,  , drop = FALSE] 
      dd <- t(x) %*% centroids
      dist2centroids <- scale(dd, dd0, FALSE)
      pred <- dimnames(dist2centroids)[[2]][unlist(apply(dist2centroids,1,which.max))]
      
      list( pred.class=factoall(as.data.frame(cbind(dist2centroids,pred=pred,class=y12,train=dimnames(x12)[[1]] %in% dimnames(x1)[[1]]))),
            centroids=centroids,      
            centroid.overall=centroid.overall,
            sd=sd,
            dd0=dd0,
            threshol=threshold)
}
