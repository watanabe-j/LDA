#############################################################
## Main functions
#############################################################

## Function for linear discriminant analysis
lda.c <- function(x, 		# data.frame with columns for variable
        group, 
        prior = table(group) / length(group), 
        sub = 1:ncol(x), 	# subset of variables specified by clumns in x
        aic = TRUE,		# Extract AIC?
        CV = FALSE,		# Conduct leave-one-out cross validation?
        ...) {
    BLS <- function(x) {
        sum(diag(x)) + 0.5 * sum(diag(crossprod(x))) + 0.5 * sum(diag(x))^2
    }
    x <- x
    g <- group
    n <- nrow(x)		# number of objects
    p <- ncol(x)		# number of parameters in full model
    q <- dim(x[, sub, drop = FALSE])[2]			# number of parameters in the sub model
    lev <- lev1 <- levels(g)	# levels of group vector
    k <- length(lev)		# number of groups
    counts <- as.vector(table(g))
    prior <- prior
    sub <- sub
    means <- t(matrix(unlist(by(x, g, colMeans)), p))
    if(is.null(colnames(x))) colnames(x) <- paste("v", 1:p, sep="")
    vname <- colnames(x)
    dimnames(means) <- list(lev, vname)
    means.sub <- means[, sub, drop=FALSE]
    if (!missing(prior)) {		# check for prior
        if (any(prior < 0) || round(sum(prior), 5) != 1) 
            stop("invalid 'prior'")
        if (length(prior) != nlevels(g)) 
            stop("'prior' is of incorrect length")
        prior <- prior[counts > 0L]
    }
    if (any(counts == 0L)) {
        empty <- lev[counts == 0L]
        warning(sprintf(ngettext(length(empty), "group %s is empty", 
            "groups %s are empty"), paste(empty, collapse = " ")), 
            domain = NA)
        lev1 <- lev[counts > 0L]
        g <- factor(g, levels = lev1)
        counts <- as.vector(table(g))
    }
    names(prior) <- names(counts) <- lev1
    W <- matrix(0, p, p)
    for(i in 1:k){
        if(counts[i] - 1 > 1) W <- W + by(x, g, cov)[[i]] * (counts[i] - 1)
    }				# within group cov matrix (multiplied by (n - k) for AIC)
    W1 <- as.matrix(W[sub, sub])
    M <- colSums(counts / sum(counts) * means)	# total mean vector
    D <- t(t(means) - M)
    B <- matrix(0, p, p)
    for(i in 1:k){
        B <- B + (D[i, ] %*% t(D[i, ]) * counts[i])
    }											# between group cov matrix
    T <- W + B									# total cov mat
    if(aic==TRUE) {
        T221 <- T[-sub, -sub] - T[-sub, sub] %*% solve(T[sub, sub]) %*% T[sub, -sub]
        B1 <- as.matrix(B[sub, sub])    
        LL <- n * log(det(W1 / n)) + n * log(det(T221 / n)) + n * p * (1 + log(2 * pi))
        K <- (q * k + p - q + p * (p + 1) / 2)
      ## Bias correction terms for MAIC and HAIC
        invW <- solve(W)
        invW1 <- solve(W1)
        C0 <- invW %*% B %*% solve(diag(p) + invW %*% B)
        C1 <- invW1 %*% B1 %*% solve(diag(q) + invW1 %*% B1)
        he01 <- (1 - p / n) ^ -1 * sum(diag(solve(diag(p) + invW %*% B))) - (1 - p / n) ^ -1 * (p - k + 1)
        he02 <- (1 - p / n) ^ -2 * sum(diag(crossprod(solve(diag(p) + invW %*% B)))) - (1 - p/n) ^ -2 * (p - k + 1)
        he11 <- (1 - q / n) ^ -1 * sum(diag(solve(diag(q) + invW1 %*% B1))) - (1 - q / n) ^ -1 * (q - k + 1)
        he12 <- (1 - q / n) ^ -2 * sum(diag(crossprod(solve(diag(q) + invW1 %*% B1)))) - (1 - q / n) ^ -2 * (q - k + 1)
        K0 <- 2 * (k + 1) * he01 - he02 - he01 ^ 2
        K1 <- 2 * (k + 1) * he11 - he12 - he11 ^ 2
        b0hd <- 2 * K + q * (q + k + 1) * (q + 2 * k + 1) / (n - q - k - 1) +
                        (p + 2) * (p - k + 1) * (p + k + 2) / (n - p - 2) -
                        (q + 2) * (q - k + 1) * (q + k + 2) / (n - q - 2)
        aic <- LL + 2 * K
        bic <- LL + log(n) * K
        maic <- LL + 2 * K - 2 * BLS(C0) + 2 * BLS(C1)
        haic <- LL + b0hd + (n - k - 1) / (n - p - 2) * K0 - (n - k - 1) / (n - q - 2) * K1
    } else {
        LL <- NULL
        bic <- NULL
        maic <- NULL
        haic <- NULL
    }
    a <- -2 * (n - k) * solve(W1) %*% t(means.sub)
    a0 <- rowSums(means.sub %*% solve(W1) * means.sub) * (n - k)
    c.function <- rbind(a, a0)
    temp <- matrix(0, nr = k, nc = k)
    temp1 <- row(temp)
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- col(temp)
    temp2 <- temp2[upper.tri(temp2)]
    d.function <- sapply(1:length(temp1), function(i) (c.function[, temp2[i]] - c.function[, temp1[i]]) / 2)
    C <- sapply(1:length(temp1), function(i) log(prior[temp2[i]] / prior[temp1[i]]))		# constant for prior probability
    colnames(d.function) <- paste(levels(g)[temp1], levels(g)[temp2], sep=":")
    rownames(c.function) <- rownames(d.function) <- c(vname[sub], "constant")
    colnames(c.function) <- levels(g)
    sd.v <- sqrt(diag(W1)/(n - k))		# standard deviation for each variable
    scaled.coef <- d.function[-dim(d.function)[1], ] * sd.v
    score.naive <- as.matrix(cbind(x[, sub, drop=FALSE], 1)) %*% d.function		# discriminant score
    true.score <- score.naive - C	# score corrected for prior probabilities
    Lambda_0 <- det(W) / det(T)			# Wilks' lambda for full model
    Lambda_sub <- det(W1) / det(as.matrix(T[sub, sub]))	# Wilks' lambda for subset model
    Lambda_red <- Lambda_0 / Lambda_sub	# Wilks' lambda for redundancy
    if(CV) {		# Leave-one out cross validation
        sc.cv <- matrix(0, ncol = dim(d.function)[2], nrow = n)
        dimnames(sc.cv) <- list(rownames(x), colnames(d.function))
        class.cv <- factor(n, levels = lev)
        post.cv <- matrix(0, ncol = k, nrow = n)
        colnames(post.cv) <- lev
        rownames(post.cv) <- rownames(x[-sub])
        for(i in 1:n){
            ld <- Recall(x[-i, , drop=FALSE], g[-i], sub = sub, prior = prior, 
                        aic = FALSE, CV = FALSE)
            sc.cv[i, ] <- drop(as.matrix(cbind(x[, sub, drop = FALSE][i, ], 1)) %*% ld$d.function - ld$C)
            post.cv[i, ] <- predict(ld, newdata=x[i, , drop = FALSE])$posterior
            class.cv[i] <- predict(ld, newdata=x[i, , drop = FALSE])$class
        }
        CV <- list(class = class.cv, posterior = post.cv, score = sc.cv)
    }
    structure(list(prior = prior, counts = counts, means = means.sub, 
                   original.means = means, d.function = d.function,
                   c.function = c.function, C = C, score.naive = score.naive, 
                   true.score = true.score, scaled.coef = scaled.coef,
                   sub = sub, lev = lev, N = n, x = x[, sub, drop=FALSE], 
                   gr = g, W = W, W1=W1, B = B, T = T, Lambda_full = Lambda_0, 
                   Lambda_sub = Lambda_sub, Lambda_red = Lambda_red, 
                   minus2.logL = LL, AIC = aic, BIC = bic, MAIC = maic, 
                   HAIC = haic, CV = CV, call = match.call()),
              class = "lda.c")
}

## High speed version for bootstrapping
lda.ch <- function(x, group, prior = table(group) / length(group),
                   sub = 1:ncol(x), CV = FALSE, ...) {
    x <- x
    g <- group
    n <- nrow(x)
    p <- ncol(x)
    q <- dim(x[, sub, drop = FALSE])[2]
    lev <- levels(g)
    k <- length(lev)
    counts <- as.vector(table(g))
    prior <- prior
    sub <- sub
    means <- t(matrix(unlist(by(x, g, colMeans)), p))
    means.sub <- means[, sub, drop=FALSE]
    names(prior) <- names(counts) <- lev
    W <- matrix(0, p, p)
    for(i in 1:k){
        if(counts[i] - 1 > 1) W <- W + by(x, g, cov)[[i]] * (counts[i] - 1)
    }
    W1 <- as.matrix(W[sub, sub])
    a <- -2 * (n - k) * solve(W1) %*% t(means.sub)
    a0 <- rowSums(means.sub %*% solve(W1) * means.sub) * (n - k)
    c.function <- rbind(a, a0)
    temp <- matrix(0, nr = k, nc = k)
    temp1 <- row(temp)
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- col(temp)
    temp2 <- temp2[upper.tri(temp2)]
    d.function <- sapply(1:length(temp1), function(i) (c.function[, temp2[i]] - c.function[, temp1[i]]) / 2)
    if(CV) {		# Leave-one out cross validation
        class.cv <- factor(n, levels = lev)
        for(i in 1:n){
            class.cv[i] <- predict(Recall(x[-i, , drop = FALSE], g[-i],
                                          sub = sub, prior = prior, CV = FALSE),
                                   newdata=x[i, , drop=FALSE])$class
        }
        CV <- list(class = class.cv)
    }
    structure(list(prior = prior, counts = counts, means = means.sub, 
                   d.function = d.function, sub = sub, lev = lev, N = n, 
                   x = x[, sub, drop = FALSE], gr = g, W1 = W1, CV = CV,
                   call = match.call()),
              class = "lda.c")
}

predict.lda.c <- function (object, newdata, prior = object$prior, ...) 
{
    if (!inherits(object, "lda.c")) 
        stop("object not of class \"lda.c\"")
    if (missing(newdata)) {
        newdata <- eval.parent(object$call$x)
        if (!is.null(nas <- object$call$na.action)) 
            newdata <- eval(call(nas, newdata))
    }
    if (is.null(dim(newdata))) 
        dim(newdata) <- c(1L, length(newdata))
    sub <- object$sub
    if(identical(colnames(newdata), colnames(object$means))) {
        x <- as.matrix(newdata)[, , drop = FALSE]
    } else {
        x <- as.matrix(newdata)[, sub, drop = FALSE]
    }
    if (ncol(x) != ncol(object$means)) 
        stop("wrong number of variables")
    if (length(colnames(x)) > 0L && any(colnames(x) != dimnames(object$means)[[2L]])) 
        warning("variable names in 'newdata' do not match those in 'object'")
    k <- length(object$prior)
    if (!missing(prior)) {
        if (any(prior < 0) || round(sum(prior), 5) != 1) 
            stop("invalid 'prior'")
        if (length(prior) != k) 
            stop("'prior' is of incorrect length")
    }
    N <- object$N
    n <- nrow(x)
    temp <- matrix(0, nr=k, nc=k)
    temp1 <- row(temp)
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- col(temp)
    temp2 <- temp2[upper.tri(temp2)]
    C <- sapply(1:length(temp1), function(i) log(prior[temp2[i]]/prior[temp1[i]]))											# constant for prior probability
    d.function <- object$d.function
    true.score <- cbind(x,1) %*% d.function - C
    M <- object$means		# Group means
    V <- object$W1 / (N - k)		# Cov matrix
    posterior <- matrix(0, nrow = n, ncol = k)
    for(i in 1:k) {
       posterior[, i] <- exp(-0.5 * mahalanobis(x, M[i, ], V)) * prior[i]
    }
    dimnames(posterior) <- list(dimnames(x)[[1]], dimnames(M)[[1]])
    posterior <- posterior/rowSums(posterior)		# Divide by sum
    nm <- names(object$prior)
    cl <- factor(nm[max.col(posterior)], levels = object$lev)
    if(sum(is.na(cl)) > 0) {
       cl[is.na(cl)] <- factor(sign(true.score[is.na(cl)]), levels = c(1, -1), labels = object$lev)
    }
    dimnames(posterior) <- list(rownames(x), nm)
    list(class = cl, posterior = posterior, score = true.score)
}

print.lda.c <- function (x, ...) {
    Fs <- Chis <- Fr <- Chir <- FALSE
    Ls <- x$Lambda_sub
    N <- x$N
    p <- ncol(x$means)
    k <- length(x$lev) - 1
    if(!is.null(Ls)) {
        if(min(p, k) == 1) {
            df <- c(max(p, k), N + 1 - p)
            Fs <- (1 - Ls) / Ls * df[2] / df[1]
            p.val <- pf(Fs, df[1], df[2], lower.tail=FALSE)
        } else if(min(p, k) == 2) {
            df <- 2 * c(max(p, k), N + 1 - p)
            Fs <- (1 - sqrt(Ls)) / sqrt(Ls) * df[2] / df[1]
            p.val <- pf(Fs, df[1], df[2], lower.tail=FALSE)
        } else {
            df <- p * k
            Chis <- - (N + 0.5 * (k - p - 1)) * log(Ls)
            p.val <- pchisq(Chis, df, lower.tail=FALSE)
        }
        if(!identical(x$W, x$W1)) {
            Lr <- x$Lambda_red
            p2 <- nrow(x$W)
            q2 <- nrow(x$W1)
            if(min(p2, k) == 1) {
                df2 <- c(p2 - q2, N - p2 - 1)
                Fr <- (1 - Lr) / Lr * df2[2] / df2[1]
                p.valr <- pf(Fr, df2[1], df2[2], lower.tail=FALSE)
            } else if(min(p2, k) == 2) {
                df2 <- 2 * c(max(p2, k), N + 1 - p2)
                Fr <- (1 - sqrt(Lr)) / sqrt(Lr) * df2[2] / df2[1]
                p.valr <- pf(Fr, df2[1], df2[2], lower.tail=FALSE)
            } else {
                dfr <- p2 * k
                Chir <- - (N + 0.5 * (k - p2 - 1)) * log(Ls)
                p.valr <- pchisq(Chis, df2, lower.tail=FALSE)
            }
        }
    }
    if (!is.null(cl <- x$call)) {
        names(cl)[2L] <- ""
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nSubsets of variables:\n")
    print(colnames(x$means), ...)
    cat("\nPrior probabilities of groups:\n")
    print(x$prior, ...)
    cat("\nGroup means:\n")
    print(x$means, ...)
    cat("\nCoefficients of linear discriminants:\n")
    print(x$d.function, ...)
    if(Chis) {
        cat("\nTest for significance:\n Wilks' Lambda:    ", 
        Ls, "\n Chi2 statistic:   ",
        Chis, "\n Degree of freedom:",
        df, "\n p-value:          ",
        p.val, "\n")
    } 
    if(Fs) {
        cat("\nTest for significance:\n Wilks' Lambda:    ", 
        Ls, "\n F statistic:      ",
        Fs, "\n Degree of freedom:",
        df, "\n p-value:          ",
        p.val, "\n")
    }
    if(Chir) {
        cat("\nTest for redundancy of the omitted variable(s):\n Wilks' Lambda:    ", 
        Lr, "\n Chi2 statistic:   ",
        Chir, "\n Degree of freedom:",
        df2, "\n p-value:          ",
        p.valr, "\n")
    } 
    if(Fr) {
        cat("\nTest for redundancy of the omitted variable(s):\n Wilks' Lambda:    ", 
        Lr, "\n F statistic:      ",
        Fr, "\n Degree of freedom:",
        df2, "\n p-value:          ",
        p.valr, "\n")
    }
    cat("\n")
    invisible(x)
}


## Function for linear discriminant analysis with PC scores.
lda.pc <- function(x, 		# data.frame with columns for variable
                   group, 
                   prior = table(group) / length(group), 
                   sub = 1:ncol(x), 	# subset of variables specified by clumns in x
                   spc = seq_len(ncol(x[, sub, drop = FALSE])), 	# subset of PCs
                   aic = TRUE,		# Extract AIC?
                   CV = FALSE,		# Examine leave-one-out cross validation?
                   ...) {
    # function to transform data into PC scores
    pcc <- function(x, R, M){
        result <- sweep(as.matrix(x), 2, M) %*% R
        colnames(result) <- colnames(R)
        return(result)
    }
    BLS <- function(x) {
        sum(diag(x)) + 0.5 * sum(diag(crossprod(x))) + 0.5 * sum(diag(x))^2
    }
    x <- x
    g <- group
    n <- nrow(x)				# number of objects
    p <- ncol(x)				# number of parameters in full model
    q <- dim(x[, sub, drop = FALSE])[2]			# number of parameters in the sub model
    lev <- lev1 <- levels(g)	# levels of group vector
    k <- length(lev)		# number of groups
    counts <- as.vector(table(g))
    prior <- prior
    sub <- sub
    means <- t(matrix(unlist(by(x, g, colMeans)), p))
    if(is.null(colnames(x))) colnames(x) <- paste("v", 1:p, sep="")
    vname <- colnames(x)
    dimnames(means) <- list(lev, vname)
    means.sub <- means[, sub, drop = FALSE]
    if (!missing(prior)) {		# check for prior
        if (any(prior < 0) || round(sum(prior), 5) != 1) 
            stop("invalid 'prior'")
        if (length(prior) != nlevels(g)) 
            stop("'prior' is of incorrect length")
        prior <- prior[counts > 0L]
    }
    if (any(counts == 0L)) {
        empty <- lev[counts == 0L]
        warning(sprintf(ngettext(length(empty), "group %s is empty", 
                                 "groups %s are empty"),
                        paste(empty, collapse = " ")), 
                domain = NA)
        lev1 <- lev[counts > 0L]
        g <- factor(g, levels = lev1)
        counts <- as.vector(table(g))
    }
    names(prior) <- names(counts) <- lev1
    W <- matrix(0, p, p)
    for(i in 1:k){
        if(counts[i] - 1 > 1) W <- W + by(x, g, cov)[[i]] * (counts[i] - 1)
    }		# within group cov matrix (multiplied by (n-k) for AIC)
    W1 <- as.matrix(W[sub, sub])
    M <- colSums(counts / sum(counts) * means)	# total mean vector
    D <- t(t(means) - M)
    B <- matrix(0, p, p)
    for(i in 1:k){
        B <- B + (D[i, ] %*% t(D[i, ]) * counts[i])
    }		# between group cov matrix
    T <- W + B	# total cov mat
    if(aic==TRUE) {
        T221 <- T[-sub, -sub] - T[-sub, sub] %*% solve(T[sub, sub]) %*% T[sub, -sub]
        B1 <- as.matrix(B[sub, sub])    
        LL <- n * log(det(W1 / n)) + n * log(det(T221 / n)) + n * p * (1 + log(2 * pi))
        K <- (q * k + p - q + p * (p + 1) / 2)
      ## Bias correction terms for MAIC and HAIC
        invW <- solve(W)
        invW1 <- solve(W1)
        C0 <- invW %*% B %*% solve(diag(p) + invW %*% B)
        C1 <- invW1 %*% B1 %*% solve(diag(q) + invW1 %*% B1)
        he01 <- (1 - p / n) ^ -1 * sum(diag(solve(diag(p) + invW %*% B))) - (1 - p / n) ^ -1 * (p - k + 1)
        he02 <- (1 - p / n) ^ -2 * sum(diag(crossprod(solve(diag(p) + invW %*% B)))) - (1 - p / n) ^ -2 * (p - k + 1)
        he11 <- (1 - q / n) ^ -1 * sum(diag(solve(diag(q) + invW1 %*% B1))) - (1 - q / n) ^ -1 * (q - k + 1)
        he12 <- (1 - q / n) ^ -2 * sum(diag(crossprod(solve(diag(q) + invW1 %*% B1)))) - (1 - q / n) ^ -2 * (q - k + 1)
        K0 <- 2 * (k + 1) * he01 - he02 - he01 ^ 2
        K1 <- 2 * (k + 1) * he11 - he12 - he11 ^ 2
        b0hd <- 2 * K + q * (q + k + 1) * (q + 2 * k + 1) / (n - q - k - 1) +
                        (p + 2) * (p - k + 1) * (p + k + 2) / (n - p - 2) -
                        (q + 2) * (q - k + 1) * (q + k + 2) / (n - q - 2)
        aic <- LL + 2 * K
        bic <- LL + log(n) * K
        maic <- LL + 2 * K - 2 * BLS(C0) + 2 * BLS(C1)
        haic <- LL + b0hd + (n - k - 1) / (n - p - 2) * K0 - (n - k - 1) / (n - q - 2) * K1
    } else {
        LL <- NULL
        bic <- NULL
        maic <- NULL
        haic <- NULL
    }
    a <- -2 * (n - k) * solve(W1) %*% t(means.sub)
    a0 <- rowSums(means.sub %*% solve(W1) * means.sub) * (n - k)
    c.function <- rbind(a, a0)
    temp <- matrix(0, nr = k, nc = k)
    temp1 <- row(temp)
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- col(temp)
    temp2 <- temp2[upper.tri(temp2)]
    d.function <- sapply(1:length(temp1), function(i) (c.function[, temp2[i]] - c.function[, temp1[i]]) / 2)
    C <- sapply(1:length(temp1), function(i) log(prior[temp2[i]] / prior[temp1[i]]))		# constant for prior probability
    colnames(d.function) <- paste(levels(g)[temp1], levels(g)[temp2], sep=":")
    rownames(c.function) <- rownames(d.function) <- c(colnames(x[, sub, drop=FALSE]), "constant")
    colnames(c.function) <- levels(g)
    sd.v <- sqrt(diag(W1) / (n - k))	# standard deviation for each variable
    scaled.coef <- d.function[-dim(d.function)[1], ] * sd.v
    score.naive <- as.matrix(cbind(x[, sub, drop=FALSE], 1)) %*% d.function		# discriminant score
    true.score <- score.naive - C	# score corrected for prior probabilities
  ## Discrimination with PC
    eva <- eigen(W1 / (n - k))$values	# eigenvalue of pooled cov matrix
    eve <- eigen(W1 / (n - k))$vectors	# eivenvector of pooled cov matrix
    if(eve[1, 1] < 0) eve[, 1] <- -eve[, 1]	#positivize first eigenvector
    dimnames(eve) <- list(colnames(x[, sub, drop = FALSE]), paste("PC", 1:dim(eve)[2], sep=""))
    pcs <- pcc(x[, sub, drop = FALSE], eve, M[sub])	# PC scores
    colnames(pcs) <- paste("PC", 1:q, sep="")
    means.pc <- t(matrix(unlist(by(pcs, g, colMeans)), q))	# group mean PC scores
    dimnames(means.pc) <- list(lev, colnames(pcs))
    means.pc.sub <- means.pc[, spc, drop = FALSE]
    W.pc <- matrix(0, q, q)			# within-group cov matrix of PC
    for(i in 1:k){
        if(counts[i] - 1 > 1) W.pc <- W.pc + by(pcs, g, cov)[[i]] * (counts[i] - 1)
    }		# within group cov matrix (multiplied by (n-k) for AIC)
    W1.pc <- as.matrix(W.pc[spc, spc])
    M.pc <- colSums(counts / sum(counts) * means.pc)	# total mean vector of PC
    D.pc <- t(t(means.pc) - M.pc)
    B.pc <- matrix(0, q, q)
    for(i in 1:k){
        B.pc <- B.pc + (D.pc[i, ] %*% t(D.pc[i, ]) * counts[i])
    }				# between group cov matrix of PC
    T.pc <- W.pc + B.pc		# total cov mat
   #if(!aic == FALSE) {
   #    T221.pc <- T.pc[-spc, -spc] - T.pc[-spc, spc] %*% solve(T.pc[spc, spc]) %*% T.pc[spc, -spc]
   #    aic.pc <- n * log(det(W1.pc / n)) + n * log(det(T221.pc / n)) + n * p * (1 + log(2 * pi)) + 2 * (q * k + p - q + p * (p + 1) / 2)
   #}
    a.pc <- -2 * (n - k) * solve(W1.pc) %*% t(pcc(means.sub, eve[, spc, drop = FALSE], 0))
    a0.pc <- rowSums(means.pc.sub %*% solve(W1.pc) * means.pc.sub) * (n - k)
    c.function.pc <- rbind(a.pc, a0.pc)
    d.function.pc <- sapply(1:length(temp1), function(i) (c.function.pc[, temp2[i]] - c.function.pc[, temp1[i]]) / 2)
    colnames(d.function.pc) <- paste(levels(g)[temp1], levels(g)[temp2], sep=":")
    rownames(c.function.pc) <- rownames(d.function.pc) <- c(colnames(pcs[, spc, drop = FALSE]), "constant")
    colnames(c.function.pc) <- levels(g)
    scaled.coef.pc <- d.function.pc[-dim(d.function.pc)[1], ] * sqrt(eva[spc])
    score.naive.pc <- as.matrix(cbind(pcs[, spc, drop=FALSE], 1)) %*% d.function.pc		# discriminant score
    true.score.pc <- score.naive.pc - C	# score corrected for prior probabilities
  ## Back transformation
    a.or <- eve[, spc, drop=FALSE] %*% a.pc
    a0.or <- -diag(means.sub %*% a.or)/2
    c.function.or <- rbind(a.or, a0.or)
    d.function.or <- sapply(1:length(temp1), function(i) (c.function.or[, temp2[i]] - c.function.or[, temp1[i]]) / 2)
    colnames(d.function.or) <- paste(levels(g)[temp1], levels(g)[temp2], sep=":")
    rownames(c.function.or) <- rownames(d.function.or) <- c(colnames(x[, sub, drop=FALSE]), "constant")
    colnames(c.function.or) <- levels(g)    
    scaled.coef.or <- d.function.or[-dim(d.function.or)[1], ] * sd.v
    score.naive.or <- as.matrix(cbind(x[, sub, drop=FALSE], 1)) %*% d.function.or		# discriminant score
    true.score.or <- score.naive.or - C	# score corrected for prior probabilities
    Lambda_0 <- det(W.pc) / det(T.pc)		# Wilks' lambda for model including full set of PCs
    Lambda_sub <- det(W1.pc) / det(as.matrix(T.pc[spc, spc]))	# Wilks' lambda for subset model
    Lambda_red <- Lambda_0 / Lambda_sub	# Wilks' lambda for redundancy of omitted PCs
    angle <- drop(acos(round((t(d.function.or[-(q + 1), ]) %*% d.function[-(q + 1), ]) / (sqrt(sum(d.function.or[-(q + 1), ]^2) * sum(d.function[-(q + 1)]^2))), 15)) / pi * 180)
    if(CV) {		# Leave-one out cross validation
        sc.cv <- matrix(0, ncol = dim(d.function)[2], nrow = n)
        dimnames(sc.cv) <- list(rownames(x), colnames(d.function))
        class.cv <- factor(n, levels = lev)
        post.cv <- matrix(0, ncol = k, nrow = n)
        colnames(post.cv) <- lev
        rownames(post.cv) <- rownames(x[-sub])
        for(i in 1:n){
            ld <- Recall(x[-i, , drop=FALSE], g[-i], sub = sub, spc = spc, 
                         prior = prior, aic = FALSE, CV = FALSE)
            sc.cv[i, ] <- drop(as.matrix(cbind(x[, sub, drop = FALSE][i, ],1)) %*% ld$d.function.or - ld$C)
            post.cv[i, ] <- predict(ld, newdata=x[i, , drop = FALSE])$posterior
            class.cv[i] <- predict(ld, newdata=x[i, , drop = FALSE])$class
        }
        CV <- list(class = class.cv, posterior = post.cv, score = sc.cv)
    }
      structure(list(prior = prior, counts = counts, means = means.sub, 
                     original.means = means, d.function=d.function,
                     c.function=c.function, C = C, score.naive = score.naive, 
                     true.score = true.score, scaled.coef = scaled.coef,
                     sub = sub, spc = spc, lev = lev, N = n,
                     x = x[, sub, drop=FALSE], gr = g, W = W, W1 = W1,
                     B = B, T = T, Lambda_full = Lambda_0, 
                     Lambda_sub = Lambda_sub, Lambda_red = Lambda_red, 
                     minus2.logL = LL, AIC = aic, BIC = bic, MAIC = maic, 
                     HAIC = haic, CV = CV, e.vector = eve, e.value = eva, 
                     W.pc = W.pc, W1.pc = W1.pc, pc.score = pcs, 
                     means.pc = means.pc.sub, original.means.pc = means.pc, 
                     d.function.pc = d.function.pc, 
                     c.function.pc = c.function.pc, 
                     score.naive.pc = score.naive.pc, 
                     true.score.pc = true.score.pc, 
                     scaled.coef.pc = scaled.coef.pc, 
                     d.function.or = d.function.or, 
                     c.function.or = c.function.or, 
                     score.naive.or = score.naive.or, 
                     true.score.or = true.score.or, 
                     scaled.coef.or = scaled.coef.or, angle = angle,
                     call = match.call()),
                class = "lda.pc")
}
#
## High speed version for bootstrapping
lda.pch <- function(x, 
                    group, 
                    prior = table(group) / length(group), 
                    sub = 1:ncol(x), 	# subset of variables
                    spc = seq_len(ncol(x[, sub, drop = FALSE])), 	# subset of PCs
                    CV = FALSE,	# Conduct leave-one-out cross validation?
                    ...) {
    # function to transform data into PC scores
    pcc <- function(x, R, M){
        result <- sweep(as.matrix(x), 2, M) %*% R
        colnames(result) <- colnames(R)
        return(result)
    }
    x <- x
    g <- group
    n <- nrow(x)			# number of objects
    p <- ncol(x)			# number of parameters in full model
    q <- dim(x[, sub, drop = FALSE])[2]	# number of parameters in the sub model
    lev <- levels(g)		# levels of group vector
    k <- length(lev)		# number of groups
    counts <- as.vector(table(g))
    prior <- prior
    sub <- sub
    means <- t(matrix(unlist(by(x, g, colMeans)), p))
    dimnames(means) <- list(lev, colnames(x))
    means.sub <- means[, sub, drop=FALSE]
    names(prior) <- names(counts) <- lev
    W <- matrix(0, p, p)
    for(i in 1:k){
        if(counts[i] - 1 > 1) W <- W + by(x, g, cov)[[i]] * (counts[i] - 1)
    }		# within group cov matrix (multiplied by (n-k) for AIC)
    W1 <- as.matrix(W[sub, sub])
    M <- colSums(counts / sum(counts) * means)	# total mean vector
    # for true (unscaled) discriminant function, from disc function
    a <- -2 * (n - k) * solve(W1) %*% t(means.sub)
    a0 <- rowSums(means.sub %*% solve(W1) * means.sub) * (n - k)
    c.function <- rbind(a, a0)
    temp <- matrix(0, nr = k, nc = k)
    temp1 <- row(temp)
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- col(temp)
    temp2 <- temp2[upper.tri(temp2)]
    d.function <- sapply(1:length(temp1), function(i) (c.function[, temp2[i]] - c.function[, temp1[i]]) / 2)
    eve <- eigen(W1 / (n - k))$vectors		# eivenvector of pooled cov matrix
    if(eve[1, 1] < 0) eve[,1] <- -eve[, 1]	# positivize first eigenvector
    dimnames(eve) <- list(colnames(x[, sub, drop = FALSE]), paste("PC", 1:dim(eve)[2], sep=""))
    pcs <- pcc(x[, sub, drop = FALSE],  eve, M[sub])				# PC scores
    colnames(pcs) <- paste("PC", 1:q, sep="")
    means.pc <- t(matrix(unlist(by(pcs, g, colMeans)), q))	# group mean PC scores
    dimnames(means.pc) <- list(lev, colnames(pcs))
    means.pc.sub <- means.pc[, spc, drop = FALSE]
    W.pc <- matrix(0, q, q)		# within-group cov matrix of PC
    for(i in 1:k){
        if(counts[i] - 1 > 1) W.pc <- W.pc + by(pcs, g, cov)[[i]] * (counts[i] - 1)
    }		# within group cov matrix (multiplied by (n-k) for AIC)
    W1.pc <- as.matrix(W.pc[spc, spc])
    M.pc <- colSums(counts / sum(counts) * means.pc)	# total mean vector of PC
    a.pc <- -2 * (n - k) * solve(W1.pc) %*% t(pcc(means.sub, 
                                                  eve[, spc, drop = FALSE], 0))
    a0.pc <- rowSums(means.pc.sub %*% solve(W1.pc) * means.pc.sub) * (n - k)
    c.function.pc <- rbind(a.pc, a0.pc)
    d.function.pc <- sapply(1:length(temp1), function(i) (c.function.pc[, temp2[i]] - c.function.pc[, temp1[i]]) / 2)
    if(CV) {		# Leave-one out cross validation
        class.cv <- factor(n, levels = lev)
        for(i in 1:n){
            class.cv[i] <- predict(Recall(x[-i, , drop = FALSE], g[-i],
                                          sub = sub, spc = spc,
                                          prior = prior, CV = FALSE),
                                   newdata = x[i, , drop=FALSE])$class
        }
        CV <- list(class=class.cv)
    }
      structure(list(prior = prior, counts = counts, means = means.sub, 
                     original.means = means, C = C, sub = sub,
                     x = x[, sub, drop = FALSE], spc = spc, lev = lev, N = n, 
                     gr = g, CV = CV, W1.pc = W1.pc, e.vector = eve, 
                     means.pc = means.pc.sub, d.function.pc = d.function.pc, 
                     call = match.call()),
                class = "lda.pc")
}

predict.lda.pc <- function (object, newdata, prior = object$prior, ...) 
{
    # function to transform data into PC scores
    pcc <- function(x, R, M){
        result <- sweep(as.matrix(x), 2, M) %*% R
        colnames(result) <- colnames(R)
        return(result)
    }
    if (!inherits(object, "lda.pc")) 
        stop("object not of class \"lda.pc\"")
    if (missing(newdata)) {
        newdata <- eval.parent(object$call$x)
        if (!is.null(nas <- object$call$na.action)) 
            newdata <- eval(call(nas, newdata))
    }
    if (is.null(dim(newdata))) 
        dim(newdata) <- c(1L, length(newdata))
    sub <- object$sub
    spc <- object$spc
    if(identical(colnames(newdata), colnames(object$means))) {
        x <- as.matrix(newdata)[, , drop=FALSE]
    } else {
        x <- as.matrix(newdata)[, sub, drop=FALSE]
    }
    if (ncol(x) != ncol(object$means)) 
        stop("wrong number of variables")
    if (length(colnames(x)) > 0L && any(colnames(x) != dimnames(object$means)[[2L]])) 
        warning("variable names in 'newdata' do not match those in 'object'")
    k <- length(object$prior)
    if (!missing(prior)) {
        if (any(prior < 0) || round(sum(prior), 5) != 1) 
            stop("invalid 'prior'")
        if (length(prior) != k) 
            stop("'prior' is of incorrect length")
    }
    N <- object$N
    n <- nrow(x)
    temp <- matrix(0, nr = k, nc = k)
    temp1 <- row(temp)
    temp1 <- temp1[upper.tri(temp1)]
    temp2 <- col(temp)
    temp2 <- temp2[upper.tri(temp2)]
    C <- sapply(1:length(temp1), function(i) log(prior[temp2[i]] / prior[temp1[i]]))			# constant for prior probability
    x <- pcc(x, object$e.vector[, spc], colMeans(object$x))
    d.function <- object$d.function.pc
    true.score <- cbind(x, 1) %*% d.function - C
    M <- object$means.pc		# Group means
    V <- object$W1.pc / (N - k)		# Cov matrix
    posterior <- matrix(0, nrow = n, ncol = k)
    for(i in 1:k) {
        posterior[, i] <- exp(-0.5 * mahalanobis(x, M[i, ], V)) * prior[i]
    }
    dimnames(posterior) <- list(dimnames(x)[[1]], dimnames(M)[[1]])
    posterior <- posterior / rowSums(posterior)		# Divide by sum
    nm <- names(object$prior)
    cl <- factor(nm[max.col(posterior)], levels = object$lev)
    if(sum(is.na(cl)) > 0) {
       cl[is.na(cl)] <- factor(sign(true.score[is.na(cl)]), levels=c(1, -1), labels = object$lev)
    }
    dimnames(posterior) <- list(rownames(x), nm)
    list(class = cl, posterior = posterior, score = true.score)
}

print.lda.pc <- function (x, ...) {
    Fs <- Chis <- Fr <- Chir <- FALSE
    Ls <- x$Lambda_sub
    N <- x$N
    p <- ncol(x$means.pc)
    k <- length(x$lev) - 1
    if(!is.null(Ls)) {
        if(min(p, k) == 1) {
            df <- c(max(p, k), N + 1 - p)
            Fs <- (1 - Ls) / Ls * df[2] / df[1]
            p.val <- pf(Fs, df[1], df[2], lower.tail = FALSE)
        } else if(min(p, k) == 2) {
            df <- 2 * c(max(p, k), N + 1 - p)
            Fs <- (1 - sqrt(Ls)) / sqrt(Ls) * df[2] / df[1]
            p.val <- pf(Fs, df[1], df[2], lower.tail = FALSE)
        } else {
            df <- p * k
            Chis <- - (N + 0.5 * (k - p - 1)) * log(Ls)
            p.val <- pchisq(Chis, df, lower.tail = FALSE)
        }
        if(!identical(x$W, x$W1)) {
            Lr <- x$Lambda_red
            p2 <- nrow(x$W)
            q2 <- nrow(x$W1)
            if(min(p2,k) == 1) {
                df2 <- c(p2 - q2, N - p2 - 1)
                Fr <- (1 - Lr) / Lr * df2[2] / df2[1]
                p.valr <- pf(Fr, df2[1], df2[2], lower.tail = FALSE)
            } else if(min(p2, k) == 2) {
                df2 <- 2 * c(max(p2, k), N + 1 - p2)
                Fr <- (1 - sqrt(Lr)) / sqrt(Lr) * df2[2] / df2[1]
                p.valr <- pf(Fr, df2[1], df2[2], lower.tail = FALSE)
            } else {
                dfr <- p2 * k
                Chir <- - (N + 0.5 * (k - p2 - 1)) * log(Ls)
                p.valr <- pchisq(Chis, df2, lower.tail = FALSE)
            }
        }
    }
    if (!is.null(cl <- x$call)) {
        names(cl)[2L] <- ""
        cat("Call:\n")
        dput(cl, control = NULL)
    }
    cat("\nSubsets of variables:\n")
    print(colnames(x$means), ...)
    cat("\nSubsets of PCSs:\n")
    print(colnames(x$means.pc), ...)
    cat("\nPrior probabilities of groups:\n")
    print(x$prior, ...)
    cat("\nGroup means:\n")
    print(x$means, ...)
    cat("\nCoefficients of linear discriminants:\n")
    print(x$d.function.or, ...)
    if(Chis) {
        cat("\nTest for significance:\n Wilks' Lambda:    ", 
        Ls, "\n Chi2 statistic:   ",
        Chis, "\n Degree of freedom:",
        df, "\n p-value:          ",
        p.val, "\n")
    } 
    if(Fs) {
        cat("\nTest for significance:\n Wilks' Lambda:    ", 
        Ls, "\n F statistic:      ",
        Fs, "\n Degree of freedom:",
        df, "\n p-value:          ",
        p.val, "\n")
    }
    if(Chir) {
        cat("\nTest for redundancy of the omitted PC(s):\n Wilks' Lambda:    ", 
        Lr, "\n Chi2 statistic:   ",
        Chir, "\n Degree of freedom:",
        df2, "\n p-value:          ",
        p.valr, "\n")
    } 
    if(Fr) {
        cat("\nTest for redundancy of the omitted PC(s):\n Wilks' Lambda:    ", 
        Lr, "\n F statistic:      ",
        Fr, "\n Degree of freedom:",
        df2, "\n p-value:          ",
        p.valr, "\n")
    }
    cat("\n")
    invisible(x)
}



#############################################################
## Functions for comparisons of models
#############################################################

## Function for extracting AICs and bootstrap error rates from all possible
## submodels (with p > 0) of linear discriminant analysis.
all.aic <- function(data,     # data frame of discriminator
                    group,    # group variable: specify as iris[5], not iris[,5]
                    CV = FALSE, # Cross-validation error rate
                    bcv = FALSE,# Bootstrap cross-validation
                    B = FALSE,  # Bootstrap replication number, default FALSE
                    fixed = FALSE,  # Sample size of each group fixed for Error?
                    seed = FALSE,   # seed to be passed to Error function
                    save = FALSE,
                    filename = "Temp_AIC") {
    BinConv <- function(nv) {
        n <- 2^nv                # how to extract independent variables
        bincomb <- matrix(FALSE, nrow = n, ncol = nv)            # from bincombinations, package e1071
        for (j in 1:nv) {
            bincomb[, j] <- rep(c(rep(FALSE, n / 2^j), rep(TRUE, n / 2^j)), length = n)
        }
        bincomb <- bincomb[rowSums(bincomb) > 0, ]		# exclude the case with no variables
        return(bincomb)
    }
    nv <- ncol(data)
    vname <- colnames(data)       # names of variables used in print method
    if (is.null(vname)) {
          vname <- colnames(data) <- paste("v", 1:nv, sep="")
    }
    group <- factor(as.matrix(group))
    ok <- complete.cases(data, group)   # exclude cases with missing values
    data <- data[ok, ]
    group <- group[ok]
    n <- nrow(data)
    bincomb <- BinConv(nv)
    nr <- nrow(bincomb)
    AIC <- numeric(nr)                  # vector to store AICs
    BIC <- numeric(nr)
    MAIC <- numeric(nr)
    HAIC <- numeric(nr)
    Lambda <- numeric(nr)
    if(CV) e.CV <- numeric(nr)
    if(B) {
        e.BCV <- numeric(nr)            # vector to store BCV error rate
        e.BC.BS <- numeric(nr)          # vector to store BC error rate
        e.632 <- numeric(nr)            # vector to store 632 error rate
        e.632plus <- numeric(nr)        # vector to store 632plus error rate
        if(save) {
            cat("\ne.CV    e.BCV   e.BC.BS e.632   e.632+    formula\n")
        }
    }
    for (i in 1:nr) {
        a <- lda.c(data, group, sub = which(bincomb[i, ]), CV = CV)   # define with lda.c function elsewhere
        if(B) b <- Errors(data[, which(bincomb[i, ]), drop = FALSE], group, fun = lda.ch, B, bcv = bcv, fixed = fixed, seed = seed)
        AIC[i] <- a$AIC               # store AIC
        BIC[i] <- a$BIC
        MAIC[i] <- a$MAIC
        HAIC[i] <- a$HAIC
        Lambda[i] <- a$Lambda_sub
        if(CV) e.CV[i] <- mean(a$CV$class != group)
        if(B) {
            e.BCV[i] <- b$BCV.error.rate
            e.BC.BS[i] <- b$BC.BS.error.rate
            e.632[i] <- b$BS632.error.rate
            e.632plus[i] <- b$BS632plus.error.rate
        }
        if(save && B) {
            try(write.table(cbind(e.CV, e.BCV, e.BC.BS, e.632, e.632plus), 
                            file = paste(filename, ".csv", sep=""),
                            sep=","), silent = TRUE)
            cat(sprintf("%1.5f %1.5f %1.5f %1.5f %1.5f   ~ %s\n", e.CV[i], 
                        e.BCV[i], e.BC.BS[i], e.632[i], e.632plus[i], 
                        paste(vname[as.matrix(bincomb[i, ])],
                        collapse=" + ")))
        }
    }
    deltaAIC <- AIC - min(AIC)           # delta-AIC
    deltaBIC <- BIC - min(BIC)
    deltaMAIC <- MAIC - min(MAIC)
    deltaHAIC <- HAIC - min(HAIC)
    ans <- data.frame(AIC = AIC, deltaAIC = deltaAIC, BIC = BIC,
                      deltaBIC = deltaBIC, MAIC = MAIC,
                      deltaMAIC = deltaMAIC, HAIC = HAIC,
                      deltaHAIC = deltaHAIC, Lambda = Lambda)
    if(CV) ans <- data.frame(ans, e.CV = e.CV)
    if(B) ans <- data.frame(ans, e.BCV = e.BCV, e.BC.BS = e.BC.BS,
                            e.632 = e.632, e.632plus = e.632plus)
    ans <- data.frame(ans, bincomb)
    colnames(ans)[(dim(ans)[2] - nv + 1):dim(ans)[2]] <- vname
    return(structure(list(ans = ans, name = vname, CV = CV, B = B), 
                     class = "all.aic"))
}
# print method
print.all.aic <- function(obj, method = "AIC"){
    ans <- obj$ans
    name <- obj$name
    o <- order(ans[, method], decreasing = FALSE)
    ans <- ans[o, ]
    nc <- ncol(ans)
    if(obj$CV) {
        if(obj$B) {    # Both CV and B
            cat("\n AIC          deltaAIC  BIC          deltaBIC  MAIC        deltaMAIC  HAIC        deltaHAIC    Wilks' L   e.CV       e.BCV      e.BC.BS    e.632      e.632plus formula\n")
            for (i in 1:nrow(ans)) {
                cat(sprintf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f   ~ %s\n", ans[i, 1], ans[i, 2], ans[i, 3], ans[i, 4], ans[i, 5], ans[i, 6], ans[i, 7], ans[i, 8], ans[i, 9], ans[i, 10], ans[i, 11], ans[i, 12], ans[i, 13], ans[i, 14], paste(name[as.matrix(ans[i, 15:nc])], collapse=" + ")))
            }
        } else {    # CV only
            cat("\n AIC          deltaAIC  BIC          deltaBIC  MAIC        deltaMAIC  HAIC        deltaHAIC    Wilks' L   e.CV      formula\n")
            for (i in 1:nrow(ans)) {
                cat(sprintf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f   ~ %s\n", ans[i, 1], ans[i, 2], ans[i, 3], ans[i, 4], ans[i, 5], ans[i, 6], ans[i, 7], ans[i, 8], ans[i, 9], ans[i, 10], paste(name[as.matrix(ans[i, 11:nc])], collapse=" + ")))
            }
        }
    } else {
        if(obj$B) {    # B only
            cat("\n AIC          deltaAIC  BIC          deltaBIC  MAIC        deltaMAIC  HAIC        deltaHAIC    Wilks' L   e.BCV      e.BC.BS    e.632      e.632plus formula\n")
            for(i in 1:nrow(ans)) {
                cat(sprintf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f   ~ %s\n", ans[i, 1], ans[i, 2], ans[i, 3], ans[i, 4], ans[i, 5], ans[i, 6], ans[i, 7], ans[i, 8], ans[i, 9], ans[i, 10], ans[i, 11], ans[i, 12], ans[i, 13], paste(name[as.matrix(ans[i, 14:nc])], collapse=" + ")))
            }
        } else {
            cat("\n AIC          deltaAIC  BIC          deltaBIC  MAIC        deltaMAIC  HAIC        deltaHAIC    Wilks' L   formula\n")
            for (i in 1:nrow(ans)) {
                cat(sprintf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f   ~ %s\n", ans[i, 1], ans[i, 2], ans[i, 3], ans[i, 4], ans[i, 5], ans[i, 6], ans[i, 7], ans[i, 8], ans[i, 9], paste(name[as.matrix(ans[i, 10:nc])], collapse=" + ")))
            }
        }
    }
    invisible(ans)                        # gives only sorted results
}



## Function for extracting AICs and bootstrap error rates from all possible
## submodels (with p > 1) of linear discriminant analysis using PC scores
## except PC1.
all.aic.pc <- function(data,     # data frame of discriminator
                    group,    # group variable: specify as iris[5], not iris[,5]
                    spc = -1, # subset of PCs, by default excluding PC1
                    CV = FALSE, # Cross-validation error rate
                    bcv = FALSE,# Bootstrap cross-validation
                    B = FALSE,  # Bootstrap replication number, default FALSE
                    fixed = FALSE, # Sample size of each group fixed for Error?
                    seed = FALSE, # seed to be passed to Errors.pc function
                    save = FALSE,
                    filename = "Temp_AICpc")
{
        BinConv <- function(nv) {
            n <- 2^nv                # how to extract independent variables
            bincomb <- matrix(FALSE, nrow=n, ncol=nv)            # from bincombinations, package e1071
            for (j in 1:nv) {
                bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
            }
            bincomb <- bincomb[rowSums(bincomb)>1,]		# exclude the case with no variables
            return(bincomb)
        }
        nv <- ncol(data)
        vname <- colnames(data)       # names of variables used in print method
        if (is.null(vname)) {
              vname <- colnames(data) <- paste("v", 1:nv, sep="")
        }
        group <- factor(as.matrix(group))
        ok <- complete.cases(data, group)   # exclude cases with missing values
        data <- data[ok,]
        group <- group[ok]
        n <- nrow(data)
        bincomb <- BinConv(nv)
        nr <- nrow(bincomb)
        AIC <- numeric(nr)                  # vector to store AICs
        BIC <- numeric(nr)
        MAIC <- numeric(nr)
        HAIC <- numeric(nr)
        Lambda <- numeric(nr)
        if(CV) e.CV <- numeric(nr)
        if(B) {
            e.BCV <- numeric(nr)            # vector to store BCV error rate
            e.BC.BS <- numeric(nr)          # vector to store BC error rate
            e.632 <- numeric(nr)            # vector to store 632 error rate
            e.632plus <- numeric(nr)        # vector to store 632plus error rate
            if(save) {
                cat("\ne.CV    e.BCV   e.BC.BS e.632   e.632+    formula\n")
            }
        }
        for (i in 1:nr) {
             a <- lda.pc(data, group, sub = which(bincomb[i, ]), spc = spc, CV = CV)
             if(B) b <- Errors.pc(data[, which(bincomb[i, ]), drop = FALSE], 
                                  group, spc = spc, B, bcv = bcv,
                                  fixed = fixed, seed = seed)
             AIC[i] <- a$AIC                # store AIC
             BIC[i] <- a$BIC
             MAIC[i] <- a$MAIC
             HAIC[i] <- a$HAIC
             Lambda[i] <- a$Lambda_sub
             if(CV) e.CV[i] <- 1 - sum(a$CV$class == group)/n
             if(B) {
                 e.BCV[i] <- b$BCV.error.rate
                 e.BC.BS[i] <- b$BC.BS.error.rate
                 e.632[i] <- b$BS632.error.rate
                 e.632plus[i] <- b$BS632plus.error.rate
             }
             if(save && B) {
                 try(write.table(cbind(e.CV, e.BCV, e.BC.BS, e.632, e.632plus), 
                                 file=paste(filename, ".csv", sep=""), sep=","),
                     silent=TRUE)
                 cat(sprintf("%1.5f %1.5f %1.5f %1.5f %1.5f   ~ %s\n", e.CV[i], 
                     e.BCV[i], e.BC.BS[i], e.632[i], e.632plus[i], 
                     paste(vname[as.matrix(bincomb[i, ])], collapse=" + ")))
             }
        }
        deltaAIC <- AIC - min(AIC)           # delta-AIC
        deltaBIC <- BIC - min(BIC)
        deltaMAIC <- MAIC - min(MAIC)
        deltaHAIC <- HAIC - min(HAIC)
        ans <- data.frame(AIC=AIC, deltaAIC=deltaAIC, BIC=BIC, deltaBIC=deltaBIC, MAIC=MAIC, deltaMAIC=deltaMAIC, HAIC=HAIC, deltaHAIC=deltaHAIC, Lambda=Lambda)
        if(CV) ans <- data.frame(ans, e.CV=e.CV)
        if(B) ans <- data.frame(ans, e.BCV=e.BCV, e.BC.BS=e.BC.BS, e.632=e.632, e.632plus=e.632plus)
        ans <- data.frame(ans, bincomb)
        colnames(ans)[(dim(ans)[2]-nv+1):dim(ans)[2]] <- vname
        return(structure(list(ans=ans, name=vname, CV=CV, B=B), class="all.aic"))
}


## Function for extracting discriminant coefficients of all possible
## submodels (with p > 0) of linear discriminant analysis.
all.coef <- function(data,     # data frame of discriminator
                     group,    # group memberships
                     B = FALSE,  # Bootstrap replication number, default FALSE
                     fixed = FALSE,  # Sample size of each group fixed?
                     seed = FALSE)   # seed to be passed to Coef.B function
{
    BinConv <- function(nv) {
        n <- 2^nv                # how to extract independent variables
        bincomb <- matrix(FALSE, nrow=n, ncol=nv)            # from bincombinations, package e1071
        for (j in 1:nv) {
            bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
        }
        bincomb <- bincomb[rowSums(bincomb)>0,]		# exclude the case with no variables
        return(bincomb)
    }
    nv <- ncol(data)
    vname <- colnames(data)       # names of variables used in print method
    if (is.null(vname)) {
          vname <- colnames(data) <- paste("x", 1:nv, sep="")
    }
   #group <- factor(as.matrix(group))
    ok <- complete.cases(data, group)   # exclude cases with missing values
    data <- data[ok,]
    group <- group[ok]
    n <- nrow(data)
    bincomb <- BinConv(nv)
    colnames(bincomb) <- vname
    nr <- nrow(bincomb)        
    dcoef <- matrix(NA, ncol = nv + 1, nrow = nr) # matrix to store coefs
    colnames(dcoef) <- c(vname, "constant")
    if(B) {
        CI95 <- matrix(NA, ncol = 2*(nv + 1), nrow = nr)
        if(B >= 200) {
           CI99 <- matrix(NA, ncol = 2*(nv + 1), nrow = nr)
        }
    }
    for (i in 1:nr) {
        a <- lda.c(data, group, sub=which(bincomb[i,]))   # define with lda.c function elsewhere
        if(B) b <- Coef.B(data[,which(bincomb[i,]),drop=FALSE], group, B = B, fixed=fixed, seed=seed)
        dcoef[i,] <- as.data.frame(a$d.function)[colnames(dcoef),]
        if(B) {
            CI95[i,c(2*(1:(nv+1))-1)] <- as.data.frame(t(b$CI95))[colnames(dcoef),1]
            CI95[i,c(2*(1:(nv+1)))] <- as.data.frame(t(b$CI95))[colnames(dcoef),2]
            if(B >= 200) {
                CI99[i,c(2*(1:(nv+1))-1)] <- as.data.frame(t(b$CI99))[colnames(dcoef),1]
                CI99[i,c(2*(1:(nv+1)))] <- as.data.frame(t(b$CI99))[colnames(dcoef),2]
             }
        }
    }
    dcoef[is.na(dcoef)] <- 0           # Put 0s for excluded variables
    colnames(dcoef)[1:nv] <- paste(colnames(dcoef)[1:nv], "c", sep="")
    if(B) { 
        CI95[is.na(CI95)] <- 0
        if(B >= 200) {
            CI99[is.na(CI99)] <- 0
            colnames(CI99) <- paste(rep(colnames(dcoef),each=2), c("l","u"), sep=".")
            CI99 <- data.frame(CI99, bincomb)
        } else CI99 <- NULL
        colnames(CI95) <- paste(rep(colnames(dcoef),each=2), c("l","u"), sep=".")
        CI95 <- data.frame(CI95, bincomb)
    } else CI95 <- NULL
    ans <- data.frame(dcoef, bincomb)
    if(B) return(structure(list(ans = ans, CI95 = CI95, CI99 = CI99, 
                                name = vname, B = B), class="all.coef"))
    else return(structure(list(ans = ans, name = vname, B = B),
                          class = "all.coef"))
}
# print method
print.all.coef <- function(obj) {
    ans <- obj$ans
    name <- obj$name
    p <- length(name)
    nc <- ncol(ans)
    cat(" ", paste(colnames(ans[1:(p + 1)]), collapse = "       "), "formula\n", sep="  ")
    for (i in 1:nrow(ans)) {
        cat(sprintf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f   ~ %s\n", ans[i, 1], ans[i, 2], ans[i, 3], ans[i, 4], ans[i, 5], ans[i, 6], ans[i, 7], ans[i, 8], paste(name[as.matrix(ans[i, 9:nc])], collapse=" + ")))
    }
    invisible(ans)                        # gives only sorted results
}


## Function for extracting coefficients of discriminant function based on 
## PC scores, except PC1, of all possible submodels (with p > 1).
all.coef.pc <- function(data,     # data frame of discriminator
                        group,    # group membership
                        B = FALSE,  # Bootstrap replication number
                        fixed=FALSE, # Sample size of each group fixed?
                        seed=FALSE) {
    BinConv <- function(nv) {
        n <- 2^nv                # how to extract independent variables
        bincomb <- matrix(FALSE, nrow=n, ncol=nv)            # from bincombinations, package e1071
        for (j in 1:nv) {
            bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
        }
        bincomb <- bincomb[rowSums(bincomb)>1,]		# exclude the case with no variables
        return(bincomb)
    }
    nv <- ncol(data)
    vname <- colnames(data)       # names of variables used in print method
    if (is.null(vname)) {
        vname <- colnames(data) <- paste("x", 1:nv, sep="")
    }
    ok <- complete.cases(data, group)   # exclude cases with missing values
    data <- data[ok, ]
    group <- group[ok]
    n <- nrow(data)
    bincomb <- BinConv(nv)
    colnames(bincomb) <- vname
    nr <- nrow(bincomb)        
    dcoef <- matrix(NA, ncol = nv + 1, nrow = nr) # matrix to store coefs
    colnames(dcoef) <- c(vname, "constant")
    angles <- numeric(nr)
    if(B) {
        CI95 <- matrix(NA, ncol = 2 * (nv + 1), nrow = nr)
        if(B >= 200) {
            CI99 <- matrix(NA, ncol = 2*(nv + 1), nrow = nr)
        }
    }
    for (i in 1:nr) {
        a <- lda.pc(data, group, sub = which(bincomb[i, ]), spc = -1)   # define with lda.pc function elsewhere
        if(B) b <- Coef.pc.B(data[, which(bincomb[i, ]), drop=FALSE], group, spc = 2:sum(bincomb[i, ]), B = B, fixed = fixed, seed = seed)
        dcoef[i, ] <- as.data.frame(a$d.function.or)[colnames(dcoef), ]
        angles[i] <- a$angle
        if(B) {
            CI95[i, c(2 * (1:(nv + 1)) - 1)] <- as.data.frame(t(b$CI95))[colnames(dcoef), 1]
            CI95[i, c(2 * (1:(nv + 1)))] <- as.data.frame(t(b$CI95))[colnames(dcoef), 2]
            if(B >= 200) {
                CI99[i, c(2 * (1:(nv + 1)) - 1)] <- as.data.frame(t(b$CI99))[colnames(dcoef), 1]
                CI99[i, c(2 * (1:(nv + 1)))] <- as.data.frame(t(b$CI99))[colnames(dcoef), 2]
            }
        }
    }
    dcoef[is.na(dcoef)] <- 0           # Put 0s for excluded variables
    colnames(dcoef)[1:nv] <- paste(colnames(dcoef)[1:nv], "c", sep="")
    if(B) { 
        CI95[is.na(CI95)] <- 0
        if(B >= 200) {
            CI99[is.na(CI99)] <- 0
            colnames(CI99) <- paste(rep(colnames(dcoef), each=2), c("l","u"), sep=".")
            CI99 <- data.frame(CI99, bincomb)
        } else CI99 <- NULL
        colnames(CI95) <- paste(rep(colnames(dcoef), each=2), c("l","u"), sep=".")
        CI95 <- data.frame(CI95, bincomb)
    } else CI95 <- NULL
    ans <- data.frame(dcoef, angles, bincomb)
    if(B) return(structure(list(ans = ans, CI95 = CI95, CI99 = CI99,
                                name = vname, B = B), class = "all.coef"))
    else return(structure(list(ans = ans, name = vname, B = B),
                          class = "all.coef.pc"))
}
# print method
print.all.coef.pc <- function(obj) {
    ans <- obj$ans
    name <- obj$name
    p <- length(name)
    nc <- ncol(ans)
    cat(" ", paste(colnames(ans[1:(p + 1)]), collapse="       "),
        " angle", "   formula\n", sep="  ")
    for (i in 1:nrow(ans)) {
        cat(sprintf("%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f   ~ %s\n", ans[i, 1], ans[i, 2], ans[i, 3], ans[i, 4], ans[i, 5], ans[i, 6], ans[i, 7], ans[i, 8], ans[i, 9],
                    paste(name[as.matrix(ans[i, 10:nc])],
                          collapse=" + ")))
    }
    invisible(ans)                        # gives only sorted results
}


## Function to return predictions of all possible submodels with lda.c function.
all.class <- function(data,     # data frame of discriminator
                      group,    # group membership
                      newdata,
                      CV = FALSE)  # Leave-one-out cross-validation
{
    if(!missing(newdata) & CV) stop("CV option is not available with newdata")
    BinConv <- function(nv) {
        n <- 2^nv                # how to extract independent variables
        bincomb <- matrix(FALSE, nrow=n, ncol=nv)            # from bincombinations, package e1071
        for (j in 1:nv) {
            bincomb[, j] <- rep(c(rep(FALSE, n/2^j), rep(TRUE, n/2^j)), length = n)
        }
        bincomb <- bincomb[rowSums(bincomb)>0,]		# exclude the case with no variables
        return(bincomb)
    }
   nv <- ncol(data)
   vname <- colnames(data)       # names of variables used in print method
   if (is.null(vname)) {
         vname <- colnames(data) <- paste("x", 1:nv, sep="")
   }
   ok <- complete.cases(data, group)   # exclude cases with missing values
   data <- data[ok, ]
   group <- group[ok]
   if(missing(newdata)) newdata <- data        # replace newdata when missing
       n <- nrow(newdata)
       bincomb <- BinConv(nv)
       nr <- nrow(bincomb)
       result <- as.data.frame(matrix(NA, ncol=n, nrow=nr))    # matrix to store results
       if(CV) r.CV <- as.data.frame(matrix(NA, ncol=n, nrow=nr))
       for (i in 1:nr) {
           a <- lda.c(data, group, sub = which(bincomb[i, ]), CV = CV)   # define with lda.c function elsewhere
           result[i, ] <- predict(a, newdata = newdata)$class
           if(CV) r.CV[i, ] <- a$CV$class
       }
       ans <- data.frame(result)
       if(CV) ans <- data.frame(ans, r.CV)
       ans <- data.frame(ans, bincomb)
       colnames(ans) <- c(rownames(newdata), if(CV) paste(rownames(data), ".CV", sep="") else NULL, vname)
       return(structure(list(ans = ans, name = vname, CV = CV), 
                        class="all.class"))
}
# print method
print.all.class <- function(obj) {
message('Use write.table(...$ans, "....csv", sep=",") to see results')
}

## PC version of all.class
all.class.pc <- function(data,     # data frame of discriminator
                         group,    # group membership
                         newdata,
                         CV = FALSE)  # Leave-one-out cross-validation
{
    if(!missing(newdata) & CV) stop("CV option is not available with newdata")
    BinConv <- function(nv) {
        n <- 2^nv                # how to extract independent variables
        bincomb <- matrix(FALSE, nrow = n, ncol = nv)            # from bincombinations, package e1071
        for (j in 1:nv) {
            bincomb[, j] <- rep(c(rep(FALSE, n / 2^j), rep(TRUE, n / 2^j)), length = n)
        }
        bincomb <- bincomb[rowSums(bincomb)>1,]		# exclude the case with less than 2 variables
        return(bincomb)
    }
    nv <- ncol(data)
    vname <- colnames(data)       # names of variables used in print method
    if (is.null(vname)) {
          vname <- colnames(data) <- paste("x", 1:nv, sep="")
    }
    ok <- complete.cases(data, group)   # exclude cases with missing values
    data <- data[ok, ]
    group <- group[ok]
    if(missing(newdata)) newdata <- data        # replace newdata when missing
    n <- nrow(newdata)
    bincomb <- BinConv(nv)
    nr <- nrow(bincomb)
    result <- as.data.frame(matrix(NA, ncol = n, nrow = nr))    # matrix to store results
    if(CV) r.CV <- as.data.frame(matrix(NA, ncol = n, nrow = nr))
    for (i in 1:nr) {
         a <- lda.pc(data, group, sub = which(bincomb[i, ]), spc = 2:sum(bincomb[i, ]), CV = CV)   # define with lda.pc function elsewhere
         result[i, ] <- predict(a, newdata=newdata)$class
         if(CV) r.CV[i, ] <- a$CV$class
    }
    ans <- data.frame(result)
    if(CV) ans <- data.frame(ans, r.CV)
    ans <- data.frame(ans, bincomb)
    colnames(ans) <- c(rownames(newdata), if(CV) paste(rownames(data), ".CV", sep="") else NULL, vname)
    return(structure(list(ans = ans, name = vname, CV = CV), 
                     class = "all.class.pc"))
}
# print method
print.all.class.pc <- function(obj) {
    message('Use write.table(...$ans, "....csv", sep=",") to see results')
}


#############################################################
## Functions mainly used internally
#############################################################

## Internal function for all.aic to extract bootstrap estimates of error rate
Errors <- function(x, group, fun = lda.ch, B = 1000, bcv = FALSE, 
                   fixed = FALSE, seed=FALSE, ...) {
    E632_plus <- function(x, group, fun=lda.ch, e.1, e.A1, ...) {	# function for .632+
       lev <- levels(group)
       obj <- fun(x, group, ...)
       pr.cl <- predict(obj, newdata = x)$class
       p <- obj$prior
       q <- numeric(length(lev))
       for(i in 1:length(lev)){
          q[i] <- sum(pr.cl == lev[i]) / length(group)
       }
       Gm <- sum(p * (1 - q))    # expected prediction error rate: Eq. 27
       R <- if (e.1 > e.A1 & Gm > e.A1) (e.1 - e.A1) / (Gm - e.A1)
           else 0 # Relative overfitting rate: Eq. 28, 31
       ww <- (1 - 1 / exp(1))/(1 - 1 / exp(1) * R)  # Weights: Eq. 29
       (1 - ww) * e.A1 + ww * e.1  # Eq. 29
    }    
    e.A1 <- sum(predict(fun(x, group, ...), newdata = x)$class != group) / nrow(x)
    N <- nrow(x)
    B <- B
    i <- 1
    e.Ik <- matrix(NA, ncol = B, nrow = N)	# Matrix of non-inclusion
    e.Qk <- matrix(NA, ncol = B, nrow = N)	# Matrix of misclassification
    e.Ak <- matrix(NA, ncol = B, nrow = N)	# Matrix of apparent-misclassification
    if(bcv) {
        e.Ck <- matrix(NA, ncol = B, nrow = N)	# Matrix for cross-validation
    }
    if(fixed) {				# Use fixed number of subsamples for each group
        x.1 <- x[as.numeric(group) == 1, , drop=FALSE]
        x.2 <- x[as.numeric(group) == 2, , drop=FALSE]
        x <- rbind(x.1, x.2)
        group <- factor(c(group[as.numeric(group) == 1],
                         group[as.numeric(group) == 2]), labels = levels(group))
        if(seed) set.seed(seed)
        while(i <= B) {
            S.1 <- sample(1:dim(x.1)[1], replace=TRUE)
            S.2 <- sample(1:dim(x.2)[1], replace=TRUE)
            X <- rbind(x.1[S.1, , drop=FALSE], x.2[S.2, , drop=FALSE])	# BSS
            G <- group[c(S.1, dim(x.1)[1] + S.2)]
            e.Ik[, i] <- as.numeric(!rownames(x) %in% rownames(X))
            CR.B <- fun(X, G, CV = bcv, ...)
            e.Ak[, i] <- as.numeric(predict(CR.B, newdata=X)$class != G)
            e.Qk[, i] <- as.numeric(predict(CR.B, newdata=x)$class != group)
            if(bcv) { 
                e.Ck[,i] <- as.numeric(CR.B$CV$class != G)
            }
            i <- i + 1
        }
    } else {
        pri <- table(group) / length(group)
        if(seed) set.seed(seed)
        while(i <= B) {
            S <- sample(1:N, replace=TRUE)
            X <- x[S, , drop=FALSE]	# BSS
            G <- group[S]
            if(any(table(G) <= ifelse(bcv, 1, 0))) next
            e.Ik[, i] <- as.numeric(!rownames(x) %in% rownames(X))
            CR.B <- fun(X, G, CV = bcv, prior = pri, ...)
            e.Ak[, i] <- as.numeric(predict(CR.B, newdata = X)$class != G)
            e.Qk[, i] <- as.numeric(predict(CR.B, newdata = x)$class != group)
            if(bcv) { 
                e.Ck[, i] <- as.numeric(CR.B$CV$class != G)
            }
            i <- i + 1
        }
    }
    EX.e <- colSums(e.Ik) != 0 # All objects included in a set? If so, exclude it.
    e.I <- rowSums(e.Ik[, EX.e] * e.Qk[, EX.e]) / rowSums(e.Ik[, EX.e])
    names(e.I) <- rownames(x)
    e.1 <- mean(e.I, na.rm=TRUE)
    e.Boot <- mean(e.Qk - e.Ak) + e.A1		# Bias-corrected BS error rate
    e.A <- mean(e.Ak)
    e.BCV <- ifelse(bcv, mean(e.Ck), NA)
    e632 <-  (1 - 1 / exp(1)) * e.1 + (1 / exp(1)) * e.A1
    e632plus <- E632_plus(x, group, fun, e.1 = e.1, e.A1 =e.A1, ...)
    list(BCV.error.rate = e.BCV, BS.apparent.error = e.A,
         BC.BS.error.rate = e.Boot, BS632.error.rate = e632,
         BS632plus.error.rate = e632plus)
}


## Errors function modified for lda.pc, discrimination with PC scores except PC1
Errors.pc <- function(x, group, spc = -1, B = 1000, bcv = FALSE, fixed = FALSE,
                      seed = FALSE, ...) {
    E632_plus.pc <- function(x, group, spc, e.1, e.A1, ...) {	# function for .632+
       lev <- levels(group)
       obj <- lda.pch(x, group, spc = spc, ...)
       pr.cl <- predict(obj, newdata = x)$class
       p <- obj$prior
       q <- numeric(length(lev))
       for(i in 1:length(lev)){
          q[i] <- sum(pr.cl == lev[i]) / length(group)
       }
       Gm <- sum(p * (1 - q))    # expected prediction error rate: Eq. 27
       R <- if (e.1 > e.A1 & Gm > e.A1) (e.1 - e.A1) / (Gm - e.A1)
           else 0 # Relative overfitting rate: Eq. 28, 31
       ww <- (1 - 1 / exp(1))/(1 - 1 / exp(1) * R)  # Weights: Eq. 29
       (1 - ww) * e.A1 + ww * e.1  # Eq. 29
    }    
    e.A1 <- sum(predict(lda.pch(x, group, spc = spc))$class != group) / nrow(x)
    N <- nrow(x)
    B <- B
    i <- 1
    e.Ik <- matrix(NA, ncol = B, nrow = N)	# Matrix of non-inclusion
    e.Qk <- matrix(NA, ncol = B, nrow = N)	# Matrix of misclassification
    e.Ak <- matrix(NA, ncol = B, nrow = N)	# Matrix of apparent-misclassification
    e.Ck <- matrix(NA, ncol = B, nrow = N)	# Matrix for cross-validation
    if(fixed) {				# Use fixed number of subsamples for each group
        x.1 <- x[as.numeric(group) == 1, , drop = FALSE]
        x.2 <- x[as.numeric(group) == 2, , drop = FALSE]
        x <- rbind(x.1, x.2)
        group <- factor(c(group[as.numeric(group) == 1],
                         group[as.numeric(group) == 2]), labels = levels(group))
        if(seed) set.seed(seed)
        while(i <= B) {
            S.1 <- sample(1:dim(x.1)[1], replace=TRUE)
            S.2 <- sample(1:dim(x.2)[1], replace=TRUE)
            X <- rbind(x.1[S.1, , drop=FALSE], x.2[S.2, , drop=FALSE])	# BSS
            G <- group[c(S.1, dim(x.1)[1] + S.2)]
            e.Ik[, i] <- as.numeric(!rownames(x) %in% rownames(X))
            CR.B <- lda.pch(X, G, spc = spc, CV = bcv, ...)
            e.Ak[, i] <- as.numeric(predict(CR.B, newdata = X)$class != G)
            e.Qk[, i] <- as.numeric(predict(CR.B, newdata = x)$class != group)
            e.Ck[, i] <- as.numeric(CR.B$CV$class != G)
            i <- i + 1
        }
    } else {
        pri <- table(group) / length(group)
        if(seed) set.seed(seed)
        while(i <= B) {
            S <- sample(1:N, replace = TRUE)
            X <- x[S, , drop=FALSE]	# BSS
            G <- group[S]
            if(any(table(G) <= ifelse(bcv, 1, 0))) next
            e.Ik[, i] <- as.numeric(!rownames(x) %in% rownames(X))
            CR.B <- lda.pch(X, G, spc = spc, CV = bcv, prior = pri, ...)
            e.Ak[, i] <- as.numeric(predict(CR.B, newdata = X)$class != G)
            e.Qk[, i] <- as.numeric(predict(CR.B, newdata = x)$class != group)
            e.Ck[, i] <- as.numeric(CR.B$CV$class != G)
            i <- i + 1
        }
    }
    EX.e <- colSums(e.Ik) != 0 # All objects included in a set? If so, exclude it.
    e.I <- rowSums(e.Ik[, EX.e] * e.Qk[, EX.e]) / rowSums(e.Ik[, EX.e])
    names(e.I) <- rownames(x)
    e.1 <- mean(e.I)
    e.Boot <- mean(e.Qk - e.Ak) + e.A1		# Bias-corrected BS error rate
    e.A <- mean(e.Ak, na.rm = TRUE)
    e.BCV <- ifelse(bcv, mean(e.Ck), NA)
    e632 <-  (1 - 1 / exp(1)) * e.1 + (1 / exp(1)) * e.A1
    e632plus <- E632_plus.pc(x, group, spc = spc, e.1 = e.1, e.A1 = e.A1, ...)
    list(BCV.error.rate = e.BCV, BS.apparent.error = e.A,
         BC.BS.error.rate = e.Boot, BS632.error.rate = e632,
         BS632plus.error.rate = e632plus)
}


## Function to extract variance-covariance matrix, a summary (range), and
## 95% and 99% confidence intervals of coefficients of discriminant function
## of linear discriminant analysis with bootstrapping.
Coef.B <- function(x, group, sub = 1:ncol(x), B = 1000, fixed = FALSE,
                   seed = FALSE) {
    require(rrcov)
    N <- nrow(x)
    k <- length(levels(group))
    B <- B
    coef.B <- matrix(NA, ncol = B, nrow = ncol(x[sub]) + 1)
    rownames(coef.B) <- rownames(lda.c(x, group, sub = sub, AIC = FALSE)$d.function)
    i <- 1
    if(fixed) {				# Use fixed number of subsamples for each group
        x.1 <- x[as.numeric(group) == 1, , drop=FALSE]
        x.2 <- x[as.numeric(group) == 2, , drop=FALSE]
        x <- rbind(x.1, x.2)
        group <- factor(c(group[as.numeric(group) == 1],
                          group[as.numeric(group) == 2]), labels=levels(group))
        if(seed) set.seed(seed)
        while(i <= B) {
            S.1 <- sample(1:dim(x.1)[1], replace = TRUE)
            S.2 <- sample(1:dim(x.2)[1], replace = TRUE)
            X.1 <- x.1[S.1, , drop=FALSE]
            X.2 <- x.2[S.2, , drop=FALSE]
            X <- rbind(X.1, X.2)	# BSS
            G <- group[c(S.1, dim(x.1)[1] + S.2)]
            if(isSingular(Cov(t(t(X) - c(rep(colMeans(X.1),length(S.1)), rep(colMeans(X.2), length(S.2))))))) next			# next loop when singular W
            CR.B <- lda.c(X, G, sub = sub)
            coef.B[, i] <- CR.B$d.function
            i <- i + 1
        }
    } else {
        if(seed) set.seed(seed)
        while(i <= B) {
            S <- sample(1:N, replace=TRUE)
            X <- x[S, , drop=FALSE]	# BSS
            G <- group[S]
            if(any(table(G) == 0)) next
            if(isSingular(Cov(X - t(matrix(unlist(by(X, G, colMeans)),ncol = 2)[,as.numeric(G)])))) next			# next loop when singular W
            CR.B <- lda.c(X, G, sub = sub)
            coef.B[, i] <- CR.B$d.function
            i <- i + 1
        }
    }
    coef.B <- t(coef.B)
    coef.vcov <- var(coef.B)
    sum.coef <- summary(coef.B)
    CI95 <- apply(coef.B, 2, sort)[c(B * 0.025, B * 0.975),]
    rownames(CI95) <- c("lower", "upper")
    if(B >= 200) {
       CI99 <- apply(coef.B, 2, sort)[c(B * 0.005, B * 0.995),]
       rownames(CI99) <- c("lower", "upper")
    } else CI99 <- NULL
    list(vcov = coef.vcov, summary = sum.coef, CI95 = CI95, CI99 = CI99)
}


## PC version of Coef.B
Coef.pc.B <- function(x, group, sub = 1:ncol(x), spc = -1, B, fixed = FALSE, 
                      seed = FALSE) {
    require(rrcov)
    N <- nrow(x)
    k <- length(levels(group))
    B <- B
    coef.B <- matrix(NA, ncol = B, nrow = ncol(x[sub]) + 1)
    rownames(coef.B) <- rownames(lda.pc(x, group, sub = sub, spc = spc)$d.function.or)
    i <- 1
    if(fixed) {				# Use fixed number of subsamples for each group
        x.1 <- x[as.numeric(group) == 1, , drop = FALSE]
        x.2 <- x[as.numeric(group) == 2, , drop = FALSE]
        x <- rbind(x.1, x.2)
        group <- factor(c(group[as.numeric(group) == 1],
                          group[as.numeric(group)==2]), labels = levels(group))
        if(seed) set.seed(seed)
        while(i <= B) {
            S.1 <- sample(1:dim(x.1)[1], replace=TRUE)
            S.2 <- sample(1:dim(x.2)[1], replace=TRUE)
            X.1 <- x.1[S.1, , drop=FALSE]
            X.2 <- x.2[S.2, , drop=FALSE]
            X <- rbind(X.1, X.2)	# BSS
            G <- group[c(S.1, dim(x.1)[1] + S.2)]
            if(isSingular(Cov(t(t(X)-c(rep(colMeans(X.1),length(S.1)), rep(colMeans(X.2),length(S.2))))))) next			# next loop when singular W
            CR.B <- lda.pc(X, G, sub = sub, spc = spc)
            coef.B[, i] <- CR.B$d.function.or
            i <- i + 1
        }
    } else {
        if(seed) set.seed(seed)
        while(i <= B) {
            S <- sample(1:N, replace=TRUE)
            X <- x[S, , drop=FALSE]	# BSS
            G <- group[S]
            if(any(table(G) == 0)) next
            if(isSingular(Cov(X - t(matrix(unlist(by(X, G, colMeans)), ncol = 2)[, as.numeric(G)])))) next			# next loop when singular W
            CR.B <- lda.pc(X, G, sub = sub, spc = spc)
            coef.B[, i] <- CR.B$d.function.or
            i <- i + 1
        }
    }
    coef.B <- t(coef.B)
    coef.vcov <- var(coef.B)
    sum.coef <- summary(coef.B)
    CI95 <- apply(coef.B, 2, sort)[c(B * 0.025, B * 0.975),]
    rownames(CI95) <- c("lower", "upper")
    if(B >= 200) {
       CI99 <- apply(coef.B, 2, sort)[c(B * 0.005, B * 0.995),]
       rownames(CI99) <- c("lower", "upper")
    } else CI99 <- NULL
    list(vcov = coef.vcov, summary = sum.coef, CI95 = CI95, CI99 = CI99)
}
