#' Performs forward projection PLS with 
#' the cindex of given conditional variable index.
#' @param x the original data matrices of dimension $n x p$.
#' @param y the original response variable of dimension $n x 1$.
#' @param cindex the index of given conditional variable 'cindex'.
#' @param K the number of forward projection.
#' @param ratio the number of final selected variables.

fppls_beta <- function(x, y, cindex, K, ratio) {
  library(pls)
  n <- nrow(x)
  p <- ncol(x)
  colnames(x) <- paste0("x", 1:p)
  ox <- x
  oy <- y
  x <- (diag(n) - x[, cindex] %*% solve(t(x[, cindex]) %*% x[, cindex]) %*% t(x[, cindex])) %*% x[, -cindex]
  y <- (diag(n) - x[, cindex] %*% solve(t(x[, cindex]) %*% x[, cindex]) %*% t(x[, cindex])) %*% y

  A <- c()
  rmsecv <- matrix(0, length(ratio), K)
  xrank <- list()
  x1 <- list()
  y1 <- list()

  for (k in 1:K) {
    newdata <- data.frame(cbind(x, y))
    p1 <- ncol(x)
    names(newdata)[ncol(newdata)] <- "y"
    comp <- min(ncol(newdata) - 1, 10)
    cpls <- plsr(y ~ ., comp, data = newdata, validation = "CV")
    r <- rep(0, comp)
    for (j in 1:comp) {
      r[j] <- 1 - sum((as.numeric(cpls$validation$pred[, , j]) - y)^2) / sum((y - mean(y))^2)
    }
    ncomp <- which.max(r)
    cpls1 <- plsr(y ~ ., ncomp, data = newdata)
    cpls_beta <- coefficients(cpls1)
    cpls_beta_rank <- order(abs(cpls_beta), decreasing = TRUE)
    cindex1 <- cpls_beta_rank[1]

    cindexname <- colnames(x)[cindex1]
    a <- as.numeric(substr(cindexname, 2, nchar(cindexname)))
    x2 <- (diag(n) - x[, cindex1] %*% solve(t(x[, cindex1]) %*% x[, cindex1]) %*% t(x[, cindex1])) %*% x[, -cindex1]
    y <- (diag(n) - x[, cindex1] %*% solve(t(x[, cindex1]) %*% x[, cindex1]) %*% t(x[, cindex1])) %*% y
    x <- x2
    A <- unique(c(A, a))
    x1[[k]] <- x
    y1[[k]] <- y
  }

  for (k in 1:K) {
    for (i in 1:length(ratio)) {
      newdata <- data.frame(cbind(x1[[k]], y1[[k]]))
      p1 <- ncol(x1[[k]])
      names(newdata)[ncol(newdata)] <- "y"
      comp <- min(ncol(newdata) - 1, 10)
      cpls <- plsr(y ~ ., comp, data = newdata, validation = "CV")
      r <- rep(0, comp)
      for (j in 1:comp) {
        r[j] <- 1 - sum((as.numeric(cpls$validation$pred[, , j]) - y1[[k]])^2) / sum((y1[[k]] - mean(y1[[k]]))^2)
      }
      ncomp <- which.max(r)
      cpls1 <- plsr(y ~ ., ncomp, data = newdata)
      cpls_beta <- coefficients(cpls1)
      cpls_beta_rank <- order(abs(cpls_beta), decreasing = TRUE)
      cindexname <- colnames(x1[[k]])[cpls_beta_rank]
      a <- as.numeric(substr(cindexname, 2, nchar(cindexname)))
      xrank[[(k - 1) * length(ratio) + i]] <- c(cindex, A[1:k], a)
      value <- rep(0, 10)
      index <- rep_len(1L:10, n)
      for (nfolds in 1L:10) {
        xtrain <- ox[index != nfolds, ]
        ytrain <- oy[index != nfolds]
        xtest <- ox[index == nfolds, ]
        ytest <- oy[index == nfolds]

        xrank1 <- xrank[[(k - 1) * length(ratio) + i]]
        varindex <- unlist(xrank1[1:ratio[i]])

        ncorn <- data.frame(cbind(xtrain[, varindex], ytrain))

        names(ncorn)[ncol(ncorn)] <- "y"
        comp <- min(ratio[i], 10)

        corn_pls <- plsr(y ~ ., comp, data = ncorn, validation = "CV")
        r <- c()
        for (j in 1:comp) {
          r[j] <- 1 - sum((as.numeric(corn_pls$validation$pred[, , j]) - ytrain)^2) / sum((ytrain - mean(ytrain))^2)
        }
        ncomp <- which.max(r)
        corn_pls1 <- plsr(y ~ ., ncomp, data = ncorn)
        cpls_beta <- coefficients(corn_pls1)
        ncorn_test <- data.frame(xtest[, varindex])
        corn_pre <- predict(corn_pls1, ncorn_test)

        value[nfolds] <- sum((as.numeric(corn_pre) - ytest)^2)
      }

      rmsecv[i, k] <- mean(value)
    }
  }

  maxratio <- which.min(rmsecv) %% length(ratio)
  if (maxratio == 0) {
    maxratio <- length(ratio)
  }
  cond_var <- ceiling(which.min(rmsecv) / length(ratio))

  optimal_ratio <- ratio[maxratio]
  xrank1 <- xrank[[which.min(rmsecv)]]
  optimal_var <- c(xrank1[1:optimal_ratio])

  ncorn <- data.frame(cbind(ox[, optimal_var], oy))
  names(ncorn)[ncol(ncorn)] <- "y"
  comp <- min(optimal_ratio, 10)
  corn_pls <- plsr(y ~ ., comp, data = ncorn, validation = "CV")
  r <- c()
  for (j in 1:comp) {
    r[j] <- 1 - sum((as.numeric(corn_pls$validation$pred[, , j]) - oy)^2) / sum((oy - mean(oy))^2)
  }
  ncomp <- which.max(r)
  corn_pls1 <- plsr(y ~ ., ncomp, data = ncorn)

  list(rank = xrank1, var = optimal_var, model = corn_pls1, ncomp = ncomp)
}
