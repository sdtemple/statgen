'
File: inbr&knshp_est
Author: Seth Temple
Description: R functions to estimate inbreeding and coancestry in population genetics
'

### faster matrix algebra approaches ###

MJ <- function(d0, d1){
  '
  Calculates matching for paired alleles in an individual using dosage matrices
  d0 <- major allele dosage matrix
  d1 <- minor allele dosage matrix
  Return dosage matching matrix
  '
  return((d0 *  (d0 - 1) + d1 *  (d1 - 1)) / 2)
}

MJJ <- function(d, s){
  '
  Calculates SNP matching for alleles between individuals using a dosage matrix
  d0 <- allele dosage matrix
  s <- SNP id
  Return dosage matching matrix
  '
  return(((d[s,] - 1) %*% t(d[s,] - 1) + 1) / 2)
}

MSL <- function(d){
  '
  Calculates average SNP matching between individuals using a dosage matching matrix
  d - dosage matching matrix
  Return an average
  '
  n <- nrow(d)
  return(sum(d * upper.tri(d)) / n / (n - 1) * 2)
}

fBETA <- function(d0, d1){
  '
  Calculates inbreeding and coancestry coefficients based on matrix algebra
  d0 - major allele dosage matrix
  d1 - minor allele dosage matrix
  Return symmetric matrix
  '
  
  # initialize some variables
  nvariant <- nrow(d0)
  nsample <- ncol(d0)
  m <- matrix(0, nrow = nsample, ncol = nsample)
  msl <- rep(NA, nvariant)
  
  # compute and place inbreeding values
  mj <- apply(MJ(d0, d1), 2, sum)
  m <- m + diag(mj)
  
  # compute coancestry values and apply weights
  for(i in 1:nvariant){
    mjj <- MJJ(d0, i)
    msl[i] <- MSL(mjj)
    m <- mjj * upper.tri(mjj) - msl[i] + m
    m <- m * upper.tri(m, diag = TRUE)
  }
  
  # normalize
  wt <- sum(1 - msl)
  m <- m / wt
  
  m <- m + t(m) * lower.tri(m)
  return(m)
}

# does not match PLINK GCTA run
# likely due to reference population adjustment
kGCTA <- function(d1, p){
  '
  Calculates inbreeding and coancestry coefficients based on matrix algebra
  d1 - major allele dosage matrix
  p - vector of sample alleles
  Return symmetric matrix
  '
  
  # initialize some variables
  nvariant <- length(p)
  nsample <- ncol(d1)
  m <- matrix(0, nrow = nsample, ncol = nsample)
  msl <- rep(NA, nvariant)
  
  # initialize scale and shift dosage
  scale <- 1 / p / (1 - p) / 2
  x <- d1 - 2 * p
  
  for(i in 1:nvariant){
    k <- scale[i] * (x[i,] %*% t(x[i,]))
    msl[i] <- sum(k * upper.tri(k)) / nsample / (nsample - 1)
    m <- k * upper.tri(k) - msl[i] + m
    m <- m * upper.tri(m)
  }
  
  # normalize
  wt <- sum(1 - msl)
  m <- m / wt
  
  # impose symmetry
  m <- m + t(m) * lower.tri(m)
  return(m)
}

### fast iterative approaches ###

fhom <- function(a1, a2, p){
  '
  Calculates a measure of inbreeding
  a1 - a matrix of alleles
  a2 - a matrix of alleles
  p - vector of sample allele frequencies
  Return a vector of inbreeding coefficients
  '
  nsample <- ncol(a1)
  nvariant <- length(p)
  f <- rep(NA, length.out = nsample)
  
  s <- 0
  for(i in 1:nvariant){
    phat <- p[i]
    s <- 2 * phat * (1 - phat) + s
  }
  
  for(j in 1:nsample){
    h <- 0
    for(k in 1:nvariant){
      x <- a1[k,j] + a2[k,j]
      h <- x * (2 - x) + h
    }
    f[j] <- 1 - h / s
  }
  return(f)
}

fUNI <- function(d1, p){
  '
  Calculates a measure of inbreeding
  d1 - major allele dosage matrix
  p - vector of sample allele frequencies
  Return a vector of inbreeding coefficients
  '
  nsample <- ncol(d1)
  nvariant <- length(p)
  f <- rep(NA, length.out = nsample)
  
  for(j in 1:nsample){
    val <- 0 
    for(k in 1:nvariant){
      x <- d1[k,j]
      phat <- p[k]
      val <- (x ^ 2 - (1 + 2 * phat) * x + 2 * phat ^ 2) / (2 * phat * (1 - phat)) + val
    }
    f[j] <- val / nvariant
  }
  return(f)
}

fGCTA <- function(d1, p){
  '
  Calculates a measure of inbreeding
  d1 - major allele dosage matrix
  p - vector of sample allele frequencies
  Return a vector of inbreeding coefficients
  '
  nsample <- ncol(d1)
  nvariant <- length(p)
  f <- rep(NA, length.out = nsample)
  
  for(j in 1:nsample){
    val <- 0
    for(k in 1:nvariant){
      x <- d1[k,j]
      phat <- p[k]
      val <- ((x - 2 * phat) ^ 2) / (2 * phat * (1-phat)) + val
    }
    f[j] <- val / nvariant - 1
  }
  return(f) 
}

### slower iterative approaches ###

# iterative approach
funi <- function(a1, a2, p){
  '
  Calculates a measure of inbreeding
  a1 - matrix of alleles
  a2 - matrix of alleles
  p - vector of sample allele frequencies
  Return a vector of inbreeding coefficients
  '
  nsample <- ncol(a1)
  nvariant <- length(p)
  f <- rep(NA, length.out = nsample)
  
  for(j in 1:nsample){
    val <- 0 
    for(k in 1:nvariant){
      x <- a1[k,j] + a2[k,j]
      phat <- p[k]
      val <- (x ^ 2 - (1 + 2 * phat) * x + 2 * phat ^ 2) / (2 * phat * (1 - phat)) + val
    }
    f[j] <- val / nvariant
  }
  return(f)
}

# iterative approach
fgcta <- function(a1, a2, p){
  '
  Calculates a measure of inbreeding
  a1 - matrix of alleles
  a2 - matrix of alleles
  p - vector of sample allele frequencies
  Return a vector of inbreeding coefficients
  '
  nsample <- ncol(a1)
  nvariant <- length(p)
  f <- rep(NA, length.out = nsample)
  
  for(j in 1:nsample){
    val <- 0
    for(k in 1:nvariant){
      x <- a1[k,j] + a2[k,j]
      phat <- p[k]
      val <- ((x - 2 * phat) ^ 2) / (2 * phat * (1-phat)) + val
    }
    f[j] <- val / nvariant - 1
  }
  return(f)
}


mj <- function(a1, a2){
  '
  Calculates matching for paired alleles in an individual
  a1 - first allele
  a2 - second allele
  '
  x0 <- -((a1 - 1) + (a2 - 1))
  x1 <- a1 + a2
  return(((x0 * (x0 - 1)) + (x1 * (x1 - 1))) / 2)
}

mjj <- function(a11, a12, a21, a22){
  '
  Calculates allele matching for two individuals
  a11 - first allele in first individual
  a12 - second allele in first individual
  a21 - first allele in second individual
  a22 - second allele in second individual
  '
  x10 <- -((a11 - 1) + (a12 - 1))
  x11 <- a11 + a12
  x20 <- -((a21 - 1) + (a22 - 1))
  x21 <- a21 + a22
  return((x10 * x20 + x11 * x21) / 4)
}

msl <- function(a1, a2, l){
  '
  Calculates average between individual matching in a population
  a1 - allele matrix
  a2 - allele matrix
  l - a SNP id in the allele matrices
  '
  nvariant <- nrow(a1)
  nsample <- ncol(a1)
  m <- 0
  for(i in 1:nsample){
    for(j in 1:nsample){
      if(i != j){
        m <- mjj(a1[l,i], a2[l,i], a1[l,j], a2[l,j]) + m
      }
    }
  }
  return(m / (nsample * (nsample - 1)))
}

# high runtime complexity
# not tested
fbeta <- function(a1, a2){
  '
  Calculates a measure of inbreeding based on allele matching
  a1 - allele matrix
  a2 - allele matrix
  '
  nsample <- ncol(a1)
  nvariant <- nrow(a1)
  ms <- rep(NA, length.out = nvariant)
  for(l in 1:nvariant){
    ms[l] <- msl(a1, a2, l)
  }
  divisor <- sum(1 - ms)
  
  b <- rep(NA, length.out = nsample)
  for(j in 1:nsample){
    bj <- 0
    for(l in 1:nvariant){
      bj <- (mj(a1[l,j], a2[l,j]) - ms[l]) + bj
    }
    b[j] <- bj / divisor
  }
  return(b)
}

# high runtime complexity
# not tested
kbeta <- function(a1, a2){
  '
  Calculates a measure of coancestry based on allele matching
  a1 - allele matrix
  a2 - allele matrix
  '
  nsample <- ncol(a1)
  nvariant <- nrow(a1)
  ms <- rep(NA, length.out = nvariant)
  for(l in 1:nvariant){
    ms[l] <- msl(a1, a2, l)
  }
  divisor <- sum(1 - ms)
  
  b <- matrix(NA, nrow = nsample, ncol = nsample)
  for(i in 1:nsample){
    for(j in 1:nsample){
      bj <- 0
      for(l in 1:nvariant){
        bj <- (mjj(a1[l,i], a2[l,i], a1[l,j], a2[l,j]) - ms[l]) + bj
      }
      b[i,j] <- bj / divisor
      b[j,i] <- b[i,j]
    }
  }
  return(b)
}
