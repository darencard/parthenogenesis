`ir` <-
  function(genotypes) {
    
    genotypes <- as.matrix(genotypes)
    
    individuals <- nrow(genotypes)
    loci <- ncol(genotypes) / 2
    ir <- array(NA, dim=c(individuals, 1))
    frequencies <- array()
    
    for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      frequencies[l] <- list(table(genotypes[, g:h]))
    }
    
    for (i in 1:individuals) {
      
      H <- 0
      N <- 0
      f <- 0
      
      for (l in 1:loci) {
        
        g <- 2 * l - 1
        h <- 2 * l
        
        if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
          N <- N + 1
          if (genotypes[i, g] == genotypes[i, h]) {
            H <- H + 1
            c <- as.character(genotypes[i, g])
            f <- f + (2 * frequencies[[l]][[c]] - 2) / (sum(frequencies[[l]]) - 2)
          }
          else {
            c <- as.character(genotypes[i, g])
            f <- f + (frequencies[[l]][[c]] - 1) / (sum(frequencies[[l]]) - 2)
            c <- as.character(genotypes[i, h])
            f <- f + (frequencies[[l]][[c]] - 1) / (sum(frequencies[[l]]) - 2)
          }
        }
      }
      
      ir[i] <- (2 * H - f) / (2 * N - f)
      
    }
    
    ir
    
  }


`sh` <-
  function(genotypes) {
    
    genotypes <- as.matrix(genotypes)
    
    individuals <- nrow(genotypes)
    loci <- ncol(genotypes) / 2
    sh <- array(NA, dim=c(individuals, 1))
    heterozygosity <- array()
    
    for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      freq <- table((genotypes[, g] == genotypes[, h]))
      if (nrow(freq) == 2) {
        heterozygosity[l] <- sweep(freq, 1, sum(freq), "/")[["FALSE"]]
      }
      else {
        if (labels(freq)[[1]][1] == TRUE) {
          heterozygosity[l] <- 1 - sweep(freq, 1, sum(freq), "/")[["TRUE"]]
        }
        else {
          heterozygosity[l] <- sweep(freq, 1, sum(freq), "/")[["FALSE"]]
        }
      }
      
    }
    
    for (i in 1:individuals) {
      
      H <- 0
      N <- 0
      mh <- 0
      
      for (l in 1:loci) {
        
        g <- 2 * l - 1
        h <- 2 * l
        
        if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
          mh <- mh + heterozygosity[l] 
          N <- N + 1
          if (genotypes[i, g] != genotypes[i, h]) {
            H <- H + 1
          }
        }
      }
      
      sh[i] <- (H / N) / (mh / N)
      
    }
    
    sh
    
  }


`hl` <-
  function(genotypes) {
    
    genotypes <- as.matrix(genotypes)
    
    individuals <- nrow(genotypes)
    loci <- ncol(genotypes) / 2
    hl <- array(NA, dim=c(individuals, 1))
    E <- array(loci)
    frequencies <- array(loci)
    
    for (l in 1:loci) {
      E[l] <- 1
      g <- 2 * l - 1
      h <- 2 * l
      frequencies[l] <- list(table(genotypes[, g:h]))
      E[l] <- 1 - sum((frequencies[[l]] / sum(frequencies[[l]]))^2)
    }
    
    for (i in 1:individuals) {
      
      sum.Eh <- 0
      sum.Ej <- 0
      
      for (l in 1:loci) {
        
        g <- 2 * l - 1
        h <- 2 * l
        
        if ((!is.na(genotypes[i, g])) && (!is.na(genotypes[i, h]))) {
          if (genotypes[i, g] == genotypes[i, h]) {
            sum.Eh <- sum.Eh + E[l]
          }
          else {
            sum.Ej <- sum.Ej + E[l]
          }
        }
      }
      
      hl[i] <- sum.Eh / (sum.Eh + sum.Ej)
      
    }
    
    hl
    
  }


`chkdata` <-
  function(genotypes) {
    
    genotypes <- as.matrix(genotypes)
    
    if (ncol(genotypes) %% 2 == 1) {
      stop("Odd number of columns in the input.")
    }
    
    individuals <- nrow(genotypes)
    loci <- ncol(genotypes) / 2
    n_alleles <- array(loci)
    
    for (l in 1:loci) {
      g <- 2 * l - 1
      h <- 2 * l
      frequencies <- list(table(genotypes[, g:h]))
      n_alleles[l] <- length(frequencies[[1]])
    }
    
    if (sum(n_alleles < 2) > 0) {
      loci_list <- ""
      for (i in 1:loci) {
        if (n_alleles[i] < 2) {
          if (loci_list == "") {
            loci_list <- i
          }
          else {
            loci_list <- paste(loci_list, i, sep = ", ")
          }
        }
      }
      stop(paste("Only one allele in loci:", loci_list))
    }
    
    write(n_alleles, file="number_of_alleles.txt", ncolumns=1)
    
    k <- 0
    ind_loci_list <- "\n"
    
    for (i in 1:individuals) {
      
      for (l in 1:loci) {
        
        g <- 2 * l - 1
        h <- 2 * l
        
        if (xor(!is.na(genotypes[i, g]), !is.na(genotypes[i, h]))) {
          ind_loci_list <- paste(ind_loci_list, i, l, "\n")
          k <- k + 1
        }
      }
    }
    
    if (k > 0) {
      warning(paste("One or more individuals are missing one allele in one or more loci:\nIndividual  Locus", ind_loci_list))
    }
  }
    

`mlh` <-
  function(data) {
    
    # Read data:
    #g <- read.table(input, na.strings = na.string)
    g <- data
    
    # Check data:
    chkdata(g)
    
    # Reserve array for results:
    h <- as.data.frame(array(dim = c(nrow(g), 4)), stringsAsFactors = FALSE)
    
    # Calculate heterozygosity estimates:
    h[, 1] <- g[, 1]
    h[, 2] <- ir(g[, 2:ncol(g)])
    h[, 3] <- sh(g[, 2:ncol(g)])
    h[, 4] <- hl(g[, 2:ncol(g)])
    colnames(h) <- c("ID", "IR", "SH", "HL")
    
    # Write results to the output file:
    #write.table(format(h, digits=n.digits), file = output, col.names = c("ID", "IR", "SH", "HL"), row.names = FALSE, quote = FALSE, sep="\t")
    
    h
    
  }

boot_het <- function(data, method, reps) {
  out <- data.frame()
  nloci <- (ncol(data)-1)/2
  for (rep in 1:reps) {
    print(paste0("Bootstrap replicate #", rep))
    j <- numeric(2 * nloci)
    random <- sample(seq(1, nloci, 2), nloci, replace=TRUE)
    for (i in 1:length(random)) {
      j[2*i-1] <- random[i]
      j[2*i] <- random[i]+1
    }
    if (method == "sh" || method == "ir" || method == "hl") {
      if (method == "sh") {
        out <- rbind(out, c(rep, t(sh(data[,-1][,j]))))
      } else if (method == "ir") {
        out <- rbind(out, c(rep, t(ir(data[,-1][,j]))))
      } else {
        out <- rbind(out, c(rep, t(hl(data[,-1][,j]))))
      }
    }
  }
  names(out) <- c("rep", as.character(data[,1]))
  return(out)
}
