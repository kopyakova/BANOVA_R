results.BANOVA.mlvNormal <- function(fit_beta, dep_var_names, dMatrice){
  combine.tables <- function(list_with_tables, list_names_first = T, list_elements_length = 1,
                             anova = F){
    combine <- function(list_with_tables, element_name = "NA"){
      row_names <- rownames(list_with_tables[[1]])
      col_names <- colnames(list_with_tables[[1]])
      
      if(is.null(row_names)){
        if (is.null(dim(list_with_tables[[1]]))){
          row_names = " "
          col_names = " "
        } else {
          dimensions <- dim(list_with_tables[[1]])
          row_names = c(1:dimensions[1])
          col_names = c(1:dimensions[2])
        }
      }
      
      n_row <- length(row_names)
      n_col <- length(col_names)
      
      if(list_names_first){
        first_vec  <- list_names
        second_vec <- row_names
      }else{
        first_vec  <- row_names
        second_vec <- list_names
      }
      new_row_names <- c()
      for (element in second_vec){
        new_row_names <- c(new_row_names, format(paste0(first_vec, "  ", element)))
      }
      
      combined_table <- matrix(NA, nrow = length(row_names)*list_length, ncol = n_col)
      if (list_names_first){
        for (i in 1:list_length){
          dv_name <- list_names[i]
          combined_table[(1+n_row*(i-1)):(n_row+n_row*(i-1)),] <- as.matrix(list_with_tables[[dv_name]])
        }
      } else {
        for (i in 1:list_length){
          dv_name <- list_names[i]
          matrix  <- as.matrix(list_with_tables[[dv_name]])
          for (j in 1:n_row){
            combined_table[1+(j-1)*list_length+(i-1),] <- matrix[j,]
          }
        }
      }
      
      if (element_name == "pass_ind"){
        combined_table <- all(combined_table == TRUE) 
      } else if (element_name == "row_indices"){
        colnames(combined_table) <- c(col_names)
        rownames(combined_table) <- new_row_names
      } else {
        combined_table <- as.data.frame(combined_table)
        colnames(combined_table) <- c(col_names)
        rownames(combined_table) <- new_row_names
      }
   
      return(combined_table)
    }
    
    list_names  <- names(list_with_tables)
    list_length <- length(list_with_tables)
    
    if (list_elements_length == 1){
      combined_tables <- combine(list_with_tables)
    } else {
      element_names <- names(list_with_tables[[1]])
      combined_tables <- list()
      for (element_name in element_names){
        temp_list <- list()
        for (name in list_names)
          temp_list[[name]] <- list_with_tables[[name]][[element_name]]
        combined_tables[[element_name]] <- combine(temp_list, element_name)
      }
      if (anova){
        class(combined_tables) <- "ancova.effect"
      }
    }
    return(combined_tables)
  }
  combine.samples <- function(list_with_samples){
    list_names  <- names(list_with_samples)
    list_length <- length(list_with_samples)
    n_row <- nrow(list_with_samples[[1]])
    n_col <- ncol(list_with_samples[[1]])
    combined_samples <- matrix(NA, nrow = n_row, ncol = n_col*list_length)
    new_colnames     <- c()
    init_val <- 1
    term_val <- n_col
    for (i in 1:list_length){
      samples <- list_with_samples[[i]]
      combined_samples[, init_val:term_val] <- samples 
      new_colnames <- c(new_colnames, paste0(colnames(samples), "_", list_names[i]))
      init_val <- term_val + 1
      term_val <- n_col*(i+1)
    }
    colnames(combined_samples) <- new_colnames
    return(combined_samples)
    
  }
  
  X <- dMatrice$X
  Z <- dMatrice$Z
  X_names <- colnames(X)
  Z_names <- colnames(Z)
  
  ##### Extract first and second level coefficients ##### 
  # Identify dimentsions
  beta1_dim <- dim(fit_beta$beta1)
  beta2_dim <- dim(fit_beta$beta2)
  L = beta2_dim[2] # number of dependent variables
  M = beta1_dim[2] # number of unique subjects
  J = beta1_dim[4] # number of subject-level effects (1st level)
  K = beta2_dim[3] # number of population-level effects (2nd level)
  n_iter =  beta1_dim[1] # number MCMC iterations / number of samples
  # Prepare names of the coefficients 
  beta1_names <- c()
  for (i in 1:J){
    for (j in 1:M){
      beta1_names <- c(beta1_names, paste("beta1_",i,"_",j, sep = ""))
    }
  }
  beta2_names <- c()
  for (i in 1:J){
    for (j in 1:K){
      beta2_names <- c(beta2_names, paste("beta2_",i,"_",j, sep = ""))
    }
  }
  ##### 
  #Extract R2 
  R2 = NULL
  if (!is.null(fit_beta$r_2)){
    R2 <- colMeans(fit_beta$r_2)
    R2 <- round(R2, 4)
  }
  
  #Extract Sigma (covariane matrix) and variances of y variables
  tau_ySq = NULL
  Sigma   = NULL
  Sigma_samples <- fit_beta$Sigma
  if (!is.null(Sigma_samples)){
    Sigma <- colMeans(Sigma_samples) # Covariane matrix
    tau_ySq <- diag(Sigma)           # Vector with variances for dep vars
  }
  
  #Extract Omega (correlation matrix) and variances of y variables
  Omega = NULL
  Omega_samples <- fit_beta$Omega
  if (!is.null(Omega_samples)){
    Omega <- colMeans(Omega_samples)
  }
  
  ##### Prepare results for each of L dependent variables ##### 
  samples_l1.list    <- list()
  samples_l2.list    <- list()
  anova.tables.list  <- list()
  coef.tables.list   <- list()
  pvalue.tables.list <- list()
  conv.list          <- list()
  cat('Constructing ANOVA/ANCOVA tables...\n')
  for (l in 1:L){
    dep_var_name <- dep_var_names[l] #column name of the l_th dependent variable
    
    # Reformat level one and level two coefficient samples sample of the l_th dependent variable
    samples_l1_param <- array(0, dim = c(n_iter, J*M), dimnames = list(NULL, beta1_names))
    for (i in 1:J){
      for (j in 1:M){
        samples_l1_param[, (i-1) * M + j] <- fit_beta$beta1[, j, l, i]
      }
    }
    samples_l2_param <- array(0, dim = c(n_iter, K*J), dimnames = list(NULL, beta2_names))
    for (i in 1:J){
      for (j in 1:K){
        samples_l2_param[, (i-1) * K + j] <- fit_beta$beta2[, l, j, i]
      }
    }
    samples_l1.list[[dep_var_name]] <- samples_l1_param
    samples_l2.list[[dep_var_name]] <- samples_l2_param
    # Result tables for of the l_th dependent variable
    anova.tables.list[[dep_var_name]] <- table.ANCOVA(samples_l1_param, X, Z, 
                                                      samples_l2_param, l1_error = tau_ySq[l])
    coef.tables.list[[dep_var_name]]  <- table.coefficients(samples_l2_param, beta2_names, 
                                                            X_names, Z_names, attr(X, 'assign') + 1, 
                                                            attr(Z, 'assign') + 1, 
                                                            samples_cutp_param = array(dim = 0))
    pvalue.tables.list[[dep_var_name]] <- table.pvalue(coef.tables.list[[dep_var_name]]$coeff_table, 
                                                       coef.tables.list[[dep_var_name]]$row_indices, 
                                                       l1_names = attr(X, 'varNames'), 
                                                       l2_names = attr(Z, 'varNames'))
    conv.list[[dep_var_name]]          <- conv.geweke.heidel(samples_l2_param, X_names, Z_names)
    class(conv.list[[dep_var_name]]) <- 'conv.diag'
  }
  ##### Combine the results #####
  combined.anova  <- combine.tables(anova.tables.list, T, 2, anova = T)
  combined.coef   <- combine.tables(coef.tables.list, T, 3)
  combined.pvalue <- combine.tables(pvalue.tables.list, T, 1)
  combined.conv   <- combine.tables(conv.list, T, 2)
  
  combined.samples.l1 <- combine.samples(samples_l1.list)
  combined.samples.l2 <- combine.samples(samples_l2.list)
  cat('Done.\n')
  return(list(combined.anova = combined.anova,
              combined.coef = combined.coef, 
              combined.pvalue = combined.pvalue,
              combined.conv = combined.conv,
              combined.samples.l1 = combined.samples.l1,
              combined.samples.l2 = combined.samples.l2, 
              samples_l1.list = samples_l1.list,
              samples_l2.list = samples_l2.list,
              covariance.matrix = Sigma,
              correlation.matrix = Omega,
              R2 = R2,
              tau_ySq = tau_ySq,
              anova.tables.list = anova.tables.list,
              coef.tables.list = coef.tables.list,
              pvalue.tables.list = pvalue.tables.list,
              conv.list = conv.list,
              num_depenent_variables = length(dep_var_names),
              names_of_dependent_variables = dep_var_names))
}

source('~/BANOVA_R/R/table.ANCOVA.R', echo=F)
source('~/BANOVA_R/R/table.coefficients.R', echo=F)
source('~/BANOVA_R/R/table.pvalue.R', echo=F)
source('~/BANOVA_R/R/pValues.R', echo=F)
source('~/BANOVA_R/R/conv.geweke.heidel.R', echo=TRUE)
source('~/BANOVA_R/R/effect.matrix.factor.R', echo=TRUE)


# > class(anova.tables.list[[dep_var_name]]$ancova_table)
# [1] "data.frame"
# > class(anova.tables.list[[dep_var_name]]$effect_table)
# [1] "data.frame"

# class(combined.anova$ancova_table)
# class(combined.anova$effect_table)
# > class(coef.tables.list[[dep_var_name]]$row_indices)
# [1] "matrix"
# > class(coef.tables.list[[dep_var_name]]$full_table)
# [1] "data.frame"
# > class(coef.tables.list[[dep_var_name]]$coeff_table)
# [1] "data.frame"
# > class(pvalue.tables.list[[dep_var_name]])
# [1] "data.frame"
# > class(conv.list[[dep_var_name]]$sol_geweke)
# [1] "data.frame"
# > class(conv.list[[dep_var_name]]$sol_heidel)
# [1] "data.frame"
# > class(conv.list[[dep_var_name]]$pass_ind)
# [1] "logical"
