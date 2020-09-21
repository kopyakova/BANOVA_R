BANOVA.multi.mediation <- function(sol_1, sol_2, xvar, mediators, individual = F){
  adapt.design.matrix <- function(dMatrice, mediator){
    d_temp <- dMatrice
    #X
    attributes(d_temp$X)$varValues[[1]] <- attributes(dMatrice$X)$varValues[[1]][, mediator]
    attr(attributes(d_temp$X)$varValues[[1]], "var_names") <- mediator
    #y
    d_temp$y <- dMatrice$y[,mediator]
    attr(d_temp$y, "names") <- attr(dMatrice$y, "dimnames")[[1]] 
    return(d_temp)
  }
  adapt.mf1 <- function(mf, mediator){
    temp_mf <- mf
    mediator_of_interest <- mf[,1][,mediator]
    
    #data frame
    temp_mf[, 1] <- mediator_of_interest
    colnames(temp_mf)[1] <- mediator
    
    #terms attribute
    rownames(attr(attr(temp_mf, "terms"),'factors'))[1] <- mediator
    attr(attr(temp_mf, 'terms'),'dataClasses')[1] <- "numeric"
    names(attr(attr(temp_mf, 'terms'),'dataClasses'))[1] <- mediator
    return(temp_mf)
  }

  #Check the class of the model with multiple dependent variables
  if(sol_2$model_name != "BANOVA.multiNormal") 
    stop('The mediator must follow the multivariate Normal distribution, use BANOVA multiNormal models instead.')
  
  #Adapt the output of BANOVA.multiNormal to fit BANOVA.mediation
  temp_solution <- list()
  temp_solution$samples_cutp_param <- sol_2$samples_cutp_param
  temp_solution$data               <- sol_2$data
  temp_solution$num_trials         <- sol_2$num_trials
  temp_solution$model_code         <- sol_2$model_code
  temp_solution$single_level       <- sol_2$single_level
  temp_solution$stan_fit           <- sol_2$stan_fit
  temp_solution$contrast           <- sol_2$contrast
  temp_solution$new_id             <- sol_2$new_id
  temp_solution$old_id             <- sol_2$old_id
  temp_solution$call               <- sol_2$old_id
  temp_solution$model_name         <- 'BANOVA.Normal'
  temp_solution$samples_l2_sigma_param <- sol_2$samples_l2_sigma_param
  class(temp_solution) <- "BANOVA"
  names(sol_2$R2)      <- names(sol_2$anova.tables.list)
  names(sol_2$tau_ySq) <- names(sol_2$anova.tables.list)
  
  results <-list()
  for (mediator in mediators){
    temp_solution$anova.table      <- sol_2$anova.tables.list[[mediator]]
    temp_solution$coef.tables      <- sol_2$coef.tables.list[[mediator]]
    temp_solution$pvalue.table     <- sol_2$pvalue.tables.list[[mediator]]
    temp_solution$conv             <- sol_2$conv.list[[mediator]]
    temp_solution$samples_l1_param <- sol_2$samples_l1.list[[mediator]]
    temp_solution$samples_l2_param <- sol_2$samples_l2.list[[mediator]]
    temp_solution$R2               <- sol_2$R2[[mediator]]
    temp_solution$tau_ySq          <- sol_2$tau_ySq[[mediator]]
    temp_solution$dMatrice         <- adapt.design.matrix(sol_2$dMatrice, mediator)
    temp_solution$mf1              <- adapt.mf1(sol_2$mf1, mediator)
    temp_solution$mf2              <- sol_2$mf2 

    sol <- BANOVA.mediation(sol_1, temp_solution, xvar=xvar, mediator=mediator,
                            individual = individual, return_effects = T)
    
    
   results[[mediator]] <- sol
  }
  
  final_result <- list()
  #######Report direct effects of the causal variable on the outcome#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0("Direct effects of the causal variable ", xvar, " on the outcome variable\n\n"))
  direct_effects_list <- results[[1]]$dir_effect
  prev_table <- 0
  for (i in 1:length(direct_effects_list)){
    new_table <- direct_effects_list[[i]]
    if (i == 2){
      cat(paste0("Simple effects of the causal variable ", xvar, "\n\n"))
    } 
    if (!identical(prev_table, new_table)){
      print(as.table(new_table), right=T)
      cat("\n")
    }
    prev_table <- new_table
  }
  final_result$dir_effect <- direct_effects_list
  
  #######Report (not print) direct effects of the mediator variables on the outcome#######
  mediator_names <- paste(mediators,  collapse=" and ")
  #cat(paste0(paste("Direct effects of mediators", mediator_names, "on the outcome variable\n")))
  for (mediator in mediators){
    final_result$m1_effects[[mediator]] <- results[[mediator]]$m1_effects
  }

  #######Report direct effects of the causal variable on mediator variables#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0(paste("Direct effects of the causal variable", xvar, "on the mediator variables\n\n")))
  for (mediator in mediators){
    cat((paste("Direct effects of", xvar, "on", mediator, "\n")))
    print(as.table(results[[mediator]]$m2_effects[[1]]), right=T)
    cat("\n")
    final_result$m2_effects[[mediator]] <- results[[mediator]]$m2_effects[[1]]
  }
  
  #######Report indirect effects of the causal variable#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0(paste("Indirect effects of the causal variable", xvar, "on the outcome variables\n\n")))
  for (mediator in mediators){
    cat((paste("Indirect effects of", xvar, "through", mediator, "\n")))
    indirect_effects_list <- results[[mediator]]$indir_effects
    prev_table <- 0
    for (i in 1:length(indirect_effects_list)){
      new_table <- indirect_effects_list[[i]]
      if (!identical(prev_table, new_table)){
        print(noquote(new_table), row.names = F, right=T)
        cat("\n")
      }
      prev_table <- new_table
    }
    
  
    final_result$indir_effects[[mediator]] <- results[[mediator]]$indir_effects
    
    # print(as.table(results[[mediator]]$indir_effects), right=T)
    # cat("\n")
  }
  
  #Report total indirect effects of the causal variable
  #Report total effects of the causal variable
  
  
  #return(results)
  #######Report total indirect effects of the causal variable#######
  cat(paste(strrep("-", 100), '\n'))
  
  #Extract indirect effects
  ind_eff_samples      <- list()
  ind_eff_sample_sizes <- c()
  for (mediator in mediators){
    ind_eff_samples <- results[[mediator]]$indirect_effects_samples
    if (length(ind_eff_samples) > 1){
      ###?
    } else {
      ind_eff_samples[[mediator]] <- as.data.frame(ind_eff_samples[[1]])
    }
    ind_eff_sample_sizes <- c(ind_eff_sample_sizes, dim(ind_eff_samples[[mediator]])[2])
  }
  common_samples <- min(ind_eff_sample_sizes)
  
  first_table <- ind_eff_samples[[1]][,1:common_samples]
  condtion    <- startsWith(colnames(first_table), "s_")
  matching_columns <- colnames(first_table)[!condtion]
  samples_columns  <- colnames(first_table)[condtion]
  n_matching_columns <- length(matching_columns)
  for (i in 1:(length(mediators)-1)){
    combined_samples <- merge(first_table[, 1:common_samples], 
                              ind_eff_samples[[i+1]][,1:common_samples], 
                              by = matching_columns)
    temp1 <- combined_samples[, paste0(samples_columns, ".x")]
    temp2 <- combined_samples[, paste0(samples_columns, ".y")]
    first_table[,(n_matching_columns+1):common_samples] <-  temp1 + temp2
  }
  colnames(combined_samples)
  
  rowMeans(first_table[(n_matching_columns+1):common_samples])
  result_table[ind,'mean'] <- round(rowMeans(), 4)
  result_table[ind,c('2.5%', '97.5%')] <- round(quantile(m_samples, probs = c(0.025, 0.975)),4)
  result_table[ind,'p.value'] <- ifelse(round(pValues(array(m_samples, dim = c(length(m_samples), 1))), 4) == 0, '<0.0001', round(pValues(array(m_samples, dim = c(length(m_samples), 1))), 4))
  
  
  
  cat(paste0(paste("Total indirect effects of the causal variable", xvar, "on the outcome variables\n\n")))
  

  #######Report total effects of the causal variable#######
  #cuculate total indirect effects of the variables
  calc.mediation.results <- function(sol, effects_of_interest){
    
  }
  return(final_result)
}
