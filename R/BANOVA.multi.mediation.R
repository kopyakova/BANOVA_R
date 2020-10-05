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
  calculate.total.indirect.effects <- function(used_tables_index){
    ind_eff_samples      <- list()
    
    #Exract relevant samples of indirect effects
    for (mediator in mediators){
      used_tables <- used_tables_index[[mediator]]
      ind_eff_samples[[mediator]] <- results[[mediator]]$indirect_effects_samples[used_tables]
    }
    
    n_tables <- length(used_tables_index[[1]])
    total_indir_eff_samples_list <- list()
    total_indir_eff_table_list <- list()
    for (i in 1:n_tables){
      ind_eff_n_samples <- c()
      for (mediator in mediators){
        n_samples <- dim(ind_eff_samples[[mediator]][[i]])[2]
        ind_eff_n_samples <- c(ind_eff_n_samples, n_samples)
      }
      common_samples <- min(ind_eff_n_samples)
      
      #Combine the tables
      combined_table <- ind_eff_samples[[1]][[i]][, 1:common_samples]
      smpl_indicator <- startsWith(colnames(combined_table), "s_")
      smpl_columns   <- colnames(combined_table)[smpl_indicator]
      num_non_smpl_columns <- sum(!smpl_indicator) + 1
     
      common_columns    <- colnames(combined_table)[!smpl_indicator]
      common_columns    <- common_columns[-1] #drop the intercept
      
      for (n in 2:n_mediators){
        common_columns <- intersect(common_columns, colnames(ind_eff_samples[[n]][[i]])[!smpl_indicator])
      }
      num_common_columns <- length(common_columns)
    
      for (n in 2:n_mediators){
        temp_table1 <- combined_table[, 1:common_samples]
        temp_table2 <- ind_eff_samples[[n]][[i]][, 1:common_samples]
        combined_samples <- merge(temp_table1[, c(common_columns, smpl_columns)], 
                                  temp_table2[, c(common_columns, smpl_columns)], 
                                  by = common_columns)
        temp1 <- combined_samples[, paste0(smpl_columns, ".x")]
        temp2 <- combined_samples[, paste0(smpl_columns, ".y")]
        combined_table[, num_non_smpl_columns:common_samples] <-  temp1 + temp2
      }
      total_indirect_effects_samples <- combined_table[, c(smpl_columns)]
      
      total_indirect_effects <- data.frame(matrix(NA, nrow = nrow(combined_table), ncol = num_common_columns+4))
      colnames(total_indirect_effects) <- c(common_columns, "mean","2.5%","97.5%","p.value")
      
      total_indirect_effects[, 1:num_common_columns] <- combined_table[, common_columns]
      total_indirect_effects[, "mean"] <- round(apply(total_indirect_effects_samples, 1, mean), 4)
      quantiles    <- round(apply(total_indirect_effects_samples, 1, quantile, 
                                  probs = c(0.025, 0.975), type = 3, na.rm = FALSE), 4)
      total_indirect_effects[, "2.5%"]  <- quantiles["2.5%",]
      total_indirect_effects[, "97.5%"] <- quantiles["97.5%",]
      p_values <- round(pValues(t(total_indirect_effects_samples)), 4)
      total_indirect_effects[, "p.value"] <-  apply(p_values, 1,
                                                    FUN = function(x) ifelse(x == 0, '<0.0001', x))
      
      temp_effects <- cbind(combined_table[, common_columns], total_indirect_effects_samples)
      colnames(temp_effects) <- c(common_columns, colnames(total_indirect_effects_samples))
      total_indir_eff_samples_list[[i]] <- temp_effects
      total_indir_eff_table_list[[i]]   <- total_indirect_effects
    }
    return(list(total_indirect_effects = total_indir_eff_table_list,
           total_indirect_effects_samples = total_indir_eff_samples_list))
  }
  
  combine_direct_and_indirect_effects <- function(dir_eff_samples_df, total_indir_eff_samples){
    n_indir_eff_tables <- length(total_indir_eff_samples)
    total_eff_table_list <- list()
    for (i in 1:n_indir_eff_tables){
      dir_eff_samples <- dir_eff_samples_df
      indir_eff_samples_df <- data.frame(total_indir_eff_samples[[i]])
      
      smpl_indicator <- startsWith(colnames(indir_eff_samples_df), "s_")
      smpl_col_names <- colnames(indir_eff_samples_df)[smpl_indicator]
      
      common_samples     <- min(ncol(dir_eff_samples), ncol(indir_eff_samples_df))
      
      colnames_dir_eff   <- colnames(dir_eff_samples)
      colnames_indir_eff <- colnames(indir_eff_samples_df)
      common_columns     <- intersect(colnames_dir_eff, colnames_indir_eff)
      n_common_columns   <- length(common_columns)
      smpl_index         <- (n_common_columns+1):common_samples
      
      
      colnames(dir_eff_samples)[1:common_samples] <- colnames_indir_eff[1:common_samples]
      
      total_effect_samples <- merge(dir_eff_samples[, 1:common_samples], 
                                    indir_eff_samples_df[, 1:common_samples], 
                                    by = common_columns)
      
      temp1 <- total_effect_samples[, paste0(smpl_col_names, ".x")]
      temp2 <- total_effect_samples[, paste0(smpl_col_names, ".y")]
      
      total_effect_samples[, smpl_index] <-  temp1 + temp2
      
      #total_effects_samples <- total_effect_samples[, smpl_index]
      
      total_effects <- data.frame(matrix(NA, nrow = nrow(total_effect_samples), ncol = n_common_columns+4))
      colnames(total_effects) <- c(common_columns, "mean","2.5%","97.5%","p.value")
      
      total_effects[, 1:n_common_columns] <- total_effect_samples[, common_columns]
      total_effects[, "mean"] <- round(apply(total_effect_samples[, smpl_index], 1, mean), 4)
      quantiles    <- round(apply(total_effect_samples[, smpl_index], 1, quantile, 
                                  probs = c(0.025, 0.975), type = 3, na.rm = FALSE), 4)
      total_effects[, "2.5%"]  <- quantiles["2.5%",]
      total_effects[, "97.5%"] <- quantiles["97.5%",]
      p_values <- round(pValues(t(total_effect_samples[, smpl_index])), 4)
      total_effects[, "p.value"] <-  apply(p_values, 1,
                                           FUN = function(x) ifelse(x == 0, '<0.0001', x))
      
      total_eff_table_list[[i]]   <- total_effects
    }
    return(total_eff_table_list)
  }
  
  print.result <- function(list_with_results, extra_title = NULL, final_result, list_name,
                           extra_list_name = NULL, skip_n_last_cols = 0, return_table_index = F){
    check.if.tables.are.identical <- function(table1, table2){
      table1_rwnames <- rownames(table1)
      table2_rwnames <- rownames(table2)
      rownames(table1) <- NULL
      rownames(table2) <- NULL
      return(identical(table1, table2))
    }
    print.table <- function(new_table, skip_n_last_cols, extra_title, prev_table = 0){
      check.table.columns <- function(new_table, skip_n_last_cols){
        remove_table <- F
        n_columns <- ncol(new_table)
        if (is.null(ncol(new_table))){
          #if the table is a vector
          n_cols_to_check <- length(new_table) - skip_n_last_cols
          if (n_cols_to_check != 0){
            for (j in 1:n_cols_to_check){
              if (all(new_table[j] == '0') || all(new_table[j] == 0))
                remove_table <- T
            }
          }
        } else {
          #otherwise 
          n_cols_to_check <- n_columns - skip_n_last_cols
          for (j in 1:n_cols_to_check){
            if (all(new_table[, j] == '0') || all(new_table[, j] == 0))
              remove_table <- T
          }
        }
        return(remove_table)
      }
      
      remove_table <- check.table.columns(new_table, skip_n_last_cols)
      if (!remove_table){
        if (!check.if.tables.are.identical(prev_table, new_table)){
          if(!is.null(extra_title)){
            cat(extra_title)
          }
          print(noquote(new_table), row.names = F, right=T)
          cat("\n")
        }
        
      }
      return(list(remove_table = remove_table, prev_table = new_table))
    }
    
    
    list_length <- length(list_with_results)
    temp_list   <- list()
    prev_table  <- 0
    if (list_length > 1){
      counter    <- 1
      used_tables_index <- c()
      for (i in 1:list_length){
        new_table <- list_with_results[[i]]
        temp_res  <- print.table(new_table, skip_n_last_cols, extra_title, prev_table)
        remove_table <- temp_res$remove_table
        if (!remove_table){
          if (!check.if.tables.are.identical(prev_table, new_table)){
            temp_list[[counter]] <- new_table
            used_tables_index <- c(used_tables_index, i)
            counter <- counter + 1
            prev_table   <- temp_res$prev_table
          }
        }
      }
    } else {
      new_table    <- list_with_results[[1]]
      temp_res     <- print.table(new_table, skip_n_last_cols, extra_title, prev_table)
      remove_table <- temp_res$remove_table
      if (!remove_table){
        temp_list[[1]] <- new_table
        used_tables_index <- 1
      }
    } 
    
    temp_list_length <- length(temp_list)
    if (temp_list_length > 1){
      counter <- 1
      for (i in 1:temp_list_length){
        if(!is.null(extra_list_name)){
          final_result[[list_name]][[extra_list_name]][[counter]] <- temp_list[[i]]
        } else {
          final_result[[list_name]][[counter]] <- temp_list[[i]]
        }
        counter <- counter + 1
      }
    } else {
      if(!is.null(extra_list_name)){
        final_result[[list_name]][[extra_list_name]] <- temp_list[[1]]
      } else {
        final_result[[list_name]] <- temp_list[[1]]
      }
    }
    if (return_table_index){
      return(list(final_result = final_result, used_tables_index = used_tables_index))
    } else {
      return(final_result)
    }
  }

  #Check the class of the model with multiple dependent variables
  if(sol_2$model_name != "BANOVA.multiNormal") 
    stop('The mediator must follow the multivariate Normal distribution, use BANOVA multiNormal models instead.')
  
  n_mediators <- length(mediators)
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

    sol <- BANOVA.mediation(sol_1, sol_2 = temp_solution, xvar=xvar, mediator=mediator,
                            individual = individual, return_effects = T)
    
   results[[mediator]] <- sol
  }

  final_result <- list()
  #######Report direct effects of the causal variable on the outcome#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0("Direct effects of the causal variable ", xvar, " on the outcome variable\n\n"))
  final_result <- print.result(list_with_results = results[[1]]$dir_effect, 
                               final_result = final_result, 
                               list_name = "dir_effect", skip_n_last_cols = 3)
  
  #######Report direct effects of the mediator variables on the outcome#######
  cat(paste(strrep("-", 100), '\n'))
  mediator_names <- paste(mediators,  collapse=" and ")
  cat(paste0(paste("Direct effects of mediators", mediator_names, "on the outcome variable\n")))
  for (mediator in mediators){
    final_result <- print.result(results[[mediator]]$m1_effects, final_result = final_result, 
                                 extra_title = paste0("Direct effects of ", mediator, "\n"),
                                 list_name = "m1_effects", extra_list_name = mediator,
                                 skip_n_last_cols = 3)
  }

  #######Report direct effects of the causal variable on mediator variables#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0(paste("Direct effects of the causal variable", xvar, "on the mediator variables\n\n")))
  for (mediator in mediators){
    final_result <- print.result(list_with_results = results[[mediator]]$m2_effects, final_result = final_result, 
                                 extra_title = paste("Direct effects of", xvar, "on", mediator, "\n"),
                                 list_name = "m2_effects", extra_list_name = mediator,
                                 skip_n_last_cols = 3)
  }
  
  #######Report indirect effects of the causal variable#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0(paste("Indirect effects of the causal variable", xvar, "on the outcome variables\n\n")))
  used_tables_index <- list()
  for (mediator in mediators){
    temp_result <- print.result(list_with_results = results[[mediator]]$indir_effects, 
                                 final_result = final_result, 
                                 extra_title = paste("Indirect effects of", xvar, "via", mediator, "\n"),
                                 list_name = "indir_effects", extra_list_name = mediator,
                                 skip_n_last_cols = 3, return_table_index = T)
    final_result <- temp_result$final_result
    used_tables_index[[mediator]] <- temp_result$used_tables_index
  }

  #######Report total indirect effects of the causal variable#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0(paste("Total indirect effects of the causal variable", xvar, 
                   "on the outcome variables\n\n")))
  total_indirect_effects_results <- calculate.total.indirect.effects(used_tables_index)
  n_tables_with_indirect_effects <- length(used_tables_index[[1]])
  if (n_tables_with_indirect_effects > 1){
    for (i in 1:n_tables_with_indirect_effects){
      table <- total_indirect_effects_results$total_indirect_effects[[i]]
      print(noquote(table), row.names = F, right=T)
      cat('\n')
    }
    final_result$total_indir_effect <- total_indirect_effects_results$total_indirect_effects
  } else{
      table <- total_indirect_effects_results$total_indirect_effects[[1]]
      print(noquote(table), row.names = F, right=T)
      cat('\n')
      
      final_result$total_indir_effects <- table
  }

 #######Report total effects of the causal variable#######
  cat(paste(strrep("-", 100), '\n'))
  cat(paste0(paste("Total effects of the causal variable", xvar, 
                   "on the outcome variables\n\n")))
  total_indir_eff_samples <- total_indirect_effects_results$total_indirect_effects_samples
  
  dir_eff_list <- results[[1]]$direct_effects_samples
  dir_eff_list <- dir_eff_list[[length(dir_eff_list)]] #select the last element

  dir_eff_samples     <- dir_eff_list$samples
  dim_dir_eff_samples <- dim(dir_eff_samples)
  dir_eff_samples_df  <- data.frame(matrix(dir_eff_samples, 
                                          dim_dir_eff_samples[1]*dim_dir_eff_samples[2],
                                          dim_dir_eff_samples[3]))
  dir_effect_names <- dir_eff_list$index_name
  temp_colnames <- colnames(dir_effect_names)
  if ("(Intercept)" %in% colnames(dir_effect_names)){
    dir_effect_names <- dir_effect_names[,!colnames(dir_effect_names) %in% "(Intercept)"]
    if (is.null(dim(dir_effect_names))){
      dir_effect_names <- as.data.frame(t(as.matrix(dir_effect_names)))
      colnames(dir_effect_names) <- temp_colnames[!temp_colnames %in% "(Intercept)"]
    }
  }
 
 dir_eff_samples_df <- cbind(dir_effect_names, dir_eff_samples_df)
 total_effects <- combine_direct_and_indirect_effects(dir_eff_samples_df, total_indir_eff_samples)
 n_tables_with_total_effects <- length(total_effects)
 if (n_tables_with_total_effects > 1){
   for (i in 1:n_tables_with_total_effects){
     table <- total_effects[[i]]
     print(noquote(table), row.names = F, right=T)
     cat('\n')
   }
   final_result$total_effects <- total_effects
 } else{
   table <- total_effects[[1]]
   print(noquote(table), row.names = F, right=T)
   cat('\n')
   
   final_result$total_effects <- table
 }


  return(final_result)
}
