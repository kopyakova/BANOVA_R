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

    cat(paste0("\nResults for the ", mediator, " mediator\n"))
    sol <- BANOVA.mediation(sol_1, temp_solution, xvar=xvar, mediator=mediator,
                            individual = individual)
    print(sol)
    results[[mediator]] <- sol
  }
  return(results)
}
