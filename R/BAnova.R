BAnova <-
function(x){
  if (is.list(x$anova.table)){
    for (i in 1:length(x$anova.table)){
      cat('\nChoice: ', i, '\n')
      print(x$anova.table[[i]])
    }
  }else{
    print(x$anova.table)
  }
}
