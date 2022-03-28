isotopic_mass <- c("H"=1.007825,
                   "C"=12.0,
                   "N"=14.003074,
                   "O"=15.994915,
                   "F"=18.000938,
                   "P"=30.973762,
                   "S"=31.972071)
exact_mass <- function(formula, db=isotopic_mass){
  vector <- unlist(strsplit(formula, split=""))
  mass=0
  for(i in 1:length(vector)){
    if(vector[i] %in% names(db)){
      if(i == length(vector)){
        mass = mass + db[[vector[i]]]
        return(mass)
      }
      options (warn = -1)
      test <- as.numeric(vector[i+1])
      options (warn = 1)
      if(is.na(test) == T){
        mass = mass + db[[vector[i]]]
      }else{
        mass = mass + db[[vector[i]]] * test
      }
    }else{
      next
    }
  }
  return(mass)
}
