## ---IBM Function-------------------------------------------------------------------- ##

initialise_population <- function(n_steps, n_people, prop_sub, prop_min = 0){
    
    if (prop_sub < 0 || prop_sub > 1) stop("prop_sub is out of the range 0 - 1")
    if (prop_sub + prop_min > 1) stop("prop people over 1")
    
    n_sub <- round(prop_sub * n_people)
    n_min <- round(prop_min * n_people)
    n_clin <- n_people - (n_sub + n_min)
    I <- matrix(data = NA, nrow = n_steps, ncol = n_people)
    I[1,] <- c(rep(x = "m", times = n_min), rep(x = "s", times = n_sub), rep(x = "c", times = n_clin))
    return(I)
}

sort_risks <- function(n_steps, n_people, params, method){
    if (isFALSE(dplyr::setequal(names(params), c("min_out","min_sub",  "sub_min",  "sub_clin", "clin_sub","mort")))) stop("incorrect parameter names - include all of min_sub, sub_min, sub_clin, clin_sub")
    
    if (isFALSE(method %in% c("median","fixed","random","one"))) stop("incorrect method - choose one of one, median, fixed, or random")
    
    params <- as.data.frame(params)
    
    if (method == "median") {
        R = median_risks(n_steps, n_people, params)
    } else if (method == "fixed") {
        R = fixed_risks(n_steps, n_people, params)
    } else if (method == "random") {
        R = random_risks(n_steps, n_people, params)
    } else if (method == "one") {
        R = one_risk(n_steps, n_people, params)
    }
    
    return(R)
}

median_risks <- function(n_steps, n_people, params){
    med_params <- apply(params, MARGIN = 2, FUN = median)
    med_risks <- param_to_risk(med_params)
    R <- list(recover = matrix(data = med_risks[["min_out"]], nrow = n_steps, ncol = n_people),
              min_sub = matrix(data = med_risks[["min_sub"]], nrow = n_steps, ncol = n_people),
              sub_min = matrix(data = med_risks[["sub_min"]], nrow = n_steps, ncol = n_people),
              sub_clin = matrix(data = med_risks[["sub_clin"]], nrow = n_steps, ncol = n_people),
              clin_sub = matrix(data = med_risks[["clin_sub"]], nrow = n_steps, ncol = n_people),
              death = matrix(data = med_risks[["mort"]], nrow = n_steps, ncol = n_people))
    
    return(R)
}

fixed_risks <- function(n_steps, n_people, params){
    risks <- param_to_risk(params)
    R <- list(recover = matrix(data = NA, nrow = n_steps, ncol = n_people),
              min_sub = matrix(data = NA, nrow = n_steps, ncol = n_people),
              sub_min = matrix(data = NA, nrow = n_steps, ncol = n_people),
              sub_clin = matrix(data = NA, nrow = n_steps, ncol = n_people),
              clin_sub = matrix(data = NA, nrow = n_steps, ncol = n_people),
              death = matrix(data = NA, nrow = n_steps, ncol = n_people))
    
    for (i in 1:n_people) {
        R$recover[,i] = rep(base::sample(x = risks$min_out, size  = 1),n_steps)
        R$min_sub[,i] = rep(base::sample(x = risks$min_sub, size  = 1),n_steps)
        R$sub_min[,i] = rep(base::sample(x = risks$sub_min, size  = 1),n_steps)
        R$sub_clin[,i] = rep(base::sample(x = risks$sub_clin, size  = 1), n_steps)
        R$clin_sub[,i] = rep(base::sample(x = risks$clin_sub, size  = 1), n_steps)
        R$death[,i] = rep(base::sample(x = risks$mort, size  = 1),n_steps)
    }
    
    return(R)
}

random_risks <- function(n_steps, n_people, params){
    risks <- param_to_risk(params)
    R <- list(recover = matrix(data = NA, nrow = n_steps, ncol = n_people),
              min_sub = matrix(data = NA, nrow = n_steps, ncol = n_people),
              sub_min = matrix(data = NA, nrow = n_steps, ncol = n_people),
              sub_clin = matrix(data = NA, nrow = n_steps, ncol = n_people),
              clin_sub = matrix(data = NA, nrow = n_steps, ncol = n_people),
              death = matrix(data = NA, nrow = n_steps, ncol = n_people))
    for (i in 1:n_people) {
        R$recover[,i] = base::sample(x = risks$min_out, size  = n_steps, replace = TRUE)
        R$min_sub[,i] = base::sample(x = risks$min_sub, size  = n_steps, replace = TRUE)
        R$sub_min[,i] = base::sample(x = risks$sub_min, size  = n_steps, replace = TRUE)
        R$sub_clin[,i] = base::sample(x = risks$sub_clin, size  = n_steps, replace = TRUE)
        R$clin_sub[,i] = base::sample(x = risks$clin_sub, size  = n_steps, replace = TRUE)
        R$death[,i] = base::sample(x = risks$mort, size  = n_steps, replace = TRUE)
    }
    
    return(R)
}

one_risk <- function(n_steps, n_people, params) {
    risks <- param_to_risk(params)
    R <- list(recover = matrix(data = rep(base::sample(x = risks$min_out, size  = 1),n_steps), nrow = n_steps, ncol = n_people),
              min_sub = matrix(data = rep(base::sample(x = risks$min_sub, size  = 1),n_steps), nrow = n_steps, ncol = n_people),
              sub_min = matrix(data = rep(base::sample(x = risks$sub_min, size  = 1),n_steps), nrow = n_steps, ncol = n_people),
              sub_clin = matrix(data = rep(base::sample(x = risks$sub_clin, size  = 1),n_steps), nrow = n_steps, ncol = n_people),
              clin_sub = matrix(data = rep(base::sample(x = risks$clin_sub, size  = 1),n_steps), nrow = n_steps, ncol = n_people),
              death = matrix(data = rep(base::sample(x = risks$mort, size  = 1),n_steps), nrow = n_steps, ncol = n_people))
    
    return(R)
}

param_to_risk <- function(param){
    param_new <- 1 - exp(-param/12)
    return(param_new)
}

update_population <- function(X, t, n_people, risks, treatment){
    
    for (i in 1:n_people) {
       if (X[t - 1,i] == "m") {
           X[t,i] <-  update_minimal(i,t - 1,risks)
       } else if (X[t - 1,i] == "s") {
           X[t,i] <-  update_subclinical(i,t - 1,risks)
       } else if (X[t - 1,i] == "c") {
           X[t,i] <-  update_clinical(i,t - 1,risks, treatment)
       } else {
           X[t,i] <-  "-"
       }
    }
    return(X)
}

update_minimal <- function(person, time, risks){
    min_out <-  risks$recover[time, person]
    min_sub <-  risks$min_sub[time, person]
    prob <- runif(1,0,1)
    
    if (prob < min_out) {
        return("r")
    } else if (prob < min_sub + min_out) {
        return("s")
    } else {
        return("m")
    }
}

update_subclinical <- function(person, time, risks){
    sub_min <- risks$sub_min[time, person]
    sub_clin <- risks$sub_clin[time, person]
    
    prob <- runif(1,0,1)
    
    if (prob < sub_min) {
        return("m")
    } else if (prob < sub_min + sub_clin) {
        return("c")
    } else {
        return("s")
    }
}

update_clinical <- function(person, time, risks, treatment){
    treat <- param_to_risk(treatment)
    death <- risks$death[time, person]
    clin_sub <- risks$clin_sub[time, person]
    
    prob <- runif(1,0,1)
    
    if (prob < death) {
        return("d")
    } else if (prob < death + clin_sub) {
        return("s")
    } else {
        prob_treat <- runif(1,0,1)
        if (prob_treat < treat) {
            return("t")
        } else {
           return("c") 
        }
        
    }
}

full_simulation <- function(n_steps, n_people, prop_sub, params, method, treatment = 0, prop_min = 0){
    pop <- initialise_population(n_steps, n_people, prop_sub, prop_min)
    risks <- sort_risks(n_steps, n_people, params, method)
    
    for (i in seq(2,n_steps)) {
        pop <- update_population(pop,i,n_people,risks, treatment)
    }
    return(pop)
}


 



                           