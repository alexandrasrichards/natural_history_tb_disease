model TB_model {
    param mort
    param min_out
    
    param min_sub
    param sub_min
    param sub_clin
    param clin_sub
    
    param phi
    
    state I 
    state t_med                         
    
    state Minmin                        (has_output = 0) 
    state Submin  
    state Clinmin  
    
    state Minsub  
    state Subsub                        
    state Clinsub  
    
    state Minclin  
    state Subclin  
    state Clinclin                      (has_output = 0) 
    
    state Minmin_sub                    (has_output = 0) 
    state Submin_sub 
    
    state Minsub_min  
    state Subsub_min                    (has_output = 0) 
    state Clinsub_min                   (has_output = 0) 
    
    state Minsub_clin                   (has_output = 0) 
    state Subsub_clin                   (has_output = 0) 
    state Clinsub_clin 
    
    state Subclin_sub 
    state Clinclin_sub                  (has_output = 0) 
    
    state Minmin_clin                   (has_output = 0) 
    state Submin_clin                   (has_output = 0) 
    state Clinmin_clin  
    
    state Minclin_min  
    state Subclin_min                   (has_output = 0) 
    state Clinclin_min                  (has_output = 0) 
    
    state MinminI                       (has_output = 0) 
    state SubminI                       (has_output = 0) 
    state ClinminI                      (has_output = 0) 
    state InfminI                           
    
    state MininfI  
    state SubinfI                       (has_output = 0) 
    state ClininfI                      (has_output = 0) 
    state InfinfI                       (has_output = 0) 
    
    state Minmin_infI                   (has_output = 0) 
    state Submin_infI                   (has_output = 0) 
    state Clinmin_infI                  (has_output = 0) 
    state Infmin_infI  
    
    state Mininf_minI  
    state Subinf_minI                   (has_output = 0) 
    state Clininf_minI                  (has_output = 0) 
    state Infinf_minI                   (has_output = 0) 
    
    state sub_over_clin_p
    state min_over_infect_p
    state duration
    
    dim series(32)
    
    obs sub_over_clin_obs[series]
    obs min_over_infect_obs[series]
    obs duration_obs[series]
    
    obs Subminobs[series]
    obs Minsubobs[series]
    obs Clinminobs[series]
    obs Minclinobs[series]
    
    obs Submin_subobs[series]
    obs Clinmin_clinobs[series]
    obs Minclin_minobs[series]
    
    obs InfminIobs[series]
    obs MininfIobs[series]
        
    obs Infmin_infIobs[series]
    obs Mininf_minIobs[series]
    
    input Submininput[series]
    input Minsubinput[series]
    input Clinmininput[series]
    input Minclininput[series]
    
    input Submin_subinput[series]
    input Clinmin_clininput[series]
    input Minclin_mininput[series]
    
    input InfminIinput[series]
    input MininfIinput[series]
        
    input Infmin_infIinput[series]
    input Mininf_minIinput[series]
    
    sub initial {
        Minmin <- 1 
        Submin <- 0
        Clinmin <- 0
        
        Minsub <- 0
        Subsub <- 1
        Clinsub <- 0
        
        Minclin <- 0
        Subclin <- 0
        Clinclin <- 1
        
        Minmin_sub <- 1
        Submin_sub <- 0
        
        Minsub_min <- 0
        Subsub_min <- 1
        Clinsub_min <- 0
        
        Minsub_clin <- 0
        Subsub_clin <- 1
        Clinsub_clin <- 0
        
        Subclin_sub <- 0
        Clinclin_sub <- 1
        
        Minmin_clin <- 1
        Submin_clin <- 0
        Clinmin_clin <- 0
        
        Minclin_min <- 0
        Subclin_min <- 0
        Clinclin_min <- 1
        
        MinminI <- 1
        SubminI <- 0
        ClinminI <- 0
        InfminI <- 0
        
        MininfI <- 0
        SubinfI <- phi
        ClininfI <- (1 - phi)
        InfinfI <- 1
        
        Minmin_infI <- 1
        Submin_infI <- 0
        Clinmin_infI <- 0
        Infmin_infI <- 0
        
        Mininf_minI <- 0
        Subinf_minI <- phi
        Clininf_minI <- (1 - phi)
        Infinf_minI <- 1
        
        I <- 1
    }
    
    sub parameter {
        min_out ~ uniform(0,12)
        min_sub ~ uniform(0,12)
        sub_min ~ uniform(0,12)
        sub_clin ~ uniform(0,12)
        clin_sub ~ uniform(0,12)
        phi ~ uniform(0,1)
        
        mort ~ gaussian(0.389, 0.028)
    }
    
    sub transition (delta = 0.1) {
        ode {
            dMinmin/dt = - min_sub * Minmin + sub_min * Submin - min_out * Minmin
            dSubmin/dt = min_sub * Minmin - sub_min * Submin - sub_clin * Submin + clin_sub * Clinmin
            dClinmin/dt = sub_clin * Submin - clin_sub * Clinmin - mort * Clinmin
            
            dMinsub/dt = - min_sub * Minsub + sub_min * Subsub - min_out * Minsub
            dSubsub/dt = min_sub * Minsub - sub_min * Subsub - sub_clin * Subsub + clin_sub * Clinsub
            dClinsub/dt = sub_clin * Subsub - clin_sub * Clinsub - mort * Clinsub
            
            dMinclin/dt = - min_sub * Minclin + sub_min * Subclin - min_out * Minclin
            dSubclin/dt = min_sub * Minclin - sub_min * Subclin - sub_clin * Subclin + clin_sub * Clinclin
            dClinclin/dt = sub_clin * Subclin - clin_sub * Clinclin - mort * Clinclin
            
            dMinmin_sub/dt = - min_sub * Minmin_sub - min_out * Minmin_sub
            dSubmin_sub/dt = min_sub * Minmin_sub 
            
            dMinsub_min/dt = sub_min * Subsub_min
            dSubsub_min/dt = - sub_min * Subsub_min - sub_clin * Subsub_min + clin_sub * Clinsub_min
            dClinsub_min/dt = sub_clin * Subsub_min - clin_sub * Clinsub_min - (mort * Clinsub_min)
            
            dMinsub_clin/dt = - min_sub * Minsub_clin + sub_min * Subsub_clin - min_out * Minsub_clin
            dSubsub_clin/dt = min_sub * Minsub_clin - sub_min * Subsub_clin - sub_clin * Subsub_clin
            dClinsub_clin/dt = sub_clin * Subsub_clin
            
            dSubclin_sub/dt = clin_sub * Clinclin_sub
            dClinclin_sub/dt = - clin_sub * Clinclin_sub - (mort * Clinclin_sub)
            
            dMinmin_clin/dt = - min_sub * Minmin_clin + sub_min * Submin_clin - min_out * Minmin_clin
            dSubmin_clin/dt = min_sub * Minmin_clin - sub_min * Submin_clin - sub_clin * Submin_clin
            dClinmin_clin/dt = sub_clin * Submin_clin 
            
            dMinclin_min/dt = sub_min * Subclin_min
            dSubclin_min/dt = - sub_min * Subclin_min - sub_clin * Subclin_min + clin_sub * Clinclin_min 
            dClinclin_min/dt = sub_clin * Subclin_min - clin_sub * Clinclin_min - (mort * Clinclin_min)
            
            
            dMinminI/dt = - min_sub * MinminI + sub_min * SubminI - min_out * MinminI
            dSubminI/dt = min_sub * MinminI - sub_min * SubminI - sub_clin * SubminI + clin_sub * ClinminI
            dClinminI/dt = sub_clin * SubminI - clin_sub * ClinminI - mort * ClinminI
            
            dMininfI/dt = - min_sub * MininfI + sub_min * SubinfI - min_out * MininfI
            dSubinfI/dt = min_sub * MininfI - sub_min * SubinfI - sub_clin * SubinfI + clin_sub * ClininfI
            dClininfI/dt = sub_clin * SubinfI - clin_sub * ClininfI - mort * ClininfI
            
            dMinmin_infI/dt = - min_sub * Minmin_infI - min_out * Minmin_infI 
            dSubmin_infI/dt = min_sub * Minmin_infI - sub_clin * Submin_infI + clin_sub * Clinmin_infI
            dClinmin_infI/dt = sub_clin * Submin_infI - clin_sub * Clinmin_infI
            
            dMininf_minI/dt = sub_min * Subinf_minI
            dSubinf_minI/dt = - sub_min * Subinf_minI - sub_clin * Subinf_minI + clin_sub * Clininf_minI
            dClininf_minI/dt = sub_clin * Subinf_minI - clin_sub * Clininf_minI - (mort * Clininf_minI)
        }
        
        InfminI <- SubminI + ClinminI
        InfinfI <- SubinfI + ClininfI
        Infmin_infI <- Submin_infI + Clinmin_infI
        Infinf_minI <- Subinf_minI + Clininf_minI
        
        sub_over_clin_p <- (clin_sub + mort) / sub_clin
        min_over_infect_p <- ((sub_min + sub_clin) * (clin_sub + mort) - (clin_sub * sub_clin))/(min_sub * (clin_sub + mort + sub_clin))
        duration <- -2 * log(2) / log(Subsub + Clinsub)
    }
    
    sub observation {
        Subminobs[series] ~ binomial(size = Submininput[series], prob = Submin)
        Minsubobs[series] ~ binomial(size = Minsubinput[series], prob = Minsub)
        Clinminobs[series] ~ binomial(size = Clinmininput[series], prob = Clinmin)
        Minclinobs[series] ~ binomial(size = Minclininput[series], prob = Minclin)
        
        Submin_subobs[series] ~ binomial(size = Submin_subinput[series], prob = Submin_sub)
        Clinmin_clinobs[series] ~ binomial(size = Clinmin_clininput[series], prob = Clinmin_clin)
        Minclin_minobs[series] ~ binomial(size = Minclin_mininput[series], prob = Minclin_min)
        
        InfminIobs[series] ~ binomial(size = InfminIinput[series], prob = InfminI)
        MininfIobs[series] ~ binomial(size = MininfIinput[series], prob = MininfI)
        
        Infmin_infIobs[series] ~ binomial(size = Infmin_infIinput[series], prob = Infmin_infI)
        Mininf_minIobs[series] ~ binomial(size = Mininf_minIinput[series], prob = Mininf_minI)

        sub_over_clin_obs[series] ~ gaussian(mean = sub_over_clin_p, std = 0.25)
        min_over_infect_obs[series] ~ gaussian(mean = min_over_infect_p, std = 0.5)
        duration_obs[series] ~ truncated_gaussian(mean = duration, std = 0.5)
    }
    
    sub proposal_parameter {
        min_out ~ truncated_gaussian(mean = min_out, std = 0.01, lower = 0, upper = 12)
        min_sub ~ truncated_gaussian(mean = min_sub, std = 0.01, lower = 0, upper = 12)
        sub_min ~ truncated_gaussian(mean = sub_min, std = 0.01, lower = 0, upper = 12)
        sub_clin ~ truncated_gaussian(mean = sub_clin, std = 0.01, lower = 0, upper = 12)
        clin_sub ~ truncated_gaussian(mean = clin_sub, std = 0.01, lower = 0, upper = 12)
        phi ~ truncated_gaussian(mean = phi, std = 0.01, lower = 0, upper = 1)
        
        mort ~ truncated_gaussian(mean = mort, std = 0.01, lower = 0)
    }
}
