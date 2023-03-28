
fit_plot <- function(transition, posterior, data, title){
    
    data <-  data[which(data$Transition == transition),]
    
    posterior <-  posterior[which(posterior$var == transition),]
    colour <- ifelse(transition %in% (c("Submin","Minsub","Clinmin","Minclin","Clinsub","Subclin","InfminI","MininfI")),
                     "dodgerblue4",
                     "tomato4")
                     
    plot <- ggplot() +
        geom_ribbon(data = posterior, aes(x = time, ymin = `2.5%`, ymax = `97.5%`), fill = colour, alpha = 0.5) +
        geom_line(data = posterior, aes(x = time, y = Median), colour = colour) +
        xlim(c(0,15)) +
        ylim(c(0,1)) +
        ggtitle(title) +
        xlab("") +
        ylab("") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
    
    if (nrow(data) > 0) {
        
        plot <- plot +
            geom_pointrange(data = data, aes(x = times, y = value,ymin = error_low, ymax = error_upp), colour = colour, shape = 20, size = 0.25)#  + 
        #theme(legend.position="bottom") +
        #geom_text(data = data, aes(label = record.id,x = times, y = value),hjust = 0, vjust = 0, size = 2.5)
        
        
    } 
    
    return(plot)
}


total_fit_plots <- function(data, posterior){
    fit_layout <- png::readPNG(("Figures/fit-model-structure.png"))
    min_sub <- png::readPNG(("Figures/fit-min-sub.png"))
    sub_clin <- png::readPNG(("Figures/fit-sub-clin.png"))
    min_clin <- png::readPNG(("Figures/fit-min-clin.png"))
    min_inf <- png::readPNG(("Figures/fit-min-inf.png" ))
    
    fit_layout <- ggplot() + background_image(fit_layout) + theme_minimal()
    min_sub <- ggplot() + background_image(min_sub) + theme_minimal()
    sub_clin <- ggplot() + background_image(sub_clin) + theme_minimal()
    min_clin <- ggplot() + background_image(min_clin) + theme_minimal()
    min_inf <- ggplot() + background_image(min_inf) + theme_minimal()
    
    data$value = data$n.of.y_weight/data$d.of.y_weight
    data$error_low = as.double(unlist(binom.confint(data$n.of.y_weight,data$d.of.y_weight,conf.level = 0.95, method = "wilson")["lower"]))
    data$error_upp = as.double(unlist(binom.confint(data$n.of.y_weight,data$d.of.y_weight,conf.level = 0.95, method = "wilson")["upper"]))
    
    
    min_sub_cs <-  fit_plot("Submin", posterior, data, "Minimal to Subclinical")    
    sub_min_cs <-  fit_plot("Minsub", posterior, data, "Subclinical to Minimal")  
    sub_clin_cs <-  fit_plot("Clinsub", posterior, data, "Subclinical to Clinical")  
    clin_sub_cs <-  fit_plot("Subclin", posterior, data, "Clinical to Subclinical")  
    min_clin_cs <-  fit_plot("Clinmin", posterior, data, "Minimal to Clinical")  
    clin_min_cs <-  fit_plot("Minclin", posterior, data, "Clinical to Minimal") 
    min_inf_cs <- fit_plot("InfminI", posterior, data, "Minimal to Infectious")
    inf_min_cs <- fit_plot("MininfI",posterior, data, "Infectious to Minimal")
    
    min_sub_tte <-  fit_plot("Submin_sub", posterior, data, "Minimal to Subclinical")    
    sub_min_tte <-  fit_plot("Minsub_min", posterior, data, "Subclinical to Minimal")  
    sub_clin_tte <-  fit_plot("Clinsub_clin", posterior, data, "Subclinical to Clinical")  
    clin_sub_tte <-  fit_plot("Subclin_sub", posterior, data, "Clinical to Subclinical")  
    min_clin_tte <-  fit_plot("Clinmin_clin", posterior, data, "Minimal to Clinical")  
    clin_min_tte <-  fit_plot("Minclin_min", posterior, data, "Clinical to Minimal") 
    min_inf_tte <- fit_plot("Infmin_infI", posterior, data, "Minimal to Infectious")
    inf_min_tte <- fit_plot("Mininf_minI",posterior, data, "Infectious to Minimal")
    
    
    layout <- "AAAAAAAAAAAAAAAAAAA
            FFFFGGGGBBBHHHHIIII
            JJJJKKKKCCCLLLLMMMM
            NNNNOOOODDDPPPPQQQQ
            RRRRSSSSEEETTTTUUUU"
    
    fit_plot <- (fit_layout + 
                     min_sub +
                     sub_clin +
                     min_clin +
                     min_inf +
                     min_sub_tte + sub_min_tte + 
                     min_sub_cs + sub_min_cs +
                     sub_clin_tte + clin_sub_tte +
                     sub_clin_cs + clin_sub_cs +
                     min_clin_tte + clin_min_tte +
                     min_clin_cs + clin_min_cs +
                     min_inf_tte + inf_min_tte +
                     min_inf_cs + inf_min_cs +
                     plot_layout(design = layout)) +
        theme_minimal()
    
    x_lab <- "Time (years)"
    y_lab <- "Proportion moved from initial to final state"
    
    x_lab <- 
        ggplot() + 
        annotate(geom = "text", x = 1, y = 1, label = x_lab) +
        coord_cartesian(clip = "off") +
        theme_void()
    
    y_lab <- 
        ggplot() + 
        annotate(geom = "text", x = 1, y = 1, label = y_lab, angle = 90) +
        coord_cartesian(clip = "off") +
        theme_void()
    
    fit_plot <- (fit_plot / x_lab + plot_layout(heights = c(1, 0.01)))
    fit_plot <- (y_lab + fit_plot + plot_layout(widths = c(1, 100)))
    
    return(fit_plot)
}

