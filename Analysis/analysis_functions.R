### ---IBM analysis------------------------------------------------------------------ ###
summarise_simulation <- function(n_steps, sim) {
    S <- matrix(data = NA, nrow = n_steps, ncol = 6)
    colnames(S) <- c("r","m","s","c","t","d")
    
    S[,1] <-  apply(sim, 1, function(x) table(x)["r"])
    S[,2] <-  apply(sim, 1, function(x) table(x)["m"])
    S[,3] <-  apply(sim, 1, function(x) table(x)["s"])
    S[,4] <-  apply(sim, 1, function(x) table(x)["c"])
    S[,5] <-  apply(sim, 1, function(x) table(x)["t"])
    S[,6] <-  apply(sim, 1, function(x) table(x)["d"])
    
    S[which(is.na(S))] <- 0
    
    for (i in seq(2,n_steps)) { # create cumulative count for r, t, and d
        S[i,1] <- S[i - 1,1] + S[i,1]
        S[i,5] <- S[i - 1,5] + S[i,5]
        S[i,6] <- S[i - 1,6] + S[i,6]
    }
    
    return(S)
}

find_median_disease <- function(sim_sum, n_people) {
    M <- sim_sum[,1] + sim_sum[,5] + sim_sum[,6]
    med <- min(which(M >= n_people/2)) - 1
    
    return(med)
}

find_median_infectious <- function(sim_sum, n_people) {
    M <- sim_sum[,1] + sim_sum[,2] + sim_sum[,5] + sim_sum[,6]
    med <- min(which(M >= n_people/2)) - 1
    
    return(med)
}

bar_chart_disease <- function(sim_sum, n_people) {
    med <- find_median_disease(sim_sum, n_people)
    sim_sum <- sim_sum[,c("m","s","c")]
    sim_melt = as.data.frame(melt(sim_sum))
    names(sim_melt) <- c("times","state","value")
    sim_melt$times <- sim_melt$times - 1
    pal <- c("c" = '#F8766D',
             "s" = '#00BFC4',
             "m" = '#C77CFF')
    
    p <- ggplot(sim_melt, aes(x = times, fill = factor(state, levels = c("c","s","m")))) +
        geom_bar(aes(y = value), stat = "identity") +
        geom_hline(yintercept = n_people/2,linetype = "dashed") +
        geom_linerange(x = med,ymin = 0, ymax = n_people/2, linetype = "dashed") +
        theme_bw() +
        xlab("Time (months)") +
        scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) * n_people,
                           labels = c(0,10,20,30,40,50,60,70,80,90,100),
                           name = "Percentage of cohort") +
        scale_fill_manual(values = pal, labels = c("Clinical","Subclinical","Minimal")) +
        scale_colour_manual(values = pal) +
        theme(legend.position = "right", legend.title = element_blank()) +
        annotate(
            geom = "text",
            label = paste("Median time:\n", med, " months", sep = ""),
            x = 75,
            y = 0.75*n_people
        )
    
    
    return(p)
    
}

bar_chart_infectious <- function(sim_sum, n_people) {
    med <- find_median_infectious(sim_sum, n_people)
    sim_sum <- sim_sum[,c("s","c")]
    sim_melt <-  as.data.frame(melt(sim_sum))
    names(sim_melt) <- c("times","state","value")
    sim_melt$state = factor(sim_melt$state, levels = c("c","s"))
    sim_melt$times <- sim_melt$times - 1
    pal <- c("c" = '#F8766D',
             "s" = '#00BFC4',
             "m" = '#C77CFF')
    
    p <- ggplot(sim_melt, aes(x = times, fill = state)) +
        geom_bar(aes(y = value), stat = "identity") +
        geom_hline(yintercept = n_people/2,linetype = "dashed") +
        geom_linerange(x = med,ymin = 0, ymax = n_people/2, linetype = "dashed") +
        theme_bw() +
        xlab("Time (months)") +
        scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) * n_people,
                           labels = c(0,10,20,30,40,50,60,70,80,90,100),
                           name = "Percentage of cohort") +
        scale_fill_manual(values = pal, labels = c("Clinical","Subclinical","Minimal")) +
        scale_colour_manual(values = pal) +
        theme(legend.position = "right", legend.title = element_blank()) +
        annotate(
            geom = "text",
            label = paste("Median time\n", med, " months", sep = ""),
            x = 75,
            y = 0.75 * n_people
        ) 
    return(p)
}

time_clinical_before_event <- function(sim, event = "d") {
    if (event %in% c("d","death")) {
        e = "d"
    } else if (event %in% c("t","treat","treatment")) {
        e = "t"
    } else {
        e = "d"
        warning("invalid event given, calculating time before death")
    }
    
    events <- which(apply(sim, 2, function(x) e %in% x))
    I <- sim[,events]
    
    n_events <- length(events)
    symp <- rep(0,n_events)
    for (i in seq(1,n_events)) {
        count = which(I[,i] == e)
        while (I[count - 1,i] == "c") {
            symp[i] = symp[i] + 1
            count = count - 1
            if (count == 1) {
                break()
            }
        }
        
    }
    return(symp)
}

time_clinical_before_recovery <- function(sim) {
    n_people <- ncol(sim)
    n_months <- nrow(sim)
    symp <- c()
    
    for (i in seq(1,n_people)) {
        for (j in seq(2,n_months)) {
            if (sim[j,i] == "s" && sim[j - 1,i] == "c") {
                track = j
                count = 0
                while (sim[track - 1,i] == "c" && track > 1) {
                    count = count + 1
                    track = track - 1
                }
                symp <- c(symp, count)
            }
        }
    }
    return(symp)
    
}

### ---Trajectories------------------------------------------------------------------ ###

report_trajectories <- function(simulation, n_people, years, threshold_months = 9, threshold_changes = 2) {
    
    if (n_people > ncol(simulation)) {
        n_people = ncol(simulation)
        warning("fewer people in simulation than requested, analysing on full simulation size")
    }
    
    if (years > (12 * nrow(simulation))) {
        years = floor(nrow(simulation))
        warning("fewer years in simulation than requested, analysing on full years in simulation")
    }
    
    if (threshold_months < 0 || threshold_months > 12) stop("provide threshold_months in range 0 - 12, or use default of 9")
    
    if (threshold_changes < 0 || threshold_changes > 12) stop("provide threshold_changes in range 0 - 12, or use default of 2")
    
    Tr <- matrix(data = NA ,nrow = years + 1, ncol = n_people)
    
    
    for (p in seq(1,n_people)) {
        
        if (simulation[1,p] == "m") {
            Tr[1,p] = "minimal"
        } else if (simulation[1,p] == "s") {
            Tr[1,p] = "subclinical"
        } else if (simulation[1,p] == "c") {
            Tr[1,p] = "clinical"
        }
        
        for (y in seq(1,years)) {
            start <- (y * 12) - 10
            end <- (y * 12) + 1
                if (Tr[y,p] %in% c("treat","death","recover")) {
                    Tr[y + 1,p] = Tr[y,p]
                } else {
                    Tr[y + 1,p] <- undulation(simulation[(start:end),p], threshold_months, threshold_changes)
                }
            
        }
    }
    return(Tr)
    
}

undulation <- function(data, threshold_months, threshold_changes) {
    treat <- length(which(data == "t"))
    death <- length(which(data == "d"))
    recover <- length(which(data == "r"))
    clin <- length(which(data == "c"))
    sub <- length(which(data == "s"))
    min <- length(which(data == "m"))
    changes <- count_changes(data)
    
    if (treat == 1) {
        return("treat")
    } else if (death == 1) {
        return("death")
    } else if (recover == 1) {
        return("recover")
    } else if (clin >= threshold_months && changes <= threshold_changes) {
        return("clinical")
    } else if (sub >= threshold_months && changes <= threshold_changes) {
        return("subclinical")
    } else if (min >= threshold_months && changes <= threshold_changes) {
        return("minimal")
    } else {
        return("undulation")
    }
}

count_changes <- function(data){
    changes <- 0
    old <- data[1]
    for (i in seq(2,12)) {
        new <- data[i]
        if (new != old) {
            changes <-  changes + 1
        }
        old <- new
    }
    return(changes)
}

summarise_trajectories <- function(trajectories, years) {
    Tr <- matrix(data = NA, nrow = years + 1, ncol = 8)
    colnames(Tr) <- c("year","death","treat","clinical","undulation","subclinical","minimal","recover")
    for (y in seq(1,years + 1)) {
        Tr[y,"year"] = y - 1
        Tr[y,"death"] = length(which(trajectories[y,] == "death"))
        Tr[y,"treat"] = length(which(trajectories[y,] == "treat"))
        Tr[y,"clinical"] = length(which(trajectories[y,] == "clinical"))
        Tr[y,"undulation"] = length(which(trajectories[y,] == "undulation"))
        Tr[y,"subclinical"] = length(which(trajectories[y,] == "subclinical"))
        Tr[y,"minimal"] = length(which(trajectories[y,] == "minimal"))
        Tr[y,"recover"] = length(which(trajectories[y,] == "recover"))
    }
    
    return(Tr)
}

median_trajectories <- function(nreps, nsteps, npeople, prop_sub, params, method ,treatment, years, prop_min = 0, threshold_months = 9, threshold_changes = 2) {
    Tr <- matrix(data = 0, nrow = years + 1, ncol = 8)
    colnames(Tr) <- c("year","death","treat","clinical","undulation","subclinical","minimal","recover")
    
    for(i in seq(1,nreps)){
        IBM <- full_simulation(nsteps, npeople, prop_sub, params, method, treatment)
        traj <- report_trajectories(IBM, npeople, years, threshold_months, threshold_changes)
        Tr <- Tr + summarise_trajectories(traj,years)
    }
    
    Tr <- Tr/nreps
    
    return(Tr)
    
    
}

plot_trajectories <- function(summary, years, title) {
    pal <- c("death" = '#C0C0C0',
             "clinical" = '#177E89',#'#F8766D',
             "undulation" = '#AE0D0A',#'#7CAE00',
             "subclinical" = '#754668',#'#00BFC4',
             "minimal" = '#CBA715',#'#C77CFF',
             "recovery" = '#A9A9A9')
    summary_melt <- melt(as.data.frame(summary), id.vars = c("year"), variable.name = "state", value.name = "freq")
    p <- ggplot(data = summary_melt, aes(x = year, y = freq, fill = state)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(labels = c("Death","Clinical",
                                     "Transitional","Subclinical",
                                     "Minimal","Recovery"),
                          values = pal) +
        scale_x_continuous("Year", labels = as.character(c(0,1,2,3,4,5)), breaks = c(0,1,2,3,4,5)) +
        scale_y_continuous("Percentage", labels = as.character(c(0,10,20,30,40,50,60,70,80,90,100)),
                           breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)) +
        ggtitle(title)
    
    return(p)
}

example_trajectories <- function(full_sim, traj, state, year) {
    state_sim <- full_sim[,which(traj[year + 1,] == state)]
    example <- sample(seq(1,ncol(state_sim)),1)
    sim <- state_sim[,example]
    
    
    pal <- c("d" = '#696969',
             "t" = '#C0C0C0',
             "c" = '#F8766D',
             "undulation" = '#7CAE00',
             "s" = '#00BFC4',
             "m" = '#C77CFF',
             "e" = '#A9A9A9')
    sim <- as.data.frame(matrix(data = sim, nrow = length(sim)))
    names(sim) <- c("state")
    sim$time <- seq(0,nrow(sim) - 1)
    sim$value <- 0
    sim$value[which(sim$state == "m")] <- 30
    sim$value[which(sim$state == "c")] <- 70
    sim$value[which(sim$state == "s")] <- 50
    sim$value[which(sim$state == "t")] <- 90
    sim$value[which(sim$state == "d")] <- 90
    sim$value[which(sim$state == "e")] <- 10
    sim$value[which(sim$state == "-")] <- NA
    
    rect_data <- data.frame(starty = c(100,80,60,40,20), endy = c(80,60,40,20,0),
                            state = c("t","c","s","m","e"))
    
    sim$state <- factor(sim$state, levels = c("d","t","c","s","m","e"))
    rect_data$state <- factor(rect_data$state, levels = c("t","c","s","m","e"))
    
    
    year_data <- data.frame(startx = c(((year-1) * 12) + 1), endx = c((year * 12)),
                            starty = c(0), endy = c(100))
    
    p <- ggplot() +
        geom_rect(data = rect_data, aes(xmin = 0, xmax = 60, ymin = endy, ymax = starty, fill = state), alpha = 0.5) +
        geom_line(data = sim, aes(x = time, y = value), size = 2) +
        geom_rect(data = year_data, aes(xmin = startx, xmax = endx, ymin = starty, ymax = endy), fill = NA, colour =  "grey42") +
        xlim(c(-1,60)) +
        theme_minimal() +
        scale_fill_manual(values = pal,
                          na.translate = FALSE,
                          breaks = c("t","c","s","m","e"),
                          labels = c("death/treatment","clinical","subclinical","minimal","recovery")) +
        scale_colour_manual(values = pal,
                            na.translate = FALSE) +
        scale_y_continuous(breaks = c(10,30,50,70,90),
                           labels = c("recovery","minimal","subclinical","clinical","death/treatment"),
                           limits = c(0,100)) +
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        xlab("Time (months)") +
        ylab("Disease state") +
        theme(aspect.ratio = 1/2.5) +
        annotate(
            geom = "text",
            label = paste("Year ", year, sep = ""),
            x = ((year-1) * 12) + 6.5,
            y = 90,
            colour = "grey42"
        )
    
    return(p)
    
}



### ---Sankey------------------------------------------------------------------------ ###
sankey_treat <- function(sim) {
    start_states <- c("Clinical","Subclinical")
    end_states <- c("Death","Treatment","Clinical","Subclinical","Minimal","Recovery")
    
    I <- data.frame(id = rep(seq(1,12),2),time_point = rep(c("Start","End"),each = 12),
                    Start = rep(rep(start_states,each = 6),2), end = rep(rep(end_states, 2),2), 
                    state = c(rep(start_states,each = 6),rep(end_states, 2)),freq = rep(0, 24))
    
    I[1,6] <-  length(which(sim[,which(sim[1,] == "c")] == "d"))
    I[7,6] <- length(which(sim[,which(sim[1,] == "s")] == "d"))
    I[2,6] <- length(which(sim[,which(sim[1,] == "c")] == "t"))
    I[8,6] <- length(which(sim[,which(sim[1,] == "s")] == "t"))
    I[6,6] <- length(which(sim[,which(sim[1,] == "c")] == "r"))
    I[12,6] <- length(which(sim[,which(sim[1,] == "s")] == "r"))
    states_from_clinical <- table(sim[121,which(sim[1,] == "c")])
    states_from_subclinical <- table(sim[121,which(sim[1,] == "s")])
    
    I[3,6] <- states_from_clinical[["c"]]
    I[4,6] <- states_from_clinical[["s"]]
    I[5,6] <- states_from_clinical[["m"]]
    
    I[9,6] <- states_from_subclinical[["c"]]
    I[10,6] <- states_from_subclinical[["s"]]
    I[11,6] <- states_from_subclinical[["m"]]
    
    I[13:24,6] <- I[1:12,6]
    
    I$time_point = factor(I$time_point, levels = c("Start","End"))
    I$state = factor(I$state, levels = c("Death","Treatment","Clinical","Subclinical","Minimal","Recovery"))
    I$Start = factor(I$Start, levels = c("Death","Treatment","Clinical","Subclinical","Minimal","Recovery"))
    I$end = factor(I$end, levels = c("Death","Treatment","Clinical","Subclinical","Minimal","Recovery"))
    pal <- c("Death" = '#696969',
             "Treatment" = '#C0C0C0',
             "Clinical" = '#F8766D',
             "undulation" = '#7CAE00',
             "Subclinical" = '#00BFC4',
             "Minimal" = '#C77CFF',
             "Recovery" = '#A9A9A9')
    p <- ggplot(data = I, aes(x = time_point, split = state, value = freq, id = id)) +
        geom_parallel_sets(alpha = 0.7, aes(fill = Start), axis.width = 0.2) +
        geom_parallel_sets_axes(axis.width = 0.2, fill = "light grey") +
        geom_parallel_sets_labels(angle = 0) +
        scale_fill_manual(values = pal,
                          na.translate = FALSE) +
        scale_colour_manual(values = pal,
                            na.translate = FALSE) +
        theme_minimal(base_size = 14) +
        theme(axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank())
    
}

sankey_nat <- function(sim) {
    start_states <- c("Clinical","Subclinical")
    end_states <- c("Death","Clinical","Subclinical","Minimal","Recovery")
    
    I <- data.frame("id" = rep(seq(1,10),2),"time_point" = rep(c("Start","End"),each = 10),
                    "Start" = rep(rep(start_states,each = 5),2), "End" = rep(rep(end_states, 2),2), 
                    "state" = c(rep(start_states,each = 5),rep(end_states, 2)),freq = rep(0, 20))
    
    I[1,6] <-  length(which(sim[1:61,which(sim[1,] == "c")] == "d"))
    I[6,6] <- length(which(sim[1:61,which(sim[1,] == "s")] == "d"))
    I[5,6] <- length(which(sim[1:61,which(sim[1,] == "c")] == "r"))
    I[10,6] <- length(which(sim[1:61,which(sim[1,] == "s")] == "r"))
    states_from_clinical <- table(sim[61,which(sim[1,] == "c")])
    states_from_subclinical <- table(sim[61,which(sim[1,] == "s")])
    
    I[2,6] <- states_from_clinical[["c"]]
    I[3,6] <- states_from_clinical[["s"]]
    I[4,6] <- states_from_clinical[["m"]]
    
    I[7,6] <- states_from_subclinical[["c"]]
    I[8,6] <- states_from_subclinical[["s"]]
    I[9,6] <- states_from_subclinical[["m"]]
    
    I[11:20,6] <- I[1:10,6]
    
    I$time_point = factor(I$time_point, levels = c("Start","End"))
    I$state = factor(I$state, levels = c("Death","Clinical","Subclinical","Minimal","Recovery"))
    I$Start = factor(I$Start, levels = c("Death","Clinical","Subclinical","Minimal","Recovery"))
    I$End = factor(I$End, levels = c("Death","Clinical","Subclinical","Minimal","Recovery"))
    pal <- c("Death" = '#C0C0C0',
             "Treatment" = '#C0C0C0',
             "Clinical" = '#177E89',#'#F8766D',
             "undulation" = '#AE0D0A',#'#7CAE00',
             "Subclinical" = '#754668',#'#00BFC4',
             "Minimal" = '#CBA715',#'#C77CFF',
             "Recovery" = '#A9A9A9')
    I$time_point <- as.numeric(I$time_point)
    p <- ggplot(data = I, aes(x = time_point, split = state, value = freq, id = id)) +
        geom_parallel_sets(alpha = 0.7, aes(fill = Start), axis.width = 0.2) +
        geom_parallel_sets_axes(axis.width = 0.2, fill = "light grey") +
        geom_parallel_sets_labels(angle = 0) +
        scale_fill_manual(values = pal,
                          na.translate = FALSE) +
        scale_colour_manual(values = pal,
                            na.translate = FALSE) +
        scale_x_continuous(breaks = c(1,2), labels = c("Start","End")) +
        theme_minimal(base_size = 14) +
        theme(axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none")
    return(p)
    
}

sankey_min <- function(sim) {
    start_states <- c("Clinical","Subclinical","Minimal")
    end_states <- c("Death","Clinical","Subclinical","Minimal","Recovery")
    
    I <- data.frame("id" = rep(seq(1,15),2),"time_point" = rep(c("Start","End"),each = 15),
                    "Start" = rep(rep(start_states,each = 5),2), "End" = rep(rep(end_states, 3),2), 
                    "state" = c(rep(start_states,each = 5),rep(end_states, 3)),freq = rep(0, 30))
    
    I[1,6] <-  length(which(sim[1:61,which(sim[1,] == "c")] == "d"))
    I[6,6] <- length(which(sim[1:61,which(sim[1,] == "s")] == "d"))
    I[11,6] <- length(which(sim[1:61,which(sim[1,] == "m")] == "d"))
    I[5,6] <- length(which(sim[1:61,which(sim[1,] == "c")] == "r"))
    I[10,6] <- length(which(sim[1:61,which(sim[1,] == "s")] == "r"))
    I[15,6] <- length(which(sim[1:61,which(sim[1,] == "m")] == "r"))
    
    states_from_clinical <- table(sim[61,which(sim[1,] == "c")])
    states_from_subclinical <- table(sim[61,which(sim[1,] == "s")])
    states_from_minimal <- table(sim[61,which(sim[1,] == "m")])
    
    I[2,6] <- states_from_clinical[["c"]]
    I[3,6] <- states_from_clinical[["s"]]
    I[4,6] <- states_from_clinical[["m"]]
    
    I[7,6] <- states_from_subclinical[["c"]]
    I[8,6] <- states_from_subclinical[["s"]]
    I[9,6] <- states_from_subclinical[["m"]]
    
    I[12,6] <- states_from_minimal[["c"]]
    I[13,6] <- states_from_minimal[["s"]]
    I[14,6] <- states_from_minimal[["m"]]
    
    
    I[16:30,6] <- I[1:15,6]
    
    I$time_point = factor(I$time_point, levels = c("Start","End"))
    I$state = factor(I$state, levels = c("Death","Clinical","Subclinical","Minimal","Recovery"))
    I$Start = factor(I$Start, levels = c("Death","Clinical","Subclinical","Minimal","Recovery"))
    I$End = factor(I$End, levels = c("Death","Clinical","Subclinical","Minimal","Recovery"))
    pal <- c("Death" = '#C0C0C0',
             "Treatment" = '#C0C0C0',
             "Clinical" = '#177E89',#'#F8766D',
             "undulation" = '#AE0D0A',#'#7CAE00',
             "Subclinical" = '#754668',#'#00BFC4',
             "Minimal" = '#CBA715',#'#C77CFF',
             "Recovery" = '#A9A9A9')
    I$time_point <- as.numeric(I$time_point)
    p <- ggplot(data = I, aes(x = time_point, split = state, value = freq, id = id)) +
        geom_parallel_sets(alpha = 0.7, aes(fill = Start), axis.width = 0.2) +
        geom_parallel_sets_axes(axis.width = 0.2, fill = "light grey") +
        geom_parallel_sets_labels(angle = 0) +
        scale_fill_manual(values = pal,
                          na.translate = FALSE) +
        scale_colour_manual(values = pal,
                            na.translate = FALSE) +
        scale_x_continuous(breaks = c(1,2), labels = c("Start","End")) +
        theme_minimal(base_size = 14) +
        theme(axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none")
    return(p)
    
}

### ---Rep IBM----------------------------------------------------------------------- ###
rep_IBM_results <- function(n_reps, n_steps, n_people, prop_sub, params, method, treatment = 0, prop_min = 0) {
    if (treatment == 0) {
        results <- as.data.frame(matrix(data = NA, nrow = n_reps, ncol = 44))
        colnames(results) <- c("duration_disease",
                               "duration_infectious",
                               "symp_before_death_min",
                               "symp_before_death_med",
                               "symp_before_death_max",
                               "symp_before_recover_min",
                               "symp_before_recover_med",
                               "symp_before_recover_max",
                               "d1","d2","d3","d4","d5","d10",
                               "c1","c2","c3","c4","c5",
                               "u1","u2","u3","u4","u5",
                               "s1","s2","s3","s4","s5",
                               "m1","m2","m3","m4","m5",
                               "r1","r2","r3","r4","r5",
                               "m5pc","s5pc","u5pc","sm5pc","nc5")
    } else {
        results <- as.data.frame(matrix(data = NA, nrow = n_reps, ncol = 52))
        colnames(results) <- c("duration_disease",
                               "duration_infectious",
                               "symp_before_death_min",
                               "symp_before_death_med",
                               "symp_before_death_max",
                               "symp_before_treat_min",
                               "symp_before_treat_med",
                               "symp_before_treat_max",
                               "symp_before_recover_min",
                               "symp_before_recover_med",
                               "symp_before_recover_max",
                               "d1","d2","d3","d4","d5","d10",
                               "t1","t2","t3","t4","t5",
                               "c1","c2","c3","c4","c5",
                               "u1","u2","u3","u4","u5",
                               "s1","s2","s3","s4","s5",
                               "m1","m2","m3","m4","m5",
                               "r1","r2","r3","r4","r5",
                               "m5pc","s5pc","u5pc","sm5pc","nc5")  
    }
    
    for (i in seq(1,n_reps)) {
        temp <- full_simulation(n_steps, n_people, prop_sub, params, method, treatment, prop_min)
        temp_sum <- summarise_simulation(n_steps, temp)
        traj <- report_trajectories(temp,n_people,5,9,2)
        traj_sum <- summarise_trajectories(traj,5)
        traj1 <- traj_sum[2,]
        traj2 <- traj_sum[3,]
        traj3 <- traj_sum[4,]
        traj4 <- traj_sum[5,]
        traj5 <- traj_sum[6,]
        
        results$duration_disease[i] <- find_median_disease(temp_sum, n_people)
        results$duration_infectious[i] <- find_median_infectious(temp_sum, n_people)
        results$symp_before_death_min[i] <- min(time_clinical_before_event(temp,"d"))
        results$symp_before_death_med[i] <- median(time_clinical_before_event(temp,"d"))
        results$symp_before_death_max[i] <- max(time_clinical_before_event(temp,"d"))
        results$symp_before_recover_min[i] <- min(time_clinical_before_recovery(temp))
        results$symp_before_recover_med[i] <- median(time_clinical_before_recovery(temp))
        results$symp_before_recover_max[i] <- max(time_clinical_before_recovery(temp))
        
        results$d1[i] <- traj1[["death"]]*100/n_people
        results$d2[i] <- traj2[["death"]]*100/n_people
        results$d3[i] <- traj3[["death"]]*100/n_people
        results$d4[i] <- traj4[["death"]]*100/n_people
        results$d5[i] <- traj5[["death"]]*100/n_people
        results$d10[i] <- length(which(temp == "d"))*100/n_people
        
        results$c1[i] <- traj1[["clinical"]]*100/n_people
        results$c2[i] <- traj2[["clinical"]]*100/n_people
        results$c3[i] <- traj3[["clinical"]]*100/n_people
        results$c4[i] <- traj4[["clinical"]]*100/n_people
        results$c5[i] <- traj5[["clinical"]]*100/n_people
        
        results$u1[i] <- traj1[["undulation"]]*100/n_people
        results$u2[i] <- traj2[["undulation"]]*100/n_people
        results$u3[i] <- traj3[["undulation"]]*100/n_people
        results$u4[i] <- traj4[["undulation"]]*100/n_people
        results$u5[i] <- traj5[["undulation"]]*100/n_people
        
        results$s1[i] <- traj1[["subclinical"]]*100/n_people
        results$s2[i] <- traj2[["subclinical"]]*100/n_people
        results$s3[i] <- traj3[["subclinical"]]*100/n_people
        results$s4[i] <- traj4[["subclinical"]]*100/n_people
        results$s5[i] <- traj5[["subclinical"]]*100/n_people
        
        results$m1[i] <- traj1[["minimal"]]*100/n_people
        results$m2[i] <- traj2[["minimal"]]*100/n_people
        results$m3[i] <- traj3[["minimal"]]*100/n_people
        results$m4[i] <- traj4[["minimal"]]*100/n_people
        results$m5[i] <- traj5[["minimal"]]*100/n_people
        
        results$r1[i] <- traj1[["recover"]]*100/n_people
        results$r2[i] <- traj2[["recover"]]*100/n_people
        results$r3[i] <- traj3[["recover"]]*100/n_people
        results$r4[i] <- traj4[["recover"]]*100/n_people
        results$r5[i] <- traj5[["recover"]]*100/n_people
        
        und5 <- temp[,which(traj[6,] == "undulation")]
        sub5 <- temp[,which(traj[6,] == "subclinical")]
        min5 <- temp[,which(traj[6,] == "minimal")]
        rec5 <- temp[,which(traj[6,] == "recover")]
        
        results$u5pc[i] <- length(which(apply(und5,2,function(x) "c" %in% x)))*100/ncol(und5)
        results$m5pc[i] <- length(which(apply(min5,2,function(x) "c" %in% x)))*100/ncol(min5)
        results$r5pc[i] <- length(which(apply(rec5,2,function(x) "c" %in% x)))*100/ncol(rec5)
        if(length(sub5) == 0) {
            results$s5pc[i] <- 0
            results$sm5pc[i] <- length(which(apply(min5,1,function(x) "c" %in% x))) * 100 / (ncol(sub5) + ncol(min5))
        } else if (length(sub5) == n_steps) {
            results$s5pc[i] <- length(which(sub5 == "c"))
            results$sm5pc[i] <- length(which(sub5 == "c")) + length(which(apply(min5,1,function(x) "c" %in% x))) * 100 / (1 + ncol(min5))
        } else {
            results$s5pc[i] <- length(which(apply(sub5,2,function(x) "c" %in% x)))*100/ncol(sub5)
            results$sm5pc[i] <- length(which(apply(sub5,2,function(x) "c" %in% x))) + length(which(apply(min5,1,function(x) "c" %in% x))) * 100 / (ncol(sub5) + ncol(min5))
        }

        
        results$nc5[i] <- 100 - length(which(apply(temp,2,function(x) "c" %in% x)))*100/n_people

        if (treatment != 0) {
            results$symp_before_treat_min[i] <- min(time_clinical_before_event(temp,"t"))
            results$symp_before_treat_med[i] <- median(time_clinical_before_event(temp,"t"))
            results$symp_before_treat_max[i] <- max(time_clinical_before_event(temp,"t"))
            
            results$t1[i] <- traj1[["treat"]]*100/n_people
            results$t2[i] <- traj2[["treat"]]*100/n_people
            results$t3[i] <- traj3[["treat"]]*100/n_people
            results$t4[i] <- traj4[["treat"]]*100/n_people
            results$t5[i] <- traj5[["treat"]]*100/n_people
        }
        
    }
    results_summary <- apply(results, 2, function(x) CI_min_max(x))
    colnames(results_summary) <- colnames(results)
    
    return(results_summary)

}

CI_min_max <- function(data){
    ci <- unname(CI(data, ci = 0.95))
    ci <- unname(quantile(data,c(0.025,0.5,0.975)))
    return(c(ci[3], ci[2],ci[1], min(data), max(data)))
}
