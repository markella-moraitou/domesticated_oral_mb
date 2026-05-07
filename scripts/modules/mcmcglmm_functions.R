
############################
#### MCMCglmm FUNCTIONS ####
############################

mcmcglmm_wrapper <- function(input_data, formula, output_basename, processes) {
    rds_name <- paste0(output_basename, ".RDS")
    if (file.exists(file.path(subdir, rds_name))) {
        cat("MCMCglmm output exists. Loading...\n")
        m <- readRDS(file.path(subdir, rds_name))
    } else {
        cat("Running MCMCglmm...\n")
        set.seed(14)
        # Prep formula
        n_taxa <- length(unique(input_data$OTU))
        # Set priors
        prior = list(R = list(V = diag(1), nu = 0.002))
        start.time <- Sys.time()
        cat("Starting MCMCglmm... ")
        print(start.time)
        # Run MCMCglmm
        m <- mclapply(1:5, function(i) {
            MCMCglmm(fixed = formula,
            prior = prior,
            data = input_data,
            verbose = TRUE, # pr = TRUE, pl = TRUE,
            nitt = 50000,
            burnin = 20000,
            thin = 50)
        }, mc.cores = processes)
        end.time <- Sys.time()
        cat("Finished MCMCglmm. ")
        print(end.time)
        time.lapsed <- end.time - start.time
        cat("Time lapsed for MCMCglmm: ")
        print(time.lapsed)
        # Save output
        cat("Saving at: ", file.path(subdir, rds_name), "\n")
        saveRDS(m, file.path(subdir, rds_name))
    }
    return(m)
}

process_mcmcglmm_out <- function(mcmcglmm_output, output_basename) {
    if (file.exists(file.path(subdir, paste0(output_basename, "_results.csv")))) {
        mcmc_res <- read.csv(file.path(subdir, paste0(output_basename, "_results.csv")))
    } else {
        m <- mcmcglmm_output
        mlist <- lapply(m, function(model) model$Sol)
        mlist <- do.call(mcmc.list, mlist)
        
        # Diagnostics with gelman plot
        pdf(file=file.path(subdir, paste0(output_basename, "_gelman_plots.pdf")))
        par(mfrow=c(4,2), mar=c(2,2,1,2))
        gelman.plot(mlist, auto.layout=F, bin.width = 1)
        dev.off()
        
        # Plot first chain
        m1 = m[[1]]
        
        par(ask=FALSE)
        pdf(file=file.path(subdir, paste0(output_basename, "_plots.pdf")))
        par(mfrow = c(2,2), ask=FALSE)
        plot(m1)
        dev.off()

        # Autocorrelation
        diag(autocorr(m1$VCV)[2, , ])

        # 95% Credible interval
        HPDinterval(m1$VCV)

        # Collect results into tables
        mcmc_res <- summary(m1)$solutions %>%
            data.frame %>% rownames_to_column("term") %>%
            mutate(term = str_replace(term, "path:", "path_")) %>% # For pathway diff abund
            mutate(OTU = str_extract(term, "OTU[^:]*")) %>% # Separate OTU
            mutate(term = str_remove(term, OTU) %>% str_remove(":")) %>%  # Separate term
            mutate(OTU = str_remove_all(OTU, "OTU|:")) %>% # Remove fluff from OTU name
            mutate(OTU = str_replace(OTU, "path_", "path:")) %>% # For pathway diff abund
            # Remove OTU intercepts
            filter(term != "")

        write.csv(mcmc_res, file = file.path(subdir, paste0(output_basename, "_results.csv")), quote = FALSE, row.names = FALSE)
    }
    return(mcmc_res)
}
