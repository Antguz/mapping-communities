#' @description Spatial block cross-validation version of model_fortyp(). It keeps
#' the same 85% subsampling per iteration and the same output columns 
#' (performance_batch() reads coefficients by position), but it replaces the
#' random 10-fold CV with a k-fold spatial block CV on the training samples.
#' Out-of-fold predictions are pooled across folds and reported as dataset == "spatial_cv".
#' With data_split = NULL all plots are used for training and the testing row is skipped.
#' Class imbalance is handled with ROSE inside each training fold and on the final
#' subsample (seed = ii, as in model_fortyp), never on the held-out blocks.
#'
#' @param forest_type data.table with PLT_CN and FORTYP
#' @param ordination data.table with PLT_CN and scaled predictors (same order as forest_type)
#' @param data_split matrix from data_split_fortyp() (1 = training, 0 = testing), or NULL to use all plots
#' @param blocks vector of spatial block IDs per plot (same order as forest_type)
#' @param variables "axes", "climate-topography", or "altogether"
#' @param k number of spatial folds
#' @param n_tries random block-to-fold assignments evaluated to balance presences across folds
#' @param iterations number of repeated models
#' @param threads number of cores
#'
#' @note split_fortyp() is included at the end of this file (same as model_fortyp.R)

model_fortyp_spatial <- function(forest_type,
                                 ordination,
                                 data_split = NULL,
                                 blocks,
                                 variables,
                                 k = 5,
                                 n_tries = 100,
                                 iterations = 100,
                                 threads = 25) {

  # Get plant groups (same order as data_split_fortyp)
  if(is.null(data_split)) {
    group_names <- unique(forest_type$FORTYP)
  } else {
    group_names <- colnames(data_split)
  }

  # Get model application in batch
  apply_occurrence <- function(X,
                               forest_type,
                               ordination,
                               blocks,
                               variables,
                               k = 5,
                               n_tries = 100,
                               iterations = 100,
                               data_split,
                               group_names) {

    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(forest_type$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }

    # Group name for future reference
    group_name <- group_names[X]
    ftype <- forest_type
    ftype$occurance <- "absence"
    ftype[FORTYP == group_name, occurance := "presence"]

    if(group_name == "NA" | nrow(ftype[occurance == "presence",]) < 100) {

      return(data.table())

    } else {

      #---------------------------------------------------------------------------
      # Reshape data
      frame <- cbind(ftype[, .SD, .SDcols = 3], ordination[, .SD, .SDcols = c(2:ncol(ordination))])
      colnames(frame) <- c("occurrence", colnames(ordination[, .SD, .SDcols = c(2:ncol(ordination))]))
      frame$occurrence <- as.factor(frame$occurrence)
      frame$occurrence <- factor(frame$occurrence, levels = c("presence", "absence"))

      #---------------------------------------------------------------------------
      # Define variables

      if(variables == "axes") {
        coi <- c("occurrence", "Axis1", "Axis2", "Axis3")

      } else if(variables == "climate-topography") {
        coi <- c("occurrence", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")

      } else if(variables == "altogether") {
        coi <- c("occurrence", "Axis1", "Axis2", "Axis3", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
      }

      frame <- frame[, ..coi]

      #---------------------------------------------------------------------------
      # Training and testing (all plots for training when data_split = NULL)

      if(is.null(data_split)) {
        training_rows <- rep(TRUE, nrow(frame))
      } else {
        training_rows <- as.matrix(data_split)[, X] == 1
      }

      data_training <- frame[training_rows,]
      data_testing <- frame[!training_rows,]
      blocks_training <- blocks[training_rows]
      fortyp_collector <- data.table()

      for(ii in 1:iterations) {

        #---------------------------------------------------------------------------
        # Down sampling (same subsamples as model_fortyp)
        split <- split_fortyp(data_training$occurrence, split = 0.85, seed = ii)
        data_iteration <- data_training[split,]
        blocks_iteration <- blocks_training[split]

        #---------------------------------------------------------------------------
        # Spatial folds balanced by presences
        folds <- spatial_folds(blocks = blocks_iteration,
                               occurrence = data_iteration$occurrence,
                               k = k,
                               n_tries = n_tries,
                               seed = ii)

        #---------------------------------------------------------------------------
        # Spatial cross-validation (pooled out-of-fold predictions)

        oof_presence <- rep(NA_real_, nrow(data_iteration))

        for(f in 1:k) {

          fold_train <- data_iteration[folds != f,]

          if(sum(fold_train$occurrence == "presence") == 0) {
            next
          }

          fold_fit <- fit_fortyp(fold_train, seed = ii*100 + f)
          oof_presence[folds == f] <- 1 - predict(fold_fit,
                                                  newdata = data_iteration[folds == f,],
                                                  type = "response")

        }

        spatial_cv <- probability_frame(obs = data_iteration$occurrence[!is.na(oof_presence)],
                                        presence = oof_presence[!is.na(oof_presence)])

        spatial_cv_performance <- performance_fortyp(spatial_cv)

        #---------------------------------------------------------------------------
        # Perform model on the full subsample
        glmFit <- fit_fortyp(data_iteration, seed = ii)

        #---------------------------------------------------------------------------
        # Performance on training

        training <- probability_frame(obs = data_training$occurrence,
                                      presence = 1 - predict(glmFit,
                                                             newdata = data_training,
                                                             type = "response"))

        training_performance <- performance_fortyp(training)

        #---------------------------------------------------------------------------
        # Collect data
        training <- as.data.table(training)
        spatial_cv <- as.data.table(spatial_cv)

        # Performance and coefficients
        collect <- rbind(c(training_performance, glmFit$coefficients),
                         c(spatial_cv_performance, glmFit$coefficients))
        datasets <- c("training", "spatial_cv")
        psamples <- c(nrow(training[obs == "presence",])/nrow(training)*100,
                      nrow(spatial_cv[obs == "presence",])/nrow(spatial_cv)*100)

        #---------------------------------------------------------------------------
        # Performance on testing (only with data_split)

        if(!is.null(data_split)) {

          testing <- probability_frame(obs = data_testing$occurrence,
                                       presence = 1 - predict(glmFit,
                                                              newdata = data_testing,
                                                              type = "response"))

          testing_performance <- performance_fortyp(testing)
          testing <- as.data.table(testing)

          collect <- rbind(collect[1,], 
                           c(testing_performance, glmFit$coefficients),
                           collect[2,])
          datasets <- c("training", "testing", "spatial_cv")
          psamples <- c(psamples[1],
                        nrow(testing[obs == "presence",])/nrow(testing)*100,
                        psamples[2])

        }

        collect <- as.data.table(collect)
        rows <- length(datasets)

        # varImp
        var <- as.matrix(t(varImp(glmFit)), nrow = 1)
        var <- var[rep(1, rows), , drop = FALSE]
        colnames(var) <- paste0("imp_", colnames(var))
        collect <- cbind(collect, var)

        # Merge
        collect <- cbind(forty = rep(group_name, rows),
                         iteration = rep(ii, rows),
                         dataset = datasets,
                         psamples = psamples,
                         collect)

        # Extra columns at the end to keep coefficient positions of model_fortyp
        collect$TSS <- collect$Sens + collect$Spec - 1
        collect$min_fold_presence <- c(rep(NA, rows - 1), min(table(factor(folds, levels = 1:k)[data_iteration$occurrence == "presence"])))
        collect$n_blocks <- c(rep(NA, rows - 1), length(unique(blocks_iteration)))

        fortyp_collector <- rbind(fortyp_collector, collect)

      }

      return(fortyp_collector)

    }
  }

  # Parallel processing
  complete <- pbmclapply(X = 1:length(group_names),
                         FUN = apply_occurrence,
                         forest_type = forest_type,
                         ordination = ordination,
                         blocks = blocks,
                         variables = variables,
                         k = k,
                         n_tries = n_tries,
                         data_split = data_split,
                         group_names = group_names,
                         iterations = iterations,
                         mc.preschedule = TRUE,
                         mc.set.seed = FALSE,
                         mc.silent = FALSE,
                         mc.cores = threads,
                         mc.cleanup = TRUE,
                         mc.allow.recursive = TRUE,
                         affinity.list = NULL)

  # Report the first error from the workers
  errors <- sapply(complete, inherits, "try-error")
  if(any(errors)) {
    stop(paste0("Forest type ", group_names[which(errors)[1]], ": ", complete[[which(errors)[1]]]))
  }
  
  complete <- do.call(rbind, complete)
  complete <- complete[!is.na(dataset),]

  return(complete)

}

# Assign spatial blocks to k folds, keeping the most balanced assignment of
# presences (and samples) across folds out of n_tries random assignments
spatial_folds <- function(blocks,
                          occurrence,
                          k = 5,
                          n_tries = 100,
                          seed) {

  unique_blocks <- unique(blocks)
  
  if(length(unique_blocks) < k) {
    stop(paste0("Only ", length(unique_blocks), " spatial blocks for k = ", k, " folds; check block_size"))
  }
  block_index <- match(blocks, unique_blocks)
  block_presence <- tabulate(block_index[occurrence == "presence"], nbins = length(unique_blocks))
  block_samples <- tabulate(block_index, nbins = length(unique_blocks))

  best_score <- Inf
  best_folds <- NULL

  for(t in 1:n_tries) {

    set.seed(seed*1000 + t)
    block_fold <- sample(rep(1:k, length.out = length(unique_blocks)))

    fold_presence <- tapply(block_presence, factor(block_fold, levels = 1:k), sum)
    fold_samples <- tapply(block_samples, factor(block_fold, levels = 1:k), sum)

    # Coefficient of variation of presences and samples across folds
    score <- sd(fold_presence)/mean(fold_presence) + sd(fold_samples)/mean(fold_samples)

    if(score < best_score) {
      best_score <- score
      best_folds <- block_fold
    }
  }

  return(best_folds[block_index])

}

# Fit a binomial GLM on ROSE-balanced data to deal with class imbalance. ROSE is applied
# only to the data passed here (training folds or the final subsample), so held-out 
# spatial blocks keep the observed class distribution. As in caret, the factor levels
# c("presence", "absence") make glm() model P(absence), so P(presence) = 1 - p
fit_fortyp <- function(data, seed) {
  
  # Class balance (same as model_fortyp.R)
  rose_data <- ROSE(occurrence ~ ., data = data, seed = seed)$data
  rose_data$occurrence <- factor(rose_data$occurrence, levels = c("presence", "absence"))
  
  fit <- glm(occurrence ~ ., 
             data = rose_data, 
             family = binomial)
  
  return(fit)
  
}

# Deviance residuals on the observed plots (residuals() of a ROSE fit would return
# the residuals of the synthetic samples)
deviance_fortyp <- function(fit, data) {
  
  # glm() models P(absence)
  y <- as.numeric(data$occurrence == "absence")
  mu <- predict(fit, newdata = data, type = "response")
  mu <- pmin(pmax(mu, 1e-15), 1 - 1e-15)
  
  sign(y - mu)*sqrt(-2*(y*log(mu) + (1 - y)*log(1 - mu)))
  
}

# Observed-predicted frame as used by caret summary functions (cutoff = 0.5)
probability_frame <- function(obs, presence) {

  data.frame(obs = obs,
             pred = factor(ifelse(presence >= 0.5, "presence", "absence"),
                           levels = c("presence", "absence")),
             presence = presence,
             absence = 1 - presence)

}

# Same metrics and order as model_fortyp
performance_fortyp <- function(frame) {

  lev <- levels(frame$obs)

  c(defaultSummary(frame, lev = lev),
    twoClassSummary(frame, lev = lev),
    prSummary(frame, lev = lev),
    mnLogLoss(frame, lev = lev))

}

# Deviance residuals of a ROSE-balanced GLM per forest type (same filters as 
# model_fortyp_spatial) fitted on all plots. Used to define the spatial block size
residuals_fortyp <- function(forest_type,
                             ordination,
                             variables = "altogether",
                             min_presence = 100) {
  
  if(variables == "axes") {
    coi <- c("Axis1", "Axis2", "Axis3")
  } else if(variables == "climate-topography") {
    coi <- c("MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
  } else if(variables == "altogether") {
    coi <- c("Axis1", "Axis2", "Axis3", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
  }
  
  # Forest types modeled
  counts <- forest_type[, .N, by = FORTYP]
  group_names <- counts[N >= min_presence & FORTYP != "NA", FORTYP]
  
  residual_collector <- data.table(PLT_CN = forest_type$PLT_CN)
  
  for(i in 1:length(group_names)) {
    
    frame <- cbind(data.table(occurrence = factor(ifelse(forest_type$FORTYP == group_names[i], 
                                                         "presence", "absence"),
                                                  levels = c("presence", "absence"))),
                   ordination[, ..coi])
    
    fit <- fit_fortyp(frame, seed = i)
    residual_collector[, (make.names(group_names[i])) := deviance_fortyp(fit, frame)]
    
  }
  
  return(residual_collector)
  
}

# Parallel version of blockCV::cv_spatial_autocor(), fitting one variogram per 
# column in each core. The range is the median across columns, as in blockCV
autocor_fortyp <- function(x,
                           columns,
                           num_sample = 5000,
                           threads = 25) {
  
  apply_autocor <- function(X, x, columns, num_sample) {
    
    set.seed(X)
    autocor <- cv_spatial_autocor(x = x[, columns[X]],
                                  column = columns[X],
                                  num_sample = num_sample,
                                  plot = FALSE,
                                  progress = FALSE)
    
    return(autocor$range_table)
    
  }
  
  # Parallel processing
  complete <- pbmclapply(X = 1:length(columns), 
                         FUN = apply_autocor,
                         x = x,
                         columns = columns,
                         num_sample = num_sample,
                         mc.preschedule = TRUE, 
                         mc.set.seed = FALSE,
                         mc.silent = FALSE, 
                         mc.cores = threads,
                         mc.cleanup = TRUE, 
                         mc.allow.recursive = TRUE, 
                         affinity.list = NULL)
  
  # Report the first error from the workers
  errors <- sapply(complete, inherits, "try-error")
  if(any(errors)) {
    stop(paste0("Column ", columns[which(errors)[1]], ": ", complete[[which(errors)[1]]]))
  }
  
  range_table <- as.data.table(do.call(rbind, complete))
  setorder(range_table, range)
  
  return(list(range = median(range_table$range),
              range_table = range_table))
  
}

# Stratified subsample of presences and absences (copied from model_fortyp.R)
split_fortyp <- function(occurrence,
                         split = 0.85,
                         seed) {
  
  vec <- integer()
  
  split_frame <- data.table(sample = 1:length(occurrence), occurrence = occurrence)
  
  # Presence
  presence <- split_frame[occurrence == "presence",]
  set.seed(seed)
  train_presence <- sample(1:nrow(presence), floor(nrow(presence)*split))
  vec <- c(vec, presence$sample[train_presence])
  
  # Absence
  absence <- split_frame[occurrence == "absence",]
  set.seed(seed)
  train_absence <- sample(1:nrow(absence), floor(nrow(absence)*split))
  vec <- c(vec, absence$sample[train_absence])
  
  return(vec)
  
}
