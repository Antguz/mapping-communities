#' @description Spatial block cross-validation version of model_CWM(). It keeps
#' the same 85% subsampling per iteration and the same output columns
#' (performance_batch() reads coefficients by position), but it replaces the
#' random 10-fold CV with a k-fold spatial block CV on the training samples.
#' Out-of-fold predictions are pooled across folds and reported as dataset == "spatial_cv".
#' With data_split = NULL all plots are used for training and the testing row is skipped.
#'
#' @param ordination data.table with PLT_CN and scaled predictors (same order as aggregation)
#' @param aggregation data.table with PLT_CN and the CWM per plant trait
#' @param data_split matrix from data_split_CWM() (1 = training, 0 = testing), or NULL to use all plots
#' @param blocks vector of spatial block IDs per plot (same order as aggregation)
#' @param variables "axes", "climate-topography", or "altogether"
#' @param k number of spatial folds
#' @param n_tries random block-to-fold assignments evaluated to balance samples and trait means across folds
#' @param iterations number of repeated models
#' @param threads number of cores
#'
#' @note cwm_subsplit() and parameters() are copied at the end of this file, so it
#' does not depend on any other script

model_CWM_spatial <- function(ordination,
                              aggregation,
                              data_split = NULL,
                              blocks,
                              variables,
                              k = 5,
                              n_tries = 100,
                              iterations = 100,
                              threads = 25) {

  # Get plant traits
  group_names <- colnames(aggregation)[-1]

  # Get model application in batch
  apply_CWM <- function(X,
                        ordination,
                        aggregation,
                        blocks,
                        variables,
                        k = 5,
                        n_tries = 100,
                        iterations = 100,
                        data_split,
                        group_names) {

    #---------------------------------------------------------------------------
    # Match databases
    if(all.equal(aggregation$PLT_CN, ordination$PLT_CN) == FALSE) {
      stop("Files most match")
    }

    # Group name for future reference
    group_name <- group_names[X]

    #---------------------------------------------------------------------------
    # Reshape data
    frame <- cbind(aggregation[, .SD, .SDcols = group_name], ordination[, .SD, .SDcols = c(2:ncol(ordination))])
    colnames(frame) <- c("trait", colnames(ordination[, .SD, .SDcols = c(2:ncol(ordination))]))

    #---------------------------------------------------------------------------
    # Define variables

    if(variables == "axes") {
      coi <- c("trait", "Axis1", "Axis2", "Axis3")

    } else if(variables == "climate-topography") {
      coi <- c("trait", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")

    } else if(variables == "altogether") {
      coi <- c("trait", "Axis1", "Axis2", "Axis3", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
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
    cwm_collector <- data.table()

    for(ii in 1:iterations) {

      #-------------------------------------------------------------------------
      # Down sampling (same subsamples as model_CWM)
      sub_split <- cwm_subsplit(trait = data_training$trait,
                                split = 0.85,
                                seed = ii)
      data_iteration <- data_training[sub_split,]
      blocks_iteration <- blocks_training[sub_split]

      #-------------------------------------------------------------------------
      # Spatial folds balanced by samples and trait means
      folds <- spatial_folds_CWM(blocks = blocks_iteration,
                                 trait = data_iteration$trait,
                                 k = k,
                                 n_tries = n_tries,
                                 seed = ii)

      #-------------------------------------------------------------------------
      # Spatial cross-validation (pooled out-of-fold predictions)

      oof_pred <- rep(NA_real_, nrow(data_iteration))

      for(f in 1:k) {

        fold_fit <- glm(trait ~ .,
                        data = data_iteration[folds != f,],
                        family = gaussian)
        oof_pred[folds == f] <- predict(fold_fit,
                                        newdata = data_iteration[folds == f,])

      }

      spatial_cv_performance <- parameters_CWM(data_iteration$trait, oof_pred)

      #-------------------------------------------------------------------------
      # Perform model on the full subsample
      glmFit <- glm(trait ~ .,
                    data = data_iteration,
                    family = gaussian)

      #-------------------------------------------------------------------------
      # Performance on training

      training_performance <- parameters_CWM(data_training$trait,
                                             predict(glmFit, newdata = data_training))

      collect <- rbind(training_performance, spatial_cv_performance)
      datasets <- c("training", "spatial_cv")
      psamples <- c(nrow(data_training)/nrow(frame)*100,
                    nrow(data_iteration)/nrow(frame)*100)

      #-------------------------------------------------------------------------
      # Performance on testing (only with data_split)

      if(!is.null(data_split)) {

        testing_performance <- parameters_CWM(data_testing$trait,
                                              predict(glmFit, newdata = data_testing))

        collect <- rbind(training_performance, testing_performance, spatial_cv_performance)
        datasets <- c("training", "testing", "spatial_cv")
        psamples <- c(psamples[1],
                      nrow(data_testing)/nrow(frame)*100,
                      psamples[2])

      }

      rows <- length(datasets)

      #-------------------------------------------------------------------------
      # Collect data
      # Coefficients
      coefficients <- matrix(glmFit$coefficients, nrow = rows, ncol = length(glmFit$coefficients), byrow = TRUE)
      colnames(coefficients) <- names(glmFit$coefficients)
      collect <- cbind(collect, coefficients)

      # varImp
      var <- as.matrix(t(varImp(glmFit)), nrow = 1)
      var <- var[rep(1, rows), , drop = FALSE]
      colnames(var) <- paste0("imp_", colnames(var))
      collect <- cbind(collect, var)

      # Collect
      collect <- cbind(plant = rep(group_name, rows),
                       iteration = rep(ii, rows),
                       dataset = datasets,
                       psamples = psamples,
                       collect)

      # Extra columns at the end to keep coefficient positions of model_CWM
      collect$min_fold_samples <- c(rep(NA, rows - 1), min(tabulate(folds, nbins = k)))
      collect$n_blocks <- c(rep(NA, rows - 1), length(unique(blocks_iteration)))

      cwm_collector <- rbind(cwm_collector, collect)

    }

    return(cwm_collector)

  }

  # Parallel processing
  complete <- pbmclapply(X = 1:length(group_names),
                         FUN = apply_CWM,
                         ordination = ordination,
                         aggregation = aggregation,
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
    stop(paste0("Plant trait ", group_names[which(errors)[1]], ": ", complete[[which(errors)[1]]]))
  }

  complete <- do.call(rbind, complete)

  return(complete)

}

# Assign spatial blocks to k folds, keeping the most balanced assignment of
# samples and trait means across folds out of n_tries random assignments
spatial_folds_CWM <- function(blocks,
                              trait,
                              k = 5,
                              n_tries = 100,
                              seed) {

  unique_blocks <- unique(blocks)

  if(length(unique_blocks) < k) {
    stop(paste0("Only ", length(unique_blocks), " spatial blocks for k = ", k, " folds; check block_size"))
  }
  block_index <- match(blocks, unique_blocks)
  block_samples <- tabulate(block_index, nbins = length(unique_blocks))
  block_trait <- as.numeric(tapply(trait, factor(block_index, levels = 1:length(unique_blocks)), sum))

  best_score <- Inf
  best_folds <- NULL

  for(t in 1:n_tries) {

    set.seed(seed*1000 + t)
    block_fold <- sample(rep(1:k, length.out = length(unique_blocks)))

    fold_samples <- tapply(block_samples, factor(block_fold, levels = 1:k), sum)
    fold_mean <- tapply(block_trait, factor(block_fold, levels = 1:k), sum)/fold_samples

    # Coefficient of variation of samples and spread of trait means (in SD units) across folds
    score <- sd(fold_samples)/mean(fold_samples) + sd(fold_mean)/sd(trait)

    if(score < best_score) {
      best_score <- score
      best_folds <- block_fold
    }
  }

  return(best_folds[block_index])

}

# Residuals of a gaussian GLM per plant trait fitted on all plots.
# Used to define the spatial block size
residuals_CWM <- function(aggregation,
                          ordination,
                          variables = "altogether") {

  if(variables == "axes") {
    coi <- c("Axis1", "Axis2", "Axis3")
  } else if(variables == "climate-topography") {
    coi <- c("MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
  } else if(variables == "altogether") {
    coi <- c("Axis1", "Axis2", "Axis3", "MAT", "TS", "TAR", "MAP", "PS", "DEM", "SL", "CTI")
  }

  # Plant traits modeled
  group_names <- colnames(aggregation)[-1]

  residual_collector <- data.table(PLT_CN = aggregation$PLT_CN)

  for(i in 1:length(group_names)) {

    frame <- cbind(data.table(trait = aggregation[[group_names[i]]]),
                   ordination[, ..coi])

    fit <- glm(trait ~ ., data = frame, family = gaussian)
    residual_collector[, (make.names(group_names[i])) := residuals(fit, type = "deviance")]

  }

  return(residual_collector)

}

# Parallel version of blockCV::cv_spatial_autocor(), fitting one variogram per
# column in each core. The range is the median across columns, as in blockCV
autocor_CWM <- function(x,
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

# Stratified subsample along the trait distribution (copied from model_CWM.R)
cwm_subsplit <- function(trait,
                         split = 0.85,
                         seed) {

  test <- data.table(sample = 1:length(trait),
                     trait = trait)
  test$cut <- cut(test$trait, breaks = 10)
  test$train <- 0
  cuts_unite <- unique(test$cut)
  all <- data.frame()

  for(k in 1:length(cuts_unite)) {

    sub <- subset(test, cut == cuts_unite[k])
    set.seed(seed)
    sample_test <- sample(1:nrow(sub), ceiling(nrow(sub)*split))
    sub[sample_test, train := 1]
    all <- rbind(all, sub)

  }

  all <- all[order(sample),]
  all <- all[train == 1,]

  return(all$sample)

}

# Performance parameters (copied from parameters() in model_performance.R)
parameters_CWM <- function(obs, pred) {

  linear <- lm(obs ~ pred, na.action = na.exclude)
  linear_summary <- summary(linear)

  R2 <- 1 - sum((obs - pred)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2)
  intercept <- linear_summary$coefficients[1,1]
  slope <- linear_summary$coefficients[2,1]
  BIAS <- mean(obs-pred, na.rm = TRUE)
  RMSE <- sqrt(sum(((obs-pred)^2))/length(obs))
  perRMSE <- sqrt(sum(((obs-pred)^2))/length(obs))/(quantile(obs, 0.99) - quantile(obs, 0.01))
  names(perRMSE) <- NULL

  return(round(data.table(R2 = R2,
                          BIAS = BIAS,
                          RMSE = RMSE,
                          perRMSE = perRMSE*100,
                          intercept = intercept,
                          slope = slope), 4))

}
