#' Feature Selection with Cox Elastic Net via Stratified Bootstrap Resampling
#'
#' Performs feature selection using the Cox proportional hazards model with Elastic Net regularization
#' through stratified bootstrap resampling. Aggregates feature importance across multiple resampling rounds.
#'
#' @param data A dataframe where first two columns are survival time and status. Subsequent columns are features.
#' @param seed Random seed for reproducibility (default: 123)
#' @param n_round Number of bootstrap resampling rounds (default: 20)
#' @param topn Number of top features to retain based on average importance (default: 10)
#' @param lambda_type Regularization criterion for cv.glmnet: "lambda.min" (default) or "lambda.1se"
#'
#' @return A list containing:
#' \itemize{
#'   \item selected_matrix - Matrix with survival data and selected features
#'   \item importance_table - Dataframe of feature rankings by average importance
#'   \item coef_matrix - Raw coefficient matrix from all resampling rounds
#' }
#'
#' @details
#' Implements stratified bootstrap resampling preserving event/non-event ratios. Uses Elastic Net
#' regularization (alpha = 0.5) via glmnet. Feature importance is calculated as mean absolute coefficients
#' across resampling rounds.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(time = rnorm(100), status = rbinom(100, 1, 0.5),
#'                   matrix(rnorm(100*20), ncol=20))
#' result <- feature_selection_coxlasso(data, topn = 5)
#' }
#'
#' @export
#' @importFrom glmnet cv.glmnet
#' @importFrom stats coef
#'
feature_selection_coxlasso <- function(data, seed=123, n_round=20, topn=10, lambda_type="lambda.min") {
  # Parameter description:
  # data: Dataframe where the first two columns are survival time (time) and survival status (status)
  # seed: Random seed
  # n_round: Number of resampling rounds
  # topn: Number of top features to retain
  # lambda_type: Type of regularization parameter selection ("lambda.min" or "lambda.1se")

  require(glmnet)
  set.seed(seed)  # Initialize random seed

  # Extract survival data
  survival_time <- data[, 1]
  survival_status <- data[, 2]
  X <- as.matrix(data[, -c(1:2), drop = FALSE])

  # Initialize feature importance matrix
  coef_matrix <- matrix(0, nrow = n_round, ncol = ncol(X))
  colnames(coef_matrix) <- colnames(X)

  for (i in 1:n_round) {
    # Stratified Bootstrap resampling (preserving event ratio)
    event_indices <- which(survival_status == 1)
    nonevent_indices <- which(survival_status == 0)

    sampled_event <- sample(event_indices, size = length(event_indices), replace = TRUE)
    sampled_nonevent <- sample(nonevent_indices, size = length(nonevent_indices), replace = TRUE)
    sampled_indices <- c(sampled_event, sampled_nonevent)

    # Prepare resampled data
    X_sample <- X[sampled_indices,,drop = FALSE]
    y_sample <- cbind(time = survival_time[sampled_indices],
                      status = survival_status[sampled_indices])

    # Fit Cox LASSO model
    cv_fit <- cv.glmnet(X_sample, y_sample,
                        family = "cox",
                        alpha = 0.5,  # LASSO regularization
                        nfolds = 5,
                        grouped = TRUE)

    # Extract non-zero coefficient features
    final_coef <- as.vector(coef(cv_fit, s = cv_fit[[lambda_type]]))
    coef_matrix[i, ] <- abs(final_coef)  # Use absolute value as importance metric

    # Display progress
    if(i %% 5 == 0) message("Completed round ", i, "/", n_round, " of resampling")
  }

  # Calculate average feature importance
  avg_coef <- colMeans(coef_matrix)
  ranked_features <- names(sort(avg_coef, decreasing = TRUE))[1:topn]

  #ranked_features = ranked_features[-which(colMeans(coef_matrix)==0)]

  # Construct result data
  selected_data <- cbind(data[, 1:2], X[, ranked_features, drop = FALSE])


  colnames(selected_data)[1:2] <- c("time", "status")

  # Create importance ranking table
  importance_table <- data.frame(
    Feature = names(avg_coef),
    Importance = avg_coef,
    Rank = rank(-avg_coef, ties.method = "min")
  )[order(-avg_coef), ]
  rownames(importance_table) <- NULL

  return(list(
    selected_matrix = selected_data,
    importance_table = importance_table,
    coef_matrix = coef_matrix  # Raw coefficient matrix for further analysis
  ))
}


#' Survival Data Splitting with Stratification
#'
#' Splits survival data into training and testing sets while preserving event/non-event ratios.
#' Supports both stratified and non-stratified splitting methods.
#'
#' @param data A dataframe where the first two columns are survival time and status.
#' @param train_ratio Proportion of data to use for training (default: 0.7)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param stratify Logical indicating whether to use stratified sampling (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item train - Training set as a dataframe
#'   \item test - Testing set as a dataframe
#'   \item report - Stratification report showing event proportions
#' }
#'
#' @details
#' Automatically converts survival time to numeric and ensures proper data types.
#' Uses `rsample` for stratified splitting when possible, falling back to `caret`
#' for non-stratified splitting. Includes safeguards for small sample sizes and
#' invalid data types.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(time = rnorm(100), status = rbinom(100, 1, 0.5),
#'                   matrix(rnorm(100*20), ncol=20))
#' split_data <- surv_data_split(data, train_ratio = 0.8)
#' }
#'
#' @export
#' @importFrom dplyr mutate across where
#' @importFrom rsample initial_split training testing
#' @importFrom caret createDataPartition

surv_data_split <- function(data, train_ratio=0.7, seed = NULL, stratify=TRUE) {
  require(dplyr)

  # Force the survival time column to be numeric
  if(!is.numeric(data[,1])) {
    data[,1] <- as.numeric(data[,1])
    warning("Survival time column has been converted to numeric")
  }

  # Explicitly convert to a standard data frame
  data <- as.data.frame(data) %>%
    mutate(across(where(is.list), as.vector))  # Remove potential unconventional types

  # Rename the first two columns
  colnames(data)[1:2] <- c("time", "status")

  # Add event status factor validation
  data$status <- as.factor(data$status)

  if(stratify) {
    # Ensure minimum sample logic
    min_samples <- ceiling(1/(1 - train_ratio))  # Automatically adjust the minimum sample threshold

    # Group processing
    set.seed(seed)
    grouped_data <- split(data, data$status)

    # Ensure each group has available samples
    if(any(sapply(grouped_data, nrow) < 2)) {
      warning("Some event groups have too few samples; automatically switching to non-stratified sampling")
      stratify <- FALSE
    }
  }

  if(stratify) {
    # Use the rsample package for reliable stratified splitting
    require(rsample)
    set.seed(seed)
    split_obj <- initial_split(data,
                               strata = "status",
                               prop = train_ratio,
                               breaks = 4)  # Optimize stratified binning

    train_set <- training(split_obj)
    test_set <- testing(split_obj)

    # Shuffle the order while maintaining row consistency
    train_set <- train_set[sample(nrow(train_set)), ]
    test_set <- test_set[sample(nrow(test_set)), ]

  } else {
    # Use the caret package to create data splitting
    require(caret)
    set.seed(seed)
    train_idx <- createDataPartition(data$status,
                                     p = train_ratio,
                                     list = FALSE,
                                     times = 1)
    train_set <- data[train_idx, ]
    test_set <- data[-train_idx, ]
  }

  # Check feature data types
  type_check <- sapply(train_set, class)
  if(any(type_check == "list")) {
    stop("There are illegal list-type columns in the data; please use asi_* functions to convert")
  }

  # Output diagnostic report
  stratification_report <- rbind(
    Original = prop.table(table(data$status)),
    Training = prop.table(table(train_set$status)),
    Testing = prop.table(table(test_set$status))
  )

  message("Data splitting completed, stratification proportions:\n")
  print(round(stratification_report, 4))

  return(list(
    train = as.data.frame(train_set),
    test = as.data.frame(test_set),
    report = stratification_report
  ))
}

#' Machine Learning Pipeline for Cox Proportional Hazards Model
#'
#' A comprehensive pipeline for survival analysis using Cox proportional hazards model with Elastic Net regularization.
#' Includes data splitting, feature selection, model training, risk score calculation, and survival curve visualization.
#'
#' @param data A dataframe where the first two columns are survival time and status.
#' @param seed Random seed for reproducibility (default: 123)
#' @param train_percent Proportion of data to use for training (default: 0.7)
#' @param n_round Number of bootstrap rounds for feature selection (default: 20)
#' @param stratify Logical indicating whether to use stratified sampling (default: FALSE)
#' @param lambda_type Regularization criterion for cv.glmnet: "lambda.min" (default) or "lambda.1se"
#' @param bestcut Logical indicating whether to determine optimal cutoff for risk groups (default: FALSE)
#' @param topn Number of top features to retain in feature selection (default: 10)
#'
#' @return A list containing:
#' \itemize{
#'   \item Training set results with risk scores and groups
#'   \item Testing set results with risk scores and groups
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Data preprocessing and splitting
#' 2. Feature selection using Cox LASSO
#' 3. Model training with Elastic Net regularization
#' 4. Risk score calculation and optimal cutoff determination
#' 5. Survival curve visualization and PDF export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(time = rnorm(100), status = rbinom(100, 1, 0.5),
#'                   matrix(rnorm(100*20), ncol=20))
#' results <- machine_learning_cox(data, train_percent = 0.8, topn = 5)
#' }
#'
#' @export
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival Surv survfit
#' @importFrom survminer surv_cutpoint ggsurvplot
#' @importFrom caret trainControl
#' @importFrom patchwork plot_layout
#' @importFrom dplyr mutate

machine_learning_cox = function(data,
                                seed = 123,
                                train_percent = 0.7,
                                n_round = 20,stratify= F,
                                lambda_type="lambda.min",
                                bestcut = F,
                                topn=10){


  colnames(data)[c(1,2)] = c("time","status")


  # When the sample size is small (n<100), you can increase the minimum sample guarantee:
  xx = surv_data_split(data, train_ratio= train_percent,
                       stratify = stratify,
                       seed = seed)

  trainh <- xx$train
  testh  <- xx$test

  trainh =   trainh[,apply(trainh,2,var)!=0]
  testh =   testh[,apply(testh,2,var)!=0]

  trainh[,-c(1,2)] = scale(trainh[,-c(1,2)])
  testh[,-c(1,2)] = scale(testh[,-c(1,2)])

  # x = rfe_topn(data = trainh,topn = topn,model_type = model_type,seed = seed)
  x = feature_selection_coxlasso(data = trainh ,seed=seed,
                                 n_round= n_round, topn=topn,lambda_type =lambda_type)


  trainh = x$selected_matrix
  testh =  testh[,colnames(testh) %in% colnames(trainh)]

  trainh = trainh[,colnames(trainh) %in% colnames(testh)]

  write.csv(x$coef_matrix,"coef_matrix.csv")
  write.csv(x$importance_table,file = "improtance.csv")
  write.csv(trainh,file = "trainh.csv")
  write.csv(testh,file = "testh.csv")

  #############################################################################
  #############################################################################
  ### Step 1: Train Cox model using caret ----
  cox_ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    classProbs = FALSE,
    summaryFunction = defaultSummary
  )


  # cox_model <- coxph(cox_form, data = trainh)

  x <- as.matrix(trainh[, -c(1, 2)])  # Assuming the first two columns are time and status
  y <- Surv(as.numeric(trainh$time), as.numeric(trainh$status))


  tryCatch({
    set.seed(seed)
    # Modification 1: Explicitly set lambda.min.ratio and reduce nlambda in cv.glmnet
    cv_fit <- cv.glmnet(
      x = x,                    # Feature matrix
      y = y,                    # Surv object
      family = "cox",           # Cox model
      # type.measure = "mse",       # Using C-index
      nfolds = 5,               # 5-fold cross-validation
      alpha = 0.5,                # Lasso regularization (alpha=1)
      #lambda.min.ratio = 0.1,   # Increase the minimum lambda ratio to avoid numerical instability from very small values
      nlambda = 100,            # Reduce the number of lambdas to avoid generating too many small lambdas
      maxit = 1e5               # Increase the maximum number of iterations to ensure convergence
    )

    # View the best lambda value
    best_lambda <- cv_fit$lambda.min

    # Modification 2: Maintain parameter consistency in glmnet (no need to repeat lambda.min.ratio)
    cox_model <- glmnet(x = x,
                        y = y,
                        family = "cox",
                        lambda = best_lambda,     # Directly use the best lambda from cross-validation
                        alpha = 0.5,                # Keep alpha consistent with cv.glmnet
                        maxit = 1e5               # Same number of iterations as cv.glmnet

    )

    # View model coefficients
    # print(coef(cox_model))

    risk_score_trainh <-  predict(cox_model, newx = as.matrix(trainh[,-c(1,2)]), s = "lambda.min",type = "response")

    data_train_result = cbind(trainh[,c(1,2)],risk_score_trainh)

    colnames(data_train_result)[3] = "risk_score"

    # Predict risk scores on the test set
    risk_score_test <- predict(cox_model, newx = as.matrix(testh[,-c(1,2)]), s = "lambda.min",type = "response")


    data_test_result = cbind(testh[,c(1,2)],risk_score_test)
    colnames(data_test_result)[3] = "risk_score"

    # View risk scores
    # head(data_train_result)
    # head(data_test_result)
    data_train_result$status = as.numeric(data_train_result$status)
    data_test_result$status = as.numeric(data_test_result$status)

    data_train_result$status = data_train_result$status-1
    data_test_result$status = data_test_result$status-1
    ### Step 3: Determine the best cutoff ----
    # Need to use the original survival data to construct a Surv object
    # train_surv <- Surv(trainh$time, trainh$status)
    # Use survminer to determine the best cut-off point
    if (bestcut == T) {

      cut_point <- surv_cutpoint(
        data_train_result,
        time = "time",
        event = "status",
        variables = "risk_score"
      )
      # View the cutoff value
      # summary(cut_point)
      train_optimal_cutoff <- cut_point$cutpoint[1, "cutpoint"]
    } else{
      train_optimal_cutoff =  median(data_train_result$risk_score)
    }


    ### Step 4: Plot survival curves by group ----
    # Create grouping variable
    data_train_result <- data_train_result %>%
      mutate(risk_group = ifelse(risk_score > train_optimal_cutoff,
                                 "High Risk",
                                 "Low Risk"))

    # Create survival fit object
    fit <- survfit(Surv(time, status) ~ risk_group,data = data_train_result[,-3])

    require(patchwork)
    # Kaplan-Meier curve plotting
    p1= ggsurvplot(
      fit,
      data = data_train_result[,-3],
      pval = TRUE,            # Display log-rank p-value
      pval.method = TRUE,     # Display test method
      conf.int = TRUE,        # Display confidence intervals
      palette = c("#E41A1C", "#377EB8"),  # Custom colors
      title = "Survival Curve by Risk Group",
      xlab = "Time (days)",   # Time unit label
      #break.time.by = 100,    # Horizontal axis tick interval
      risk.table = TRUE,      # Display risk table
      risk.table.height = 0.25,  # Adjust risk table height
      ncensor.plot = FALSE,   # Do not display censoring event plot
      legend.labs = c("High Risk", "Low Risk")  # Legend labels
    )
    p1_1 = p1$plot+theme(plot.margin = margin(5, 5, 5, 5))
    p1_2 = p1$table+theme(plot.margin = margin(5, 5, 5, 5))

    # Set the proportion of the two plots, for example, set the height of p1 to twice that of p2
    combined_plot <- p1_1/p1_2 + plot_layout(heights = c(2.5, 1))
    # Export to a single plot
    ggsave("train_sur_plot.pdf", plot = combined_plot, width = 8, height = 10)


    ### Plot survival curves by group ----
    # Create grouping variable
    data_test_result <-data_test_result %>%
      mutate(risk_group = ifelse(risk_score > train_optimal_cutoff,
                                 "High Risk",
                                 "Low Risk"))

    # Create survival fit object
    surv_fit <- survfit(
      Surv(time, status) ~ risk_group,
      data = data_test_result
    )
    if (length(unique(data_test_result$risk_group))>1) {
      # Kaplan-Meier curve plotting
      p2 = ggsurvplot(
        surv_fit,
        data = data_test_result,
        pval = TRUE,            # Display log-rank p-value
        pval.method = TRUE,     # Display test method
        conf.int = TRUE,        # Display confidence intervals
        palette = c("#E41A1C", "#377EB8"),  # Custom colors
        title = "Survival Curve by Risk Group",
        xlab = "Time (days)",   # Time unit label
        break.time.by = 100,    # Horizontal axis tick interval
        risk.table = TRUE,      # Display risk table
        risk.table.height = 0.25,  # Adjust risk table height
        ncensor.plot = FALSE,   # Do not display censoring event plot
        legend.labs = c("High Risk", "Low Risk")  # Legend labels
      )

      p2_1 = p2$plot+theme(plot.margin = margin(5, 5, 5, 5))
      p2_2 = p2$table+theme(plot.margin = margin(5, 5, 5, 5))

      # Set the proportion of the two plots, for example, set the height of p1 to twice that of p2
      combined_plot <- p2_1/p2_2 + plot_layout(heights = c(2.5, 1))
      # Export to a single plot
      ggsave("test_sur_plot.pdf", plot = combined_plot, width = 8, height = 10)
    }

  },error = function(e) {
    # Print error message
    message("An error occurred: ", e$message)
    # Additional error handling logic can be added here
  })
  return(list(data_train_result,data_test_result))
}

#' Machine Learning Pipeline for Cox Proportional Hazards Model (Alternative Version)
#'
#' An alternative implementation of the survival analysis pipeline using Cox proportional hazards model with Elastic Net regularization.
#' This version performs feature selection before data splitting, which may be more suitable for certain datasets.
#'
#' @param data A dataframe where the first two columns are survival time and status.
#' @param seed Random seed for reproducibility (default: 123)
#' @param train_percent Proportion of data to use for training (default: 0.7)
#' @param n_round Number of bootstrap rounds for feature selection (default: 20)
#' @param stratify Logical indicating whether to use stratified sampling (default: FALSE)
#' @param lambda_type Regularization criterion for cv.glmnet: "lambda.min" (default) or "lambda.1se"
#' @param bestcut Logical indicating whether to determine optimal cutoff for risk groups (default: FALSE)
#' @param topn Number of top features to retain in feature selection (default: 10)
#'
#' @return A list containing:
#' \itemize{
#'   \item Training set results with risk scores and groups
#'   \item Testing set results with risk scores and groups
#' }
#'
#' @details
#' Key differences from machine_learning_cox:
#' 1. Performs feature selection on the entire dataset before splitting
#' 2. Scales data before feature selection
#' 3. May be more appropriate for smaller datasets
#'
#' The function performs the following steps:
#' 1. Data scaling and feature selection
#' 2. Data splitting
#' 3. Model training with Elastic Net regularization
#' 4. Risk score calculation and optimal cutoff determination
#' 5. Survival curve visualization and PDF export
#'
#' @examples
#' \dontrun{
#' data <- data.frame(time = rnorm(100), status = rbinom(100, 1, 0.5),
#'                   matrix(rnorm(100*20), ncol=20))
#' results <- machine_learning_cox_2(data, train_percent = 0.8, topn = 5)
#' }
#'
#' @export
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival Surv survfit
#' @importFrom survminer surv_cutpoint ggsurvplot
#' @importFrom caret trainControl
#' @importFrom patchwork plot_layout
#' @importFrom dplyr mutate
machine_learning_cox_2 = function(data,
                                  seed = 123,
                                  train_percent = 0.7,
                                  n_round = 20,
                                  stratify = F,
                                  lambda_type = "lambda.min",
                                  bestcut = F,
                                  topn = 10) {

  colnames(data)[c(1, 2)] = c("time", "status")

  data[, -c(1, 2)] = scale(data[, -c(1, 2)])

  x = feature_selection_coxlasso(data = data, seed = seed,
                                 n_round = n_round, topn = topn, lambda_type = lambda_type)

  data = x$selected_matrix

  # trainh[,-c(1,2)] = scale(trainh[,-c(1,2)])
  # testh[,-c(1,2)] = scale(testh[,-c(1,2)])

  # When the sample size is small (n < 100), the minimum sample size can be increased:
  xx = surv_data_split(data, train_ratio = train_percent,
                       stratify = stratify,
                       seed = seed)

  trainh <- xx$train
  testh <- xx$test

  trainh = trainh[, apply(trainh, 2, var) != 0]
  testh = testh[, apply(testh, 2, var) != 0]

  testh = testh[, colnames(testh) %in% colnames(trainh)]
  trainh = trainh[, colnames(trainh) %in% colnames(testh)]

  write.csv(x$coef_matrix, "coef_matrix.csv")
  write.csv(x$importance_table, file = "improtance.csv")
  write.csv(trainh, file = "trainh.csv")
  write.csv(testh, file = "testh.csv")

  #############################################################################
  #############################################################################
  ### Step 1: Train Cox model using caret ----
  cox_ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    classProbs = FALSE,
    summaryFunction = defaultSummary
  )


  # cox_model <- coxph(cox_form, data = trainh)

  x <- as.matrix(trainh[, -c(1, 2)])  # Assuming the first two columns are time and status
  y <- Surv(as.numeric(trainh$time), as.numeric(trainh$status))


  tryCatch({
    set.seed(seed)
    # Modification 1: Explicitly set lambda.min.ratio and reduce nlambda in cv.glmnet
    cv_fit <- cv.glmnet(
      x = x,                    # Feature matrix
      y = y,                    # Surv object
      family = "cox",           # Cox model
      # type.measure = "mse",       # Using C-index
      nfolds = 5,               # 5-fold cross-validation
      alpha = 0.5,              # Lasso regularization (alpha=1)
      # lambda.min.ratio = 0.1,   # Increase the minimum lambda ratio to avoid numerical instability from very small values
      nlambda = 100,            # Reduce the number of lambdas to avoid generating too many small lambdas
      maxit = 1e5               # Increase the maximum number of iterations to ensure convergence
    )

    # View the best lambda value
    best_lambda <- cv_fit$lambda.min

    # Modification 2: Maintain parameter consistency in glmnet (no need to specify lambda.min.ratio again)
    cox_model <- glmnet(x = x,
                        y = y,
                        family = "cox",
                        lambda = best_lambda,     # Use the best lambda obtained from cross-validation directly
                        alpha = 0.5,              # Keep the same alpha as in cv.glmnet
                        maxit = 1e5               # Same number of iterations as cv.glmnet

    )

    # View model coefficients
    # print(coef(cox_model))

    risk_score_trainh <- predict(cox_model, newx = as.matrix(trainh[, -c(1, 2)]), s = "lambda.min", type = "response")

    data_train_result = cbind(trainh[, c(1, 2)], risk_score_trainh)

    colnames(data_train_result)[3] = "risk_score"

    # Predict risk scores on the test set
    risk_score_test <- predict(cox_model, newx = as.matrix(testh[, -c(1, 2)]), s = "lambda.min", type = "response")


    data_test_result = cbind(testh[, c(1, 2)], risk_score_test)
    colnames(data_test_result)[3] = "risk_score"

    # View risk scores
    # head(data_train_result)
    # head(data_test_result)
    data_train_result$status = as.numeric(data_train_result$status)
    data_test_result$status = as.numeric(data_test_result$status)

    data_train_result$status = data_train_result$status - 1
    data_test_result$status = data_test_result$status - 1
    ### Step 3: Determine the best cutoff ----
    # Need to use the original survival data to construct the Surv object
    # train_surv <- Surv(trainh$time, trainh$status)
    # Use survminer to determine the best cut-off point
    if (bestcut == T) {

      cut_point <- surv_cutpoint(
        data_train_result,
        time = "time",
        event = "status",
        variables = "risk_score"
      )
      # View the cutoff value
      # summary(cut_point)
      train_optimal_cutoff <- cut_point$cutpoint[1, "cutpoint"]
    } else {
      train_optimal_cutoff = median(data_train_result$risk_score)
    }


    ### Step 4: Group and plot survival curves ----
    # Create grouping variable
    data_train_result <- data_train_result %>%
      mutate(risk_group = ifelse(risk_score > train_optimal_cutoff,
                                 "High Risk",
                                 "Low Risk"))

    # Create survival fit object
    fit <- survfit(Surv(time, status) ~ risk_group, data = data_train_result[, -3])

    require(patchwork)
    # Kaplan-Meier curve plotting
    p1 = ggsurvplot(
      fit,
      data = data_train_result[, -3],
      pval = TRUE,            # Display log-rank p-value
      pval.method = TRUE,     # Display test method
      conf.int = TRUE,        # Display confidence intervals
      palette = c("#E41A1C", "#377EB8"),  # Custom colors
      title = "Survival Curve by Risk Group",
      xlab = "Time (days)",   # Time unit label
      # break.time.by = 100,    # Horizontal axis tick interval
      risk.table = TRUE,      # Display risk table
      risk.table.height = 0.25,  # Adjust risk table height
      ncensor.plot = FALSE,   # Do not display censoring event plot
      legend.labs = c("High Risk", "Low Risk")  # Legend labels
    )
    p1_1 = p1$plot + theme(plot.margin = margin(5, 5, 5, 5))
    p1_2 = p1$table + theme(plot.margin = margin(5, 5, 5, 5))

    # Set the proportion of the two plots, for example, set the height of p1 to be twice that of p2
    combined_plot <- p1_1/p1_2 + plot_layout(heights = c(2.5, 1))
    # Export to a single plot
    ggsave("train_sur_plot.pdf", plot = combined_plot, width = 8, height = 10)


    ### Group and plot survival curves ----
    # Create grouping variable
    data_test_result <- data_test_result %>%
      mutate(risk_group = ifelse(risk_score > train_optimal_cutoff,
                                 "High Risk",
                                 "Low Risk"))

    # Create survival fit object
    surv_fit <- survfit(
      Surv(time, status) ~ risk_group,
      data = data_test_result
    )
    if (length(unique(data_test_result$risk_group)) > 1) {
      # Kaplan-Meier curve plotting
      p2 = ggsurvplot(
        surv_fit,
        data = data_test_result,
        pval = TRUE,            # Display log-rank p-value
        pval.method = TRUE,     # Display test method
        conf.int = TRUE,        # Display confidence intervals
        palette = c("#E41A1C", "#377EB8"),  # Custom colors
        title = "Survival Curve by Risk Group",
        xlab = "Time (days)",   # Time unit label
        break.time.by = 100,    # Horizontal axis tick interval
        risk.table = TRUE,      # Display risk table
        risk.table.height = 0.25,  # Adjust risk table height
        ncensor.plot = FALSE,   # Do not display censoring event plot
        legend.labs = c("High Risk", "Low Risk")  # Legend labels
      )

      p2_1 = p2$plot + theme(plot.margin = margin(5, 5, 5, 5))
      p2_2 = p2$table + theme(plot.margin = margin(5, 5, 5, 5))

      # Set the proportion of the two plots, for example, set the height of p1 to be twice that of p2
      combined_plot <- p2_1/p2_2 + plot_layout(heights = c(2.5, 1))
      # Export to a single plot
      ggsave("test_sur_plot.pdf", plot = combined_plot, width = 8, height = 10)
    }

  }, error = function(e) {
    # Print error message
    message("An error occurred: ", e$message)
    # Additional error handling logic can be added here
  })
  return(list(data_train_result, data_test_result))
}


#' Feature Selection using Random Forest with Stratified Bootstrap Resampling
#'
#' Performs feature selection using Random Forest with stratified bootstrap resampling.
#' Aggregates feature importance across multiple resampling rounds to identify top features.
#'
#' @param data A dataframe or matrix where the first column contains group labels.
#' @param seed Random seed for reproducibility (default: 123)
#' @param n_round Number of bootstrap resampling rounds (default: 20)
#' @param topn Number of top features to retain (default: 10)
#' @param importance_type Type of importance measure: 1 for mean decrease in accuracy,
#'                        2 for mean decrease in Gini index (default: 2)
#'
#' @return A list containing:
#' \itemize{
#'   \item selected_matrix - Data with selected features and group labels
#'   \item importance_table - Dataframe of feature rankings by average importance
#' }
#'
#' @details
#' Implements stratified bootstrap resampling preserving class proportions. Uses Random Forest
#' with 500 trees for stability. Feature importance is averaged across resampling rounds and
#' normalized for comparison.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(group = sample(0:1, 100, replace = TRUE),
#'                   matrix(rnorm(100*20), ncol=20))
#' result <- feature_selection_rf(data, topn = 5)
#' }
#'
#' @export
#' @importFrom randomForest randomForest importance
## RF feature selection
feature_selection_rf <- function(data, seed = 123, n_round = 20, topn = 10, importance_type = 2) {
  # Parameter description:
  # data: Data frame or matrix, with the first column as the grouping variable
  # seed: Random seed to control reproducibility
  # n_round: Number of resampling rounds
  # topn: Number of top features to retain
  # importance_type: Importance type (1 = accuracy decrease, 2 = Gini decrease)

  set.seed(seed)
  group_col <- 1
  X <- data[, -group_col, drop = FALSE]
  y <- as.factor(data[, group_col])

  importance_matrix <- matrix(0, nrow = n_round, ncol = ncol(X))
  colnames(importance_matrix) <- colnames(X)

  for (i in 1:n_round) {
    # Stratified Bootstrap sampling (maintain class proportions)
    indices_by_class <- split(1:nrow(X), y)
    sampled_indices <- unlist(lapply(indices_by_class, function(ind) {
      sample(ind, size = length(ind), replace = TRUE)
    }))

    X_sample <- X[sampled_indices, , drop = FALSE]
    y_sample <- y[sampled_indices]

    # Train model (random forest hyperparameters can be adjusted)
    rf_model <- randomForest::randomForest(
      x = X_sample,
      y = y_sample,
      importance = TRUE,
      ntree = 500  # Increase the number of trees to improve stability
    )

    # Record importance and normalize
    imp <- randomForest::importance(rf_model, type = importance_type)[, 1]
    importance_matrix[i, ] <- imp / sum(imp)  # Normalization
  }

  avg_importance <- colMeans(importance_matrix)
  ranked_features <- names(sort(avg_importance, decreasing = TRUE))[1:topn]

  # Construct the result (retain original data, only filter features)
  selected_data <- cbind(data[, group_col, drop = FALSE], X[, ranked_features, drop = FALSE])
  colnames(selected_data)[1] <- "group"

  importance_table <- data.frame(
    Feature = names(avg_importance),
    Importance = avg_importance,
    Rank = rank(-avg_importance, ties.method = "min")
  )[order(-avg_importance), ]
  rownames(importance_table) <- NULL

  return(list(
    selected_matrix = selected_data,
    importance_table = importance_table
  ))
}


#' Machine Learning Pipeline for Classification Tasks
#'
#' A comprehensive pipeline for training and evaluating multiple machine learning models for classification tasks.
#' Includes data splitting, feature selection, model training, evaluation, and visualization.
#'
#' @param data A dataframe where the first column contains group labels.
#' @param seed Random seed for reproducibility (default: 123)
#' @param group_name Column index or name of the group labels (default: 1)
#' @param train_percent Proportion of data to use for training (default: 0.7)
#' @param n_round Number of bootstrap rounds for feature selection (default: 20)
#' @param topn Number of top features to retain in feature selection (default: 10)
#'
#' @return A dataframe containing AUC values for each model (glmnet, svm, rf, nnet).
#'
#' @details
#' The function performs the following steps:
#' 1. Data preprocessing and splitting into training and testing sets.
#' 2. Feature selection using Random Forest.
#' 3. Training and evaluation of four models: glmnet (logistic regression), svm (support vector machine),
#'    rf (random forest), and nnet (neural network).
#' 4. Generation of ROC curves, confusion matrices, and model importance tables.
#' 5. Export of results, including model parameters, predictions, and visualizations.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(group = sample(0:1, 100, replace = TRUE),
#'                   matrix(rnorm(100*20), ncol=20))
#' results <- machine_learning(data, topn = 5)
#' }
#'
#' @export
#' @importFrom caret train trainControl createDataPartition confusionMatrix varImp
#' @importFrom pROC plot.roc
#' @importFrom pheatmap pheatmap
#' @importFrom xlsx write.xlsx

machine_learning = function(data, seed = 123, group_name = 1, train_percent = 0.7, n_round = 20,
                            topn = 10) {
  colnames(data)[group_name] = "group"
  data$group = factor(data$group)
  ## Split the data
  set.seed(seed)
  indexh <- createDataPartition(data[, "group"], p = train_percent, list = F)
  trainh <- data[indexh, ]
  testh <- data[-indexh, ]

  # x = rfe_topn(data = trainh,topn = topn,model_type = model_type,seed = seed)
  x = feature_selection_rf(data = trainh, seed = seed, n_round = n_round, topn = topn)

  trainh = x$selected_matrix
  testh = testh[, colnames(testh) %in% colnames(trainh)]

  write.csv(x$importance_table, file = "improtance.csv")
  write.csv(trainh, file = "trainh.csv")
  write.csv(testh, file = "testh.csv")

  data_auc = data.frame(auc = c(rep(0, 4)), row.names = c("glmnet", "svm", "rf", "nnet"))
  set.seed(seed)
  m1h <- train(group ~ .,  # The dependent variable is group, and the independent variables are all other indicators except group, represented by group ~ .
               data = trainh,  # Train using the training set
               method = "glmnet",  # glm represents generalized linear model; glmnet is regularized logistic regression to prevent overfitting
               preProcess = c("scale", "center"),  # Built-in data standardization parameters
               tuneGrid = data.frame(alpha = c(0, 0.001, 0.01, 0.05, 0.1, 0.15), lambda = c(0, 0.001, 0.01, 0.05, 0.1, 0.15)),
               # trControl is a built-in model training control parameter in caret, containing cross-validation and output type parameters, which vary for different models
               trControl = trainControl(classProbs = T,
                                        method = "repeatedcv",
                                        number = 10,
                                        verboseIter = T,
                                        repeats = 5
               ))
  glm_best_tune = m1h$bestTune

  best_lambda <- m1h$bestTune$lambda
  best_alpha <- m1h$bestTune$alpha

  # Get the coefficients of the best model
  best_coefficients <- as.data.frame(as.matrix(coef(m1h$finalModel, s = best_lambda)))
  write.csv(best_coefficients, "glmnet_model_coef.csv")
  # Prediction
  testh_pred <- predict(m1h, testh, type = "prob")
  write.csv(testh_pred, "glmnet_predict.csv")
  # Plot ROC
  pdf(file = "1、glm_model_roc.pdf", width = 8, height = 8)
  roc1 <- plot.roc(as.factor(testh$group), testh_pred[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = T,
                   main = "Area under the curve for Logistic Regression(glmnet)", print.thres = "best")
  dev.off()
  glm_roc = round(as.numeric(sub(".*:", "", roc1$auc)) / 100, digits = 3)

  data_auc["glmnet", ] = glm_roc

  sumroc = roc1$sensitivities + roc1$specificities
  cutvalue = roc1$thresholds[which(sumroc == max(sumroc))]
  glm_pred = ifelse(testh_pred[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  glm_conf_matr = as.matrix(caret::confusionMatrix(factor(glm_pred), testh$group))

  pdf(file = "1、glm_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(glm_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T,
  )
  dev.off()


  #
  # Model training based on k-fold cross-validation, linear SVM
  set.seed(seed)
  svm1 <- train(group ~ .,
                data = trainh,
                method = "svmLinear",
                metric = "Accuracy",
                tuneGrid = data.frame(C = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1)),
                preProcess = c("scale", "center"),
                trControl = trainControl(classProbs = T,
                                         method = "repeatedcv",
                                         number = 10,
                                         verboseIter = T,
                                         repeats = 5,
                ))

  svm_best_tune = svm1$bestTune

  best_C_svm <- svm1$bestTune$C

  # Get the model's feature importance evaluation
  importance <- data.frame(row.names = row.names(varImp(svm1)$importance), improtance = varImp(svm1)$importance[, 1])

  write.csv(importance, "svm_ importance.csv")

  # Prediction
  pred_svmh <- predict(svm1, testh, type = "prob")
  write.csv(pred_svmh, "svm_predict.csv")
  # Plot ROC
  pdf(file = "2、svm_model_roc.pdf", width = 8, height = 8)
  roc2 <- plot.roc(as.factor(testh$group), pred_svmh[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = TRUE,
                   main = "Area under the curve for SVM", print.thres = "best")
  dev.off()

  svm_roc = round(as.numeric(sub(".*:", "", roc2$auc)) / 100, digits = 3)

  data_auc["svm", ] = svm_roc

  sumroc2 = roc2$sensitivities + roc2$specificities
  cutvalue = roc2$thresholds[which(sumroc2 == max(sumroc2))]
  svm_pred = ifelse(testh_pred[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  svm_conf_matr = as.matrix(caret::confusionMatrix(factor(svm_pred), testh$group))

  pdf(file = "2、svm_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(svm_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T,
  )
  dev.off()

  ###### RF
  set.seed(seed)

  if (topn >= 20) {
    mtryx <- c(1, 2, 3, 4, 5, 6, 10, 20)
  } else if (topn > 5) {
    mtryx <- c(1, 2, 3, 4, 5, 6, 10)
  } else {
    mtryx <- c(1, 2, 3)
  }

  model_rfh <- train(group ~ .,
                     data = trainh,
                     method = "rf",
                     tuneGrid = data.frame(mtry = mtryx),
                     preProcess = c("scale", "center"),
                     trControl = trainControl(method = "repeatedcv",
                                              number = 10,
                                              repeats = 5,
                                              savePredictions = TRUE,
                                              verboseIter = T,
                                              classProbs = T))
  rf_best_tune = model_rfh$bestTune

  # Get the model's feature importance evaluation
  importance <- data.frame(row.names = row.names(varImp(model_rfh)$importance), improtance = varImp(model_rfh)$importance[, 1])

  write.csv(importance, "rf_ importance.csv")

  # Prediction
  pred_rfh <- predict(model_rfh, testh, type = "prob")
  write.csv(pred_rfh, "rf_predict.csv")
  # ROC curve
  pdf(file = "3、rf_model_roc.pdf", width = 8, height = 8)
  roc3 <- plot.roc(testh$group, pred_rfh[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = TRUE,
                   main = "Area under the curve for Random Forest", print.thres = "best")
  dev.off()

  rf_roc = round(as.numeric(sub(".*:", "", roc3$auc)) / 100, digits = 3)

  data_auc["rf", ] = rf_roc
  # ####
  # rf_roc = round(as.numeric(sub(".*:", "", roc3$auc))/100,digits = 3)

  sumroc3 = roc3$sensitivities + roc3$specificities
  cutvalue = roc3$thresholds[which(sumroc3 == max(sumroc3))]
  rf_pred = ifelse(pred_rfh[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  rf_conf_matr = as.matrix(caret::confusionMatrix(factor(rf_pred), testh$group))

  pdf(file = "3、rf_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(rf_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T)
  dev.off()

  #### nnet
  set.seed(seed)
  model_nneth <- caret::train(group ~ .,
                              data = trainh,
                              method = "nnet",
                              tuneGrid = data.frame(size = c(1, 3, 5, 10, 15, 20, 50), decay = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 5)),
                              preProcess = c("scale", "center"),
                              trControl = trainControl(method = "repeatedcv",
                                                       number = 10,
                                                       classProbs = T,
                                                       repeats = 5,
                                                       savePredictions = TRUE,
                                                       verboseIter = FALSE))

  nnet_bset_tune = model_nneth$bestTune

  importance <- data.frame(row.names = row.names(varImp(model_nneth)$importance),
                           improtance = varImp(model_nneth)$importance[, 1])

  write.csv(importance, "nnet_importance.csv")

  pred_ann <- predict(model_nneth, testh, type = "prob")
  write.csv(pred_ann, "nnet_predict.csv")
  # Plot ROC
  pdf(file = "4、nnet_model_roc.pdf", width = 8, height = 8)
  roc4 <- plot.roc(testh$group, pred_ann[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = TRUE,
                   print.thres = T, main = "Area under the curve for neural networks")
  dev.off()
  nnet_roc = round(as.numeric(sub(".*:", "", roc4$auc)) / 100, digits = 3)

  data_auc["nnet", ] = nnet_roc

  # ####
  # nnet_roc = round(as.numeric(sub(".*:", "", roc4$auc))/100,digits = 3)

  sumroc4 = roc4$sensitivities + roc4$specificities
  cutvalue = roc4$thresholds[which(sumroc4 == max(sumroc4))]
  nnet_pred = ifelse(pred_ann[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  nnet_conf_matr = as.matrix(caret::confusionMatrix(factor(nnet_pred), testh$group))

  pdf(file = "4、nnet_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(nnet_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T,
  )
  dev.off()

  # Plot combined ROC
  pdf(file = "5、model_roc_combine.pdf", width = 8, height = 8)
  roc1 <- plot.roc(as.factor(testh$group), testh_pred[, unique(testh$group)[1]], percent = TRUE, col = "#E41A1C", print.auc = F,
                   main = "Area under the curve for different models")
  plot(roc2, col = "#377EB8", add = T)
  plot(roc3, col = "#4DAF4A", add = T)
  plot(roc4, col = "#984EA3", add = T)

  text.legend <- c(paste0("Logistic AUC:", glm_roc), paste0("SVM AUC:", svm_roc),
                   paste0("Random Forest AUC:", rf_roc), paste0("Nnet AUC:", nnet_roc))
  col.legend <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  legend("bottomright", pch = c(20, 20, 20, 20), legend = text.legend, cex = 1.2,
         col = col.legend, bty = "n", horiz = F)
  dev.off()

  ### Export the best model parameters
  xlsx::write.xlsx(glm_best_tune, file = "param.xlsx", sheetName = "glmnet", row.names = F)
  xlsx::write.xlsx(svm_best_tune, file = "param.xlsx", sheetName = "svm", row.names = F, append = T)
  xlsx::write.xlsx(rf_best_tune, file = "param.xlsx", sheetName = "rf", row.names = F, append = T)
  xlsx::write.xlsx(nnet_bset_tune, file = "param.xlsx", sheetName = "nnet", row.names = F, append = T)

  return(data_auc)
}





#' Machine Learning Pipeline for Classification Tasks
#'
#' A comprehensive pipeline for training and evaluating multiple machine learning models for classification tasks.
#' Includes data splitting, feature selection, model training, evaluation, and visualization.
#'
#' @param data A dataframe where the first column contains group labels.
#' @param seed Random seed for reproducibility (default: 123)
#' @param group_name Column index or name of the group labels (default: 1)
#' @param train_percent Proportion of data to use for training (default: 0.7)
#' @param n_round Number of bootstrap rounds for feature selection (default: 20)
#' @param topn Number of top features to retain in feature selection (default: 10)
#'
#' @return A dataframe containing AUC values for each model (glmnet, svm, rf, nnet).
#'
#' @details
#' The function performs the following steps:
#' 1. Data preprocessing and splitting into training and testing sets.
#' 2. Feature selection using Random Forest.
#' 3. Training and evaluation of four models: glmnet (logistic regression), svm (support vector machine),
#'    rf (random forest), and nnet (neural network).
#' 4. Generation of ROC curves, confusion matrices, and model importance tables.
#' 5. Export of results, including model parameters, predictions, and visualizations.
#'
#' @examples
#' \dontrun{
#' data <- data.frame(group = sample(0:1, 100, replace = TRUE),
#'                   matrix(rnorm(100*20), ncol=20))
#' results <- machine_learning(data, topn = 5)
#' }
#'
#' @export
#' @importFrom caret train trainControl createDataPartition confusionMatrix varImp
#' @importFrom pROC plot.roc
#' @importFrom pheatmap pheatmap
#' @importFrom xlsx write.xlsx

machine_learning2 = function(data, seed = 123, group_name = 1, train_percent = 0.7, n_round = 20,
                            topn = 10) {
  colnames(data)[group_name] = "group"
  data$group = factor(data$group)

  # x = rfe_topn(data = trainh,topn = topn,model_type = model_type,seed = seed)
  x = feature_selection_rf(data = data, seed = seed, n_round = n_round, topn = topn)

  data = x$selected_matrix

  ## Split the data
  set.seed(seed)
  indexh <- createDataPartition(data[, "group"], p = train_percent, list = F)
  trainh <- data[indexh, ]
  testh <- data[-indexh, ]

  write.csv(x$importance_table, file = "improtance.csv")
  write.csv(trainh, file = "trainh.csv")
  write.csv(testh, file = "testh.csv")

  data_auc = data.frame(auc = c(rep(0, 4)), row.names = c("glmnet", "svm", "rf", "nnet"))
  set.seed(seed)
  m1h <- train(group ~ .,  # The dependent variable is group, and the independent variables are all other indicators except group, represented by group ~ .
               data = trainh,  # Train using the training set
               method = "glmnet",  # glm represents generalized linear model; glmnet is regularized logistic regression to prevent overfitting
               preProcess = c("scale", "center"),  # Built-in data standardization parameters
               tuneGrid = data.frame(alpha = c(0, 0.001, 0.01, 0.05, 0.1, 0.15), lambda = c(0, 0.001, 0.01, 0.05, 0.1, 0.15)),
               # trControl is a built-in model training control parameter in caret, containing cross-validation and output type parameters, which vary for different models
               trControl = trainControl(classProbs = T,
                                        method = "repeatedcv",
                                        number = 10,
                                        verboseIter = T,
                                        repeats = 5
               ))
  glm_best_tune = m1h$bestTune

  best_lambda <- m1h$bestTune$lambda
  best_alpha <- m1h$bestTune$alpha

  # Get the coefficients of the best model
  best_coefficients <- as.data.frame(as.matrix(coef(m1h$finalModel, s = best_lambda)))
  write.csv(best_coefficients, "glmnet_model_coef.csv")
  # Prediction
  testh_pred <- predict(m1h, testh, type = "prob")
  write.csv(testh_pred, "glmnet_predict.csv")
  # Plot ROC
  pdf(file = "1、glm_model_roc.pdf", width = 8, height = 8)
  roc1 <- plot.roc(as.factor(testh$group), testh_pred[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = T,
                   main = "Area under the curve for Logistic Regression(glmnet)", print.thres = "best")
  dev.off()
  glm_roc = round(as.numeric(sub(".*:", "", roc1$auc)) / 100, digits = 3)

  data_auc["glmnet", ] = glm_roc

  sumroc = roc1$sensitivities + roc1$specificities
  cutvalue = roc1$thresholds[which(sumroc == max(sumroc))]
  glm_pred = ifelse(testh_pred[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  glm_conf_matr = as.matrix(caret::confusionMatrix(factor(glm_pred), testh$group))

  pdf(file = "1、glm_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(glm_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T,
  )
  dev.off()


  #
  # Model training based on k-fold cross-validation, linear SVM
  set.seed(seed)
  svm1 <- train(group ~ .,
                data = trainh,
                method = "svmLinear",
                metric = "Accuracy",
                tuneGrid = data.frame(C = c(0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 1)),
                preProcess = c("scale", "center"),
                trControl = trainControl(classProbs = T,
                                         method = "repeatedcv",
                                         number = 10,
                                         verboseIter = T,
                                         repeats = 5,
                ))

  svm_best_tune = svm1$bestTune

  best_C_svm <- svm1$bestTune$C

  # Get the model's feature importance evaluation
  importance <- data.frame(row.names = row.names(varImp(svm1)$importance), improtance = varImp(svm1)$importance[, 1])

  write.csv(importance, "svm_ importance.csv")

  # Prediction
  pred_svmh <- predict(svm1, testh, type = "prob")
  write.csv(pred_svmh, "svm_predict.csv")
  # Plot ROC
  pdf(file = "2、svm_model_roc.pdf", width = 8, height = 8)
  roc2 <- plot.roc(as.factor(testh$group), pred_svmh[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = TRUE,
                   main = "Area under the curve for SVM", print.thres = "best")
  dev.off()

  svm_roc = round(as.numeric(sub(".*:", "", roc2$auc)) / 100, digits = 3)

  data_auc["svm", ] = svm_roc

  sumroc2 = roc2$sensitivities + roc2$specificities
  cutvalue = roc2$thresholds[which(sumroc2 == max(sumroc2))]
  svm_pred = ifelse(testh_pred[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  svm_conf_matr = as.matrix(caret::confusionMatrix(factor(svm_pred), testh$group))

  pdf(file = "2、svm_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(svm_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T,
  )
  dev.off()

  ###### RF
  set.seed(seed)

  if (topn >= 20) {
    mtryx <- c(1, 2, 3, 4, 5, 6, 10, 20)
  } else if (topn > 5) {
    mtryx <- c(1, 2, 3, 4, 5, 6, 10)
  } else {
    mtryx <- c(1, 2, 3)
  }

  model_rfh <- train(group ~ .,
                     data = trainh,
                     method = "rf",
                     tuneGrid = data.frame(mtry = mtryx),
                     preProcess = c("scale", "center"),
                     trControl = trainControl(method = "repeatedcv",
                                              number = 10,
                                              repeats = 5,
                                              savePredictions = TRUE,
                                              verboseIter = T,
                                              classProbs = T))
  rf_best_tune = model_rfh$bestTune

  # Get the model's feature importance evaluation
  importance <- data.frame(row.names = row.names(varImp(model_rfh)$importance), improtance = varImp(model_rfh)$importance[, 1])

  write.csv(importance, "rf_ importance.csv")

  # Prediction
  pred_rfh <- predict(model_rfh, testh, type = "prob")
  write.csv(pred_rfh, "rf_predict.csv")
  # ROC curve
  pdf(file = "3、rf_model_roc.pdf", width = 8, height = 8)
  roc3 <- plot.roc(testh$group, pred_rfh[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = TRUE,
                   main = "Area under the curve for Random Forest", print.thres = "best")
  dev.off()

  rf_roc = round(as.numeric(sub(".*:", "", roc3$auc)) / 100, digits = 3)

  data_auc["rf", ] = rf_roc
  # ####
  # rf_roc = round(as.numeric(sub(".*:", "", roc3$auc))/100,digits = 3)

  sumroc3 = roc3$sensitivities + roc3$specificities
  cutvalue = roc3$thresholds[which(sumroc3 == max(sumroc3))]
  rf_pred = ifelse(pred_rfh[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  rf_conf_matr = as.matrix(caret::confusionMatrix(factor(rf_pred), testh$group))

  pdf(file = "3、rf_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(rf_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T)
  dev.off()

  #### nnet
  set.seed(seed)
  model_nneth <- caret::train(group ~ .,
                              data = trainh,
                              method = "nnet",
                              tuneGrid = data.frame(size = c(1, 3, 5, 10, 15, 20, 50), decay = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 5)),
                              preProcess = c("scale", "center"),
                              trControl = trainControl(method = "repeatedcv",
                                                       number = 10,
                                                       classProbs = T,
                                                       repeats = 5,
                                                       savePredictions = TRUE,
                                                       verboseIter = FALSE))

  nnet_bset_tune = model_nneth$bestTune

  importance <- data.frame(row.names = row.names(varImp(model_nneth)$importance),
                           improtance = varImp(model_nneth)$importance[, 1])

  write.csv(importance, "nnet_importance.csv")

  pred_ann <- predict(model_nneth, testh, type = "prob")
  write.csv(pred_ann, "nnet_predict.csv")
  # Plot ROC
  pdf(file = "4、nnet_model_roc.pdf", width = 8, height = 8)
  roc4 <- plot.roc(testh$group, pred_ann[, unique(testh$group)[1]], percent = TRUE, col = "#1c61b6", print.auc = TRUE,
                   print.thres = T, main = "Area under the curve for neural networks")
  dev.off()
  nnet_roc = round(as.numeric(sub(".*:", "", roc4$auc)) / 100, digits = 3)

  data_auc["nnet", ] = nnet_roc

  # ####
  # nnet_roc = round(as.numeric(sub(".*:", "", roc4$auc))/100,digits = 3)

  sumroc4 = roc4$sensitivities + roc4$specificities
  cutvalue = roc4$thresholds[which(sumroc4 == max(sumroc4))]
  nnet_pred = ifelse(pred_ann[, unique(testh$group)[1]] > cutvalue, as.character(unique(testh$group)[1]), as.character(unique(testh$group)[2]))
  nnet_conf_matr = as.matrix(caret::confusionMatrix(factor(nnet_pred), testh$group))

  pdf(file = "4、nnet_conf_heatmap.pdf", width = 7, height = 7)
  pheatmap::pheatmap(nnet_conf_matr, display_numbers = T, fontsize = 15,
                     cluster_rows = F, cluster_cols = F,
                     show_colnames = T, show_rownames = T,
  )
  dev.off()

  # Plot combined ROC
  pdf(file = "5、model_roc_combine.pdf", width = 8, height = 8)
  roc1 <- plot.roc(as.factor(testh$group), testh_pred[, unique(testh$group)[1]], percent = TRUE, col = "#E41A1C", print.auc = F,
                   main = "Area under the curve for different models")
  plot(roc2, col = "#377EB8", add = T)
  plot(roc3, col = "#4DAF4A", add = T)
  plot(roc4, col = "#984EA3", add = T)

  text.legend <- c(paste0("Logistic AUC:", glm_roc), paste0("SVM AUC:", svm_roc),
                   paste0("Random Forest AUC:", rf_roc), paste0("Nnet AUC:", nnet_roc))
  col.legend <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")
  legend("bottomright", pch = c(20, 20, 20, 20), legend = text.legend, cex = 1.2,
         col = col.legend, bty = "n", horiz = F)
  dev.off()

  ### Export the best model parameters
  xlsx::write.xlsx(glm_best_tune, file = "param.xlsx", sheetName = "glmnet", row.names = F)
  xlsx::write.xlsx(svm_best_tune, file = "param.xlsx", sheetName = "svm", row.names = F, append = T)
  xlsx::write.xlsx(rf_best_tune, file = "param.xlsx", sheetName = "rf", row.names = F, append = T)
  xlsx::write.xlsx(nnet_bset_tune, file = "param.xlsx", sheetName = "nnet", row.names = F, append = T)

  return(data_auc)
}






