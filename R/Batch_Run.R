#' Batch Execution of Machine Learning Workflow with Multiple Seeds
#'
#' @description
#' Repeatedly executes machine learning workflows with different random seeds, stores results in seed-specific folders,
#' and aggregates performance metrics (AUC) across all runs. Supports glmnet, SVM, Random Forest and Neural Network models.
#'
#' @param expression_matrix Matrix or dataframe containing features (genes/proteins) as columns and samples as rows
#' @param N Number of repetitions/runs to execute
#' @param seed Initial seed for random number generator (used to generate subsequent seeds). Default 123
#' @param topn Number of top features to select. Default 10
#' @param n_round Number of rounds for internal validation. Default 20
#' @param version Algorithm version to use (1 or 2). Default 2
#' @param train_percent Proportion of data to use for training (0-1). Default 0.7
#' @param group_name Name/numeric index of the outcome variable column. Default 1
#'
#' @return Invisibly returns AUC matrix. Generates:
#' - Subdirectories for each seed
#' - total_auc_result.csv containing:
#'   - AUC values for each model (rows: glmnet/svm/rf/nnet)
#'   - Columns for each seed + average column
#'
#' @examples
#' \dontrun{
#' # Basic usage with 5 repetitions
#' batch_execute_X(expression_matrix = gene_exp, N = 5)
#'
#' # Custom configuration with 3 repetitions
#' batch_execute_X(expression_matrix = protein_data, N = 3,version = 2,
#'                topn = 15, train_percent = 0.8, seed = 456)
#' }
#'
#' @export

batch_execute_X <- function(expression_matrix,
                            N,
                            seed = NULL,
                            topn = 10,
                            n_round = 20,
                            version = 2,
                            train_percent = 0.7,
                            group_name = 1) {

  data_seed = data.frame(matrix(nrow = 4, ncol = N))
  row.names(data_seed) = c("glmnet", "svm", "rf", "nnet")


  random_seed <- seed
  model_list = list()
  # Loop to generate N random seeds and execute function X
  for (i in seq_along(random_seed)) {

    # Generate random seed
    colnames(data_seed)[i] = paste0("seed_", random_seed[i])

    # Create a folder named after the random seed
    seed_folder <- paste0(i, "、seed=", random_seed[i])
    if (!dir.exists(seed_folder)) {
      dir.create(seed_folder)
    }
    setwd(seed_folder)
    if (version == 1) {
      # Call function X and get the result
      result <- machine_learning(data = expression_matrix,
                                 seed = random_seed[i], group_name = group_name,
                                 topn = topn, n_round = n_round,
                                 train_percent = train_percent)
    }else{
      result <- machine_learning2(data = expression_matrix,
                                 seed = random_seed[i], group_name = group_name,
                                 topn = topn, n_round = n_round,
                                 train_percent = train_percent)
    }



    data_seed[, i] = result[[1]]

    model_list[[i]] = list(m1h = result[[2]],svm = result[[3]],rf = result[[4]],nnet = result[[5]])
    # Print progress information
    cat("Completed for seed", random_seed[i], "\n")
    setwd("../")

  }

  data_seed$average = apply(data_seed, 1, mean)
  write.csv(x = data_seed, file = "total_auc_result.csv")

  return(model_list)
}


#' Batch Execution of Cox Model Machine Learning Workflow
#'
#' @description
#' This function performs repeated execution of Cox model machine learning workflows using different random seeds.
#' Results are saved in seed-specific folders. Designed for stability assessment through multiple runs.
#'
#' @param data A dataframe containing survival data (time + status) and features
#' @param N Number of repetitions to execute (should match seed vector length if provided)
#' @param seed Integer or integer vector specifying random seed(s) (length should match N). Default NULL
#' @param topn Number of top features to select. Default 10
#' @param stratify Logical indicating whether to stratify sampling by event status. Default FALSE
#' @param lambda_type Lambda selection type for glmnet ("lambda.min" or "lambda.1se"). Default "lambda.min"
#' @param version Algorithm version to use (1 or 2). Default 2
#' @param n_round Number of rounds for internal validation. Default 20
#' @param bestcut Logical indicating whether to apply optimal cutpoint for features. Default FALSE
#' @param train_percent Proportion of data to use for training (0-1). Default 0.7
#'
#' @return No direct return value. Generates:
#' - Subdirectories for each seed containing:
#'   - data_train_result.csv: Training set results
#'   - data_test_result.csv: Test set results
#'
#' @examples
#' \dontrun{
#' # Run 5 repetitions with automatic seed generation
#' batch_execute_cox(data = tcga_surv, N = 5, topn = 15)
#'
#' # Run 3 repetitions with specific seeds and stratification
#' batch_execute_cox(data = clinical_data, N = 3, seed = c(123, 456, 789),
#'                  stratify = TRUE, lambda_type = "lambda.1se")
#' }
#'
#' @export
batch_execute_cox <- function(data, N, seed = NULL, topn = 10, stratify = F,
                              lambda_type = "lambda.min", version = 2,
                              n_round = 20, bestcut = F,
                              train_percent = 0.7) {

  random_seed <- seed
  # Loop to generate N random seeds and execute function X
  for (i in seq_along(random_seed)) {

    # Create a folder named after the random seed
    seed_folder <- paste0(i, "、seed=", random_seed[i])
    if (!dir.exists(seed_folder)) {
      dir.create(seed_folder)
    }
    setwd(seed_folder)

    # Call function X and get the result
    if (version == 2) {
      result <- machine_learning_cox_2(data = data, bestcut = bestcut,
                                       lambda_type = lambda_type,
                                       # seed = seed,
                                       stratify = stratify,
                                       seed = random_seed[i],
                                       topn = topn, n_round = n_round,
                                       train_percent = train_percent)
    } else {
      result <- machine_learning_cox(data = data, bestcut = bestcut,
                                     lambda_type = lambda_type,
                                     # seed = seed,
                                     stratify = stratify,
                                     seed = random_seed[i],
                                     topn = topn, n_round = n_round,
                                     train_percent = train_percent)
    }

    tryCatch({
      out1 = result[[1]]
      out2 = result[[2]]

      write.csv(x = out1, file = "data_train_result.csv")
      write.csv(x = out2, file = "data_test_result.csv")
    }, error = function(e) {
      # Print error message
      message("An error occurred: ", e$message)
      # Additional error handling logic can be added here
    })

    # Print progress information
    cat("Completed for seed", random_seed[i], "\n")
    setwd("../")
  }
}





