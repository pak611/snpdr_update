
model_eval <- function(feature.ranking, important.features, num.sig){


        #browser()


        # Define the function
        calc_metrics <- function(df) {
        #browser()
        hits <- sum((rownames(df[1:num.sig, ]) %in% important.features))
        misses <- sum(!(rownames(df[1:num.sig, ]) %in% important.features))
        true.negatives <- sum(!(rownames(df[num.sig:nrow(df), ]) %in% important.features))
        false.negatives <- sum((rownames(df[num.sig:nrow(df), ]) %in% important.features))

        tpr <- hits / length(important.features)
        fpr <- misses / length(important.features)
        fnr <- false.negatives / (nrow(df) - num.sig)
        tnr <- true.negatives / (nrow(df) - num.sig)

        precision <- tpr / (tpr + fpr)
        recall <- tpr / (tpr + fnr)

        # calculate the matthews correlation coefficient (MCC)
        mcc <- ((hits * true.negatives) - (misses * false.negatives)) / sqrt((hits + misses) * (hits + false.negatives) * (true.negatives + misses) * (true.negatives + false.negatives))

        # Return a dataframe with the results
        data.frame(precision = precision, recall = recall, mcc = mcc)
        }

        # Apply the function to each dataframe in feature.ranking
        results <- lapply(feature.ranking, calc_metrics)

        # Bind the dataframes together
        df <- do.call(rbind, results)

        # Calculate the column-wise means
        results <- colMeans(df)

        return(results)

}

model_eval2 <- function(feature.ranking, important.features, num.sig) {
        #browser()


        # Define the function
        calc_metrics <- function(df) {
        #browser()
        predicted.important <- rownames(df[df$predicted_class == 1, ])
        predicted.unimportant <- rownames(df[df$predicted_class == 0, ])
        actual.important <- important.features
        actual.unimportant <- setdiff(rownames(df), important.features)


        hits <- ifelse(sum(predicted.important %in% actual.important) == 0, 0.0001, sum(predicted.important %in% actual.important))
        misses <- ifelse(sum(predicted.important %in% actual.unimportant) == 0, 0.0001, sum(predicted.important %in% actual.unimportant))

        true.negatives <- ifelse(sum(predicted.unimportant %in% actual.unimportant) == 0, 0.0001, sum(predicted.unimportant %in% actual.unimportant))
        false.negatives <- ifelse(sum(predicted.unimportant %in% actual.important) == 0, 0.0001, sum(predicted.unimportant %in% actual.important))

        tpr <- ifelse(hits / length(important.features) == 0, 0.0001, hits / length(important.features))
        fpr <- ifelse(misses / length(important.features) == 0, 0.0001, misses / length(important.features))
        fnr <- ifelse(false.negatives / (nrow(df) - num.sig) == 0, 0.0001, false.negatives / (nrow(df) - num.sig))
        tnr <- ifelse(true.negatives / (nrow(df) - num.sig) == 0, 0.0001, true.negatives / (nrow(df) - num.sig))

        precision <- tpr / (tpr + fpr)
        recall <- tpr / (tpr + fnr)

        # calculate the matthews correlation coefficient (MCC)
        mcc <- ((hits * true.negatives) - (misses * false.negatives)) / sqrt((hits + misses) * (hits + false.negatives) * (true.negatives + misses) * (true.negatives + false.negatives))

        # Return a dataframe with the results
        data.frame(precision = precision, recall = recall, mcc = mcc)
        }

        # Apply the function to each dataframe in feature.ranking
        results <- lapply(feature.ranking, calc_metrics)

        # Bind the dataframes together
        df <- do.call(rbind, results)


        return(df)

}


