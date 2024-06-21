

calculate_auRC <- function(results_df, threshold_range) {

    #browser()
    results_class <- lapply(threshold_range, function(threshold_value) {
        results <- lapply(results_df, function(df) {
            #browser()
            df$predicted_class <- ifelse(abs(df$beta.Z.att) > threshold_value, 1, 0)
            df$actual_class <- ifelse(rownames(df) %in% seq(1, num.sig), 1, 0)
            return(df)
        })
        results_df <- do.call(rbind, results)
        return(results_df)
    })

    predicted_class <- unlist(lapply(results_class, function(df) df$predicted_class))
    actual_class <- unlist(lapply(results_class, function(df) df$actual_class))

    scores <- predicted_class
    labels <- actual_class

    return(data.frame(scores = scores, labels = labels))

}


