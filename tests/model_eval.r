
model_eval <- function(feature.ranking, important.features){


        #browser()
        total.dfs <- length(feature.ranking)
        total.features <- nrow(feature.ranking[[1]])
        hits <- sum(sapply(feature.ranking, function(df) sum(as.integer(rownames(df)) %in% important.features)))
        

        tpr <- hits / (total.dfs * length(important.features))
        fpr <- (total.dfs * length(important.features) - hits) / (total.dfs * length(important.features))
        fnr <- (total.dfs * length(important.features) - hits) / (total.dfs * length(important.features))

        precision <- tpr / (tpr + fpr)
        recall <- tpr / (tpr + fnr)

        result <- data.frame(precision = precision, recall = recall)
        return(result)

}