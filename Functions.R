classify <- function(ts, classes.ts, vs, classes.vs,  
                     classifiers = c("SVM_linear", "SVM_radial", "SVM_polynomial", "SVM_sigmoid", "PAM", "LDA", "KNN", "LASSO", "RTF"), 
                     measures = c("AUC", "train.error", "test.error", "SENS", "SPEC", "ACC"),
                     predictions_table = F, nperm = 1, 
                     opt_cutoff = T, set_cutoff = 0.5){
  
  results <- data.table(classifier = c(rep(classifiers, each = (length(measures)))), 
                        measure = rep(measures, length(classifiers)),
                        value = as.numeric(0), 
                        perm = nperm,
                        number_features = dim(ts)[1], 
                        size.ts = dim(ts)[2])
  
  result.predictions <- data.table(Filename = colnames(vs),
                                   true_class = classes.vs, 
                                   perm = nperm)
  
  print(paste("transpose"))
  ts  <- t(ts)
  vs <- t(vs)
  for (i in 1:length(classifiers)){
    
    if (classifiers[i] == "SVM_linear"){
      print(paste(classifiers[i]))
      model.temp <-  svm(ts, as.factor(classes.ts), kernel="linear", cross=10, probability = TRUE)
      pred.all <- attr(predict(model.temp, vs, probability=T), "probabilities")
      pred.dist <- pred.all[,which(colnames(pred.all)=="X2")]
      train.error <- 1- mean(predict(model.temp, ts) == classes.ts)
    }
    if (classifiers[i] == "SVM_radial"){
      print(paste(classifiers[i]))
      model.temp  <- svm(ts, as.factor(classes.ts), kernel="radial", cross=10, probability = TRUE)
      pred.all <- attr(predict(model.temp, vs, probability=T), "probabilities")
      pred.dist <- pred.all[,which(colnames(pred.all)=="CONTROL")]
      train.error <- 1- mean(predict(model.temp, ts) == classes.ts)
    }
    if (classifiers[i] == "SVM_polynomial"){
      print(paste(classifiers[i]))
      model.temp <- svm(ts, as.factor(classes.ts), kernel="polynomial", cross=10, probability = TRUE)
      pred.all <- attr(predict(model.temp, vs, probability=T), "probabilities")
      pred.dist <- pred.all[,which(colnames(pred.all)=="CONTROL")]
      train.error <- 1- mean(predict(model.temp, ts) == classes.ts)
    }
    if (classifiers[i] == "SVM_sigmoid"){
      print(paste(classifiers[i]))
      model.temp <- svm(ts, as.factor(classes.ts), kernel="sigmoid", cross=10, probability = TRUE)
      pred.all <- attr(predict(model.temp, vs, probability=T), "probabilities")
      pred.dist <- pred.all[,which(colnames(pred.all)=="CONTROL")]
      train.error <- 1- mean(predict(model.temp, ts) == classes.ts)
    }
    if (classifiers[i] == "PAM"){
      print(paste(classifiers[i]))
      model.temp <- pamr.train(list(x = as.matrix(t(ts)), threshold = 0, y = classes.ts))
      pred.all <- pamr.predict(model.temp, as.matrix(t(vs)), 0, type="posterior")
      pred.dist <- pred.all[,which(colnames(pred.all)=="CONTROL")]
      train.error <- 1- mean(pamr.predict(model.temp, as.matrix(t(ts)), 0) == classes.ts)
      
    }
    if (classifiers[i] == "LDA"){
      print(paste(classifiers[i]))
      model.temp <- lda(x = ts, grouping = classes.ts)
      model.pred <- predict(model.temp, vs)$posterior
      pred.dist <- model.pred[,which(colnames(model.pred)=="CONTROL")]
      train.error <- 1- mean(predict(model.temp)$class == classes.ts)
      
    }
    if (classifiers[i] == "KNN"){
      print(paste(classifiers[i]))
      model.temp <- knn(train = ts, test = ts, cl = as.factor(classes.ts), k = 6)
      train.error <- 1- mean(model.temp == classes.ts)
      model.temp <- knn(train = ts, test = vs, cl = as.factor(classes.ts), k = 6, prob = T)
      pred <- attr(model.temp, "prob")
      pred.dist <- ifelse(model.temp == "CASE", 1-pred, pred)
    }
    if (classifiers[i] == "LASSO"){
      print(paste(classifiers[i]))
      model.temp <- cv.glmnet(x = ts, y = as.factor(classes.ts), family = "binomial", alpha = 1) # cv to find optimal s (penaltly parameter lambda)
      train.error <- 1- mean(predict.cv.glmnet(model.temp, newx = ts, type = "class") == classes.ts)
      pred.dist <- predict(model.temp, newx = vs, type = "response", prob = T, s = model.temp$lambda.min) # prediction using optimal lambda
    }
    if (classifiers[i] == "RTF"){
      print(paste(classifiers[i]))
      model.temp <-  randomForest(x = ts, y = as.factor(classes.ts), importance = T)
      train.error <- 1- mean(predict(model.temp, newdata = ts) == classes.ts)
      pred.dist <- predict(model.temp, newdata = vs, type = "prob")[,which(colnames(pred.all)=="CONTROL")]
    }
    result.predictions$probablity_control <- pred.dist
    AUC <- performance(prediction(pred.dist,as.factor(classes.vs)),measure="auc")@y.values[[1]]
    
    if (opt_cutoff == T){
      
      position <- which((performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]] + 
                           performance(prediction(pred.dist,as.factor(classes.vs)),measure="spec")@y.values[[1]]) 
                        == max((performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]] + 
                                  performance(prediction(pred.dist,as.factor(classes.vs)),measure="spec")@y.values[[1]])))
      if (length(position) >1){
        which(performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]][position] 
              == max(performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]][position]))
        position <- position[which(performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]][position] 
                                   == max(performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]][position]))]
        if (length(position) > 1){print("error: too many maxima")}
      }
      
      cutoff <- performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@x.values[[1]][position]
      SENS <- performance(prediction(pred.dist,as.factor(classes.vs)),measure="sens")@y.values[[1]][position]
      SPEC <- performance(prediction(pred.dist,as.factor(classes.vs)),measure="spec")@y.values[[1]][position]
      test.error <- performance(prediction(pred.dist,as.factor(classes.vs)),measure="err")@y.values[[1]][position]
      ACC <- performance(prediction(pred.dist,as.factor(classes.vs)),measure="acc")@y.values[[1]][position]
      
      result.predictions[[classifiers[i]]] <- ifelse(pred.dist >= cutoff, "CONTROL", "CASE") 
    }
    
    else{
      cutoff <- set_cutoff
      result.predictions[[classifiers[i]]] <- ifelse(pred.dist >= cutoff, "CONTROL", "CASE") 
      test.error <- 1 - mean(result.predictions[[classifiers[i]]] == classes.vs)
      P <- sum(result.predictions[[classifiers[i]]] == "CASE") # all positive predictions
      TP <- sum(classes.vs[which(result.predictions[[classifiers[i]]] == "CASE")] == "CASE") # true positives
      FP <- sum(classes.vs[which(result.predictions[[classifiers[i]]] == "CASE")] == "CONTROL") # false positives
      
      N <- sum(result.predictions[[classifiers[i]]] == "CONTROL") # all negative predictions
      TN <- sum(classes.vs[which(result.predictions[[classifiers[i]]] == "CONTROL")] == "CONTROL") # true negatives
      FN <- sum(classes.vs[which(result.predictions[[classifiers[i]]] == "CONTROL")] == "CASE") # false negatives
      
      SENS <- TP/P
      ifelse(N > 0, SPEC <- TN/N, SPEC <- 0)
      ACC <- (TP+TN)/(TP+FP+FN+TN)
      
    }
    
    # get results    
    if ("AUC" %in% measures){results[classifier == classifiers[i] & measure == "AUC"]$value <- AUC}
    if ("train.error" %in% measures){results[classifier == classifiers[i] &  measure == "train.error"]$value <- train.error}
    if ("test.error" %in% measures){results[classifier == classifiers[i] &  measure == "test.error"]$value <- test.error}
    if ("SENS" %in% measures){results[classifier == classifiers[i] & measure == "SENS"]$value <- SENS}
    if ("SPEC" %in% measures){results[classifier == classifiers[i] &  measure == "SPEC"]$value <-  SPEC}    
    if ("ACC" %in% measures){results[classifier == classifiers[i] & measure == "ACC"]$value <-  ACC} 
    results$cutoff <- cutoff
  }
  
  print(paste("returning results"))
  if (predictions_table == T){
    return(list(results, result.predictions))
  } else {return(results)}
}





### Änderung Case Control 1 2
DE_genes <- function(ts, classes.ts = info.training$Condition, 
                     method = "onewayLIMMA", FC = 2, pval = 0.05, adj =T, 
                     batch = info.training$Study)
{
  if (method == "twowayLIMMA"){
    design <- model.matrix(~0 + droplevels(classes.ts) + droplevels(batch))
    colnames(design) <- as.character(c(levels(droplevels(classes.ts)), levels(droplevels(batch))[-1]))
    cm <- makeContrasts(Condition = X1-X2, levels = design) 
    fit <- lmFit(ts, design)# for lmFit genes in rows
    fit <- contrasts.fit(fit, cm)
    fit <- eBayes(fit)
    genes <- topTable(fit,sort="none",n=Inf)
    genes$FC <- logratio2foldchange(genes$logFC)
    
    if (adj == T){
      DE.ts <- rownames(genes[(genes$FC >= FC | genes$FC <= -FC) & genes$adj.P.Val < pval,])
    } else {
      DE.ts <- rownames(genes[(genes$FC >= FC | genes$FC <= -FC) & genes$P.Value < pval,])
    }
  }
  
  if (method == "onewayLIMMA"){
    design <- model.matrix(~0 + droplevels(classes.ts))
    colnames(design) <- as.character(c(levels(classes.ts)))
    cm <- makeContrasts(Condition = X1-X2, levels = design) 
    fit <- lmFit(ts, design)# for lmFit genes in rows
    fit <- contrasts.fit(fit, cm)
    fit <- eBayes(fit)
    genes <- topTable(fit,sort="none",n=Inf)
    genes$FC <- logratio2foldchange(genes$logFC)
    
    if (adj == T){
      DE.ts <- rownames(genes[(genes$FC >= FC | genes$FC <= -FC) & genes$adj.P.Val < pval,])
    } else {
      DE.ts <- rownames(genes[(genes$FC >= FC | genes$FC <= -FC) & genes$P.Value < pval,])
    }
  }
  
  if (method == "ttest"){
    DE.ts <- getDEgenes(ts[,which(classes.ts=="X1")],ts[,which(classes.ts=="X2")],fc=FC ,pval=pval, teststat="t")
  }
  
  if (method == "none"){
    DE.ts <- rownames(data)  
  }
  
  return(DE.ts)
}

library(doParallel)

permute_ts_vs_parallel <- function(nperm=10, 
                                   server = F, 
                                   classifiers. = c("SVM_linear", "SVM_radial", "SVM_polynomial", "SVM_sigmoid", "PAM", "LDA", "KNN", "LASSO", "RTF"),
                                   info = info.training, 
                                   data = data.training.rma.trimmed,
                                   size.ts.case = 50,
                                   size.ts.con = 50,
                                   size.vs.case = 500,
                                   size.vs.con = 500,
                                   dir = "E:/Stefanie/Test/", 
                                   fc = 2, pv = 0.05, adjust = T, 
                                   remove_batch_in_permutations = T, 
                                   opt_cutoff. = T,
                                   print_DEgenes = T,
                                   print_values = T){ 
  ifelse(server == T, cl <- makeCluster(12, outfile = ""), cl <- makeCluster(1, outfile = "E:/Stefanie/Test/test.txt"))
  
  registerDoParallel(cl)
  
  result <- foreach(j = 1:nperm,
                    .combine = rbind,
                    .packages = c("randomForest", "data.table", "class", "e1071", "ROCR", "multtest", "pamr", "MASS", "glmnet", "dplyr", "limma", "randomForest", "gtools"),
                    .export = c("DE_genes", "classify")
  ) %dopar% {
    print(paste("Permutation", j, sep = ":"))
    
    info$Filename <- as.character(info$Filename)
    set.seed(j)
    files.ts <- c(sample_n(info[info$Condition == "CASE",], size = size.ts.case)$Filename, sample_n(info[info$Condition == "CONTROL",], size = size.ts.con)$Filename)
    set.seed(j)
    files.vs <- c(sample_n(info[!info$Filename %in% files.ts & info$Condition == "CASE",], size = size.vs.case)$Filename, 
                  sample_n(info[!info$Filename %in% files.ts & info$Condition == "CONTROL",], size = size.vs.con)$Filename)
    
    ts. <- data[,files.ts]
    vs. <- data[,files.vs]
    
    info.ts. <- info[info$Filename %in% files.ts,]
    info.vs. <- info[info$Filename %in% files.vs,]
    
    info.ts. <- info.ts.[order(rownames(info.ts.)),]
    ts. <- ts.[,order(colnames(ts.))]
    info.vs. <- info.vs.[order(rownames(info.vs.)),]
    vs. <- vs.[,order(colnames(vs.))]
    
    classes.ts. <- info.ts.$Condition
    classes.vs. <- info.vs.$Condition
    
    if(remove_batch_in_permutations == T){
      print("batch removal")
      design <- model.matrix(~0 + as.factor(info.ts.$Study))
      colnames(design) <- levels(as.factor(info.ts.$Study))
      library(limma)
      ts. <- removeBatchEffect(ts., batch = as.factor(as.factor(info.ts.$Study)), design = design) 
    }
    
    probes <- DE_genes(ts = ts., classes.ts = classes.ts., FC = fc, pval = pv, adj = adjust)
    result <- classify(ts = ts.[probes,], classes.ts = classes.ts., vs = vs.[probes,], classes.vs = classes.vs., nperm = j, classifiers = classifiers., opt_cutoff = opt_cutoff.)
    list(result, data.frame(genes = probes, perm = j))
  }
  if(print_values == T){
  setwd(dir)
  res.values <- data.table()
  for (i in 1:nperm){
    res.values <- rbind(res.values, result[[i]])
  }
  write.table(res.values, file = paste("results.values", "Perm", nperm, "sizeTS", size.ts.case+size.ts.con, "FC", fc, "txt", sep = "."), quote = F, sep = "\t")

  rm(res.values)
  }
  if(print_DEgenes == T){
    res.DEgenes <- data.table()
  for (i in 1:nperm){
   res.DEgenes <- rbind(res.DEgenes, result[[nperm + i]])
  }
  write.table(res.DEgenes, file = paste("results.DEgenes", "Perm", nperm, "sizeTS", size.ts.case+size.ts.con,"FC", fc, "txt", sep = "."), quote = F, sep = "\t")
  }
  stopCluster(cl)
}



