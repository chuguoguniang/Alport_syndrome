XGboost使用——R_二分类_2022.06.29.
与xgboost使用——R.txt相比有如下改进：根据约登法则得最佳概率分界点，进而分类

#加载包
library(xgboost)
library(caret)
library(tidyverse)
library(skimr)
library(DataExplorer)
library(pROC)

#加载数据
setwd("D:/Alport_syndrome/Group5_clinvar_gnomAD_536T_1877F/data2_after_vep/transcript_COL4A5-201/step5_ensemble_vep_pathogenic_state_TOPMED_AF_conserve_state_MAF")
col4a5snv=read.table("16vep.tsv",header=TRUE,sep = "\t")

#数据鸟瞰
skim(col4a5snv)

#数据缺失情况
plot_missing(col4a5snv)
 figure1

#处理缺失
1.删除——使用na.omit()可删除含有缺失值的样本，但对我们的数据集而言此方法的结果剩余的snv过少
2.填充——使用randomForest::na.roughfix(),填充方式是连续型变量用中位数填充，分类型变量用众数填充
我们暂不对缺失做处理

#变量类型修正
显示列名：colnames(col4a5snv)
16个feature：EXON    cDNA_position   CDS_position    Protein_position        Amino_acids     Codons  SIFT    PolyPhen        AF             topmed_AF       jianbiyamu_conserve     lingzhangzongmu_conserve        beifangshoulei_conserve lingzhang_conserve      renke_conserve  xiabixiamu_conserve     Pathogenicity
1个结果：Pathogenicity
使用factor函数修正为因子变量，因子变量基本可理解为离散变量
for (i in c(1,2,3,4,5,6,11,12,13,14,15,16,17)){
    col4a5snv[,i] <- factor(col4a5snv[,i])
}

#数据处理后鸟瞰
skim(col4a5snv) 
table(col4a5snv$EXON)
table(col4a5snv$Pathogenicity) #因变量分布情况

#拆分数据
set.seed(1)
trains <- createDataPartition(y=col4a5snv$Pathogenicity,p=0.6,list=F)
trains2 <- sample(trains,400)
valids <- setdiff(trains,trains2)

data_train <- col4a5snv[trains2,]
data_valid <- col4a5snv[valids,]
data_test <- col4a5snv[-trains,]

#拆分后因变量分布
table(data_train$Pathogenicity)
table(data_valid$Pathogenicity)
table(data_test$Pathogenicity)

#数据准备 
dvfunc <- dummyVars(~.,data = data_train[,1:16],fullRank = T) #构建独热编码模型 ~.表示对数据的所有分类变量处理，fullRank=T排除完全共线性
data_trainx <- predict(dvfunc, newdata = data_train[,1:16])
data_trainy <- ifelse(data_train$Pathogenicity == "No", 0, 1)

data_validx <- predict(dvfunc, newdata = data_valid[,1:16]) ##dvfunc是根据data_train而来
data_validy <- ifelse(data_valid$Pathogenicity == "No", 0, 1)

data_testx <- predict(dvfunc, newdata = data_test[,1:16])
data_testy <- ifelse(data_test$Pathogenicity == "No", 0, 1)

dtrain <- xgb.DMatrix(data = data_trainx,label=data_trainy)
dvalid <- xgb.DMatrix(data = data_validx,label=data_validy)
dtest <- xgb.DMatrix(data = data_testx,label=data_testy)
watchlist <- list(train = dtrain,test = dvalid) #与提前终止有关

#训练模型
fit_xgb_cls <- xgb.train(data = dtrain,eta = 0.3,gamma = 0.001,max_depth = 2,subsample=0.7,colsample_bytree = 0.4, objective = "binary:logistic",nrounds = 1000,watchlist = watchlist,verbose = 1,print_every_n = 100,early_stopping_rounds = 200)

#模型概要
fit_xgb_cls

#变量重要性
importance_matrix <- xgb.importance(model = fit_xgb_cls)
print(importance_matrix)
xgb.plot.importance(importance_matrix = importance_matrix,measure = "Cover")

#SHAP-度量自变量对因变量贡献程度的指标
xgb.plot.shap(data = data_trainx,model = fit_xgb_cls,top_n = 5)

#预测
#训练集预测概率
trainpredprob <- predict(fit_xgb_cls,newdata = dtrain)
trainpredprob

#训练集ROC
trainroc <- roc(response = data_train$Pathogenicity, #实际类别
                         predictor = trainpredprob) #预测概率
#训练集ROC曲线
plot(trainroc,print.auc = TRUE, auc.polygon = TRUE, grid = T, max.auc.polygon = T, auc.polygon.col = "skyblue",print.thres = T, legacy.axes = T, bty = "l")

#约登法则——得最佳概率分界点
bestp <- trainroc$thresholds[which.max(trainroc$sensitivities + trainroc$specificities - 1)]
bestp

#训练集预测分类
trainpredlab <- as.factor(ifelse(trainpredprob > bestp,"Yes","No"))
trainpredlab

#训练集混淆矩阵
confusionMatrix(data = trainpredlab,reference = data_train$Pathogenicity,positive = "Yes",mode = "everything")

#测试集预测概率
testpredprob <- predict(fit_xgb_cls,newdata = dtest)

#测试集预测分类
testpredlab <- as.factor(ifelse(testpredprob > bestp,"Yes","No"))

#测试集混淆矩阵
confusionMatrix(data = testpredlab,reference = data_test$Pathogenicity,positive = "Yes",mode = "everything")

#测试集ROC
testroc <- roc(response = data_test$Pathogenicity, #实际类别
                         predictor = testpredprob) #预测概率

#训练集、测试集ROC曲线叠加
plot(trainroc,
     print.auc = TRUE,
     grid = c(0.1,0.2),
     auc.polygon = F,
     max.auc.polygon = T,
     main = "模型效果图示",
     grid.col=c("green","red") )

plot(testroc,
     print.auc = TRUE,
     print.auc.y = 0.4,
     add = T,
     col = "red")

legend("bottomright",
       legend = c("data_train","data_test"),
       col = c(par("fg"),"red"),
       lwd = 2,
       cex = 0.9)


##使用上述数据集外的点（来源于LOVD）作为额外测试集进行分类预测：
转换变量类型时，将字符型变量转为数值型变量的做法：snv3$topmed_AF <- as.numeric(snv3$topmed_AF)
中间在数据准备过程data_testlovdx <- predict(dvfunc, newdata = snv3[,1:16]) 遇到如下错误：
Error in model.frame.default(Terms, newdata, na.action = na.action, xlev = object$lvls) : 
  因子cDNA_position里出现了新的层次628, 836, 844, 871
可能的原因是：设置train和test时保证test and train 数据集要有相同的因子水平。简言之，各因子水平要同时存在test and train 里。否则你再怎么查资料，重写代码，都无用。
