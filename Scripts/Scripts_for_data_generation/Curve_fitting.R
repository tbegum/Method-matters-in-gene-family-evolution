## This script is written to test, and avoid curve overfitting 
## This example shows the test on time dependent trait acceleration model

## Loading previously stored data
## for pdup=0.5
load("timemodel_PW_pdup0.5_sig0.02_plot_data_latest.rda")

## loading functions required
source("functions_TMM.R")
set.seed(1234)

###Checking on data with trait acceleration after 70% time since old duplication
data<-processed_corr_0.5_oc0.9times
data$event<- factor(data$Event, levels=c("speciation", "duplication" ) )

##Loading libraries
library(gnm)
library(ggplot2)

##Plotting data with different degrees
plot(data$R ~ data$Mya, col = "gray", lwd = 2)
## linear model
model1<-lm(R~poly(Mya,1,raw=TRUE), data=data)
model2<-lm(R~poly(Mya,2,raw=TRUE), data=data)
model3<-lm(R~poly(Mya,3,raw=TRUE), data=data)
lines(model1$fitted.values ~ data$Mya, lwd = 2, col = "palegreen3")
lines(model2$fitted.values ~ data$Mya, lwd = 2, col = "orange")
lines(model3$fitted.values ~ data$Mya, lwd = 2, col = "steelblue")
legend(x = "topright", legend = c("Linear fit", "Quadratic fit", 
                                  "Cubic fit"), lwd = rep(3, 4), col = c( "palegreen3","orange", "steelblue"))

## Now to solve the overfitting issue
##Now for each dataset spliting it into "training" (70%) and "validation" (30%)
#split<-sample(2,nrow(data),replace=F,prob = c(0.7,0.3))
#training_set<-data[split==1,]
#test_set<-data[split==2,]
set.seed(1234)
split<-sample(nrow(data),size=(0.7*nrow(data)),replace=F)
training_set<-data[split,]
test_set<-data[-split,]

## R function to get complexity
get_complexity <- function(model) {
  length(coef(model)) - 1
}

## Now performing polynomial regressions with training data
Regression_training_fit1<-gnm(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data = training_set) 
summary(Regression_training_fit1)
get_complexity(Regression_training_fit1) # more complex model, complexity=5

Regression_training_fit2<-gnm(R~ Mya + I(Mya^2) + Event + Mya*Event, data = training_set)
summary(Regression_training_fit2)
get_complexity(Regression_training_fit2) #complexity=4

Regression_training_fit3<-gnm(R~ Mya + Event + Mya*Event, data = training_set)
summary(Regression_training_fit3)
get_complexity(Regression_training_fit3) #complexity=3

## R function for RMSE or more improved get_rmse
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

get_rmse <- function(model, data, response) {
  rmse(actual = data[, response], 
       predicted = predict(model, data))
}


## model list
model_list <- list(Regression_training_fit1, Regression_training_fit2, Regression_training_fit3)

# train RMSE for each model
#sqrt(mean((training_set$R - predict(Regression_training, training_set)) ^ 2))
train_rmse <- sapply(model_list, get_rmse, data = training_set, response = "R")

# test RMSE for each model
#sqrt(mean((test_set$R - predict(Regression_training, test_set)) ^ 2))
test_rmse <- sapply(model_list, get_rmse, data = test_set, response = "R")

## model complexity of each model
model_complexity <- sapply(model_list, get_complexity)

## model complexity of each model
degree_polynomial <- c(3,2,1)

##Plotting
plot(degree_polynomial, train_rmse, type = "b", 
     ylim = c(min(c(train_rmse, test_rmse)) - 0.02, 
              max(c(train_rmse, test_rmse)) + 0.02), 
     col = "dodgerblue", 
     xlab = "Degree of polynomial",
     ylab = "RMSE")
lines(degree_polynomial, test_rmse, type = "b", col = "darkorange")
#legend(2.7,0.12,legend=c("training_set","test_set"), pch=19, col=c("dodgerblue","darkorange"))
legend("topright",legend=c("training_set","test_set"), pch=19, col=c("dodgerblue","darkorange"))


## Specifically, we say that a model is overfitting if there exists a less complex model with lower Test RMSE. 
#Then a model is underfitting if there exists a more complex model with lower Test RMSE.
## Prediction
#prediction<-predict(Regression_training_fit3,test_set)
##plotting
#plot(test_set$Mya,prediction)
#lines(test_set$Mya,test_set$R, col="red")


################### 10 fold cross validation #########################
#split<-0.70
#trainIndex <- createDataPartition(data$Event, p=split, list=FALSE)
#data_train <- data[ trainIndex,]
#data_test <- data[-trainIndex,]
set.seed(1234)
library(caret)
# define training control
train_control <- trainControl(method="cv", number=10)

# train the model
modelt1 <- train(R~ Mya  + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt2 <- train(R~ Mya + I(Mya^2) + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt3 <- train(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data=data, trControl=train_control, method="lm")

# summarize results
print(modelt1)
print(modelt2)
print(modelt3)

dtfrm<-data.frame(model=NA,RMSE=NA,Rsquared=NA)
dtfrm<-rbind(dtfrm, data.frame(model="linear",RMSE=modelt1$results$RMSE,Rsquared=modelt1$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="quadratic",RMSE=modelt2$results$RMSE,Rsquared=modelt2$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="cubic",RMSE=modelt3$results$RMSE,Rsquared=modelt3$results$Rsquared))
dtfrm<-dtfrm[-1,]

dtfrm


################### repeated 10 fold cross validation #########################
#split<-0.70
#trainIndex <- createDataPartition(data$Event, p=split, list=FALSE)
#data_train <- data[ trainIndex,]
#data_test <- data[-trainIndex,]
set.seed(1234)
# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats=10)

# train the model
modelt1 <- train(R~ Mya  + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt2 <- train(R~ Mya + I(Mya^2) + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt3 <- train(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data=data, trControl=train_control, method="lm")

# summarize results
print(modelt1)
print(modelt2)
print(modelt3)

dtfrm<-data.frame(model=NA,RMSE=NA,Rsquared=NA)
dtfrm<-rbind(dtfrm, data.frame(model="linear",RMSE=modelt1$results$RMSE,Rsquared=modelt1$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="quadratic",RMSE=modelt2$results$RMSE,Rsquared=modelt2$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="cubic",RMSE=modelt3$results$RMSE,Rsquared=modelt3$results$Rsquared))
dtfrm<-dtfrm[-1,]

dtfrm ## This result shows that cubic model is the best fit model


########## This example shows the test on the trait jump model #################

## Loading previously stored data
## for pdup=0.5
load("Asym_jump_10000_0.5dup_0.02sig2_pw_plot_data_latest.rda")

## loading functions required
source("functions_TMM.R")
set.seed(1234)

###Checking on data with trait acceleration after 70% time since old duplication
data<-processed_corr_0.5_OC
data$event<- factor(data$Event, levels=c("speciation", "duplication" ) )

##Loading libraries
library(gnm)
library(ggplot2)

##Plotting data with different degrees
plot(data$R ~ data$Mya, col = "gray", lwd = 2)
## linear model
model1<-lm(R~poly(Mya,1,raw=TRUE), data=data)
model2<-lm(R~poly(Mya,2,raw=TRUE), data=data)
model3<-lm(R~poly(Mya,3,raw=TRUE), data=data)
lines(model1$fitted.values ~ data$Mya, lwd = 2, col = "palegreen3")
lines(model2$fitted.values ~ data$Mya, lwd = 2, col = "orange")
lines(model3$fitted.values ~ data$Mya, lwd = 2, col = "steelblue")
legend(x = "topright", legend = c("Linear fit", "Quadratic fit", 
                                  "Cubic fit"), lwd = rep(3, 4), col = c( "palegreen3","orange", "steelblue"))

## Now to solve the overfitting issue
##Now for each dataset spliting it into "training" (70%) and "validation" (30%)
#split<-sample(2,nrow(data),replace=F,prob = c(0.7,0.3))
#training_set<-data[split==1,]
#test_set<-data[split==2,]
set.seed(1234)
split<-sample(nrow(data),size=(0.7*nrow(data)),replace=F)
training_set<-data[split,]
test_set<-data[-split,]

## R function to get complexity
get_complexity <- function(model) {
  length(coef(model)) - 1
}

## Now performing polynomial regressions with training data
Regression_training_fit1<-gnm(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data = training_set) 
summary(Regression_training_fit1)
get_complexity(Regression_training_fit1) # more complex model, complexity=5

Regression_training_fit2<-gnm(R~ Mya + I(Mya^2) + Event + Mya*Event, data = training_set)
summary(Regression_training_fit2)
get_complexity(Regression_training_fit2) #complexity=4

Regression_training_fit3<-gnm(R~ Mya + Event + Mya*Event, data = training_set)
summary(Regression_training_fit3)
get_complexity(Regression_training_fit3) #complexity=3

## R function for RMSE or more improved get_rmse
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

get_rmse <- function(model, data, response) {
  rmse(actual = data[, response], 
       predicted = predict(model, data))
}


## model list
model_list <- list(Regression_training_fit1, Regression_training_fit2, Regression_training_fit3)

# train RMSE for each model
#sqrt(mean((training_set$R - predict(Regression_training, training_set)) ^ 2))
train_rmse <- sapply(model_list, get_rmse, data = training_set, response = "R")

# test RMSE for each model
#sqrt(mean((test_set$R - predict(Regression_training, test_set)) ^ 2))
test_rmse <- sapply(model_list, get_rmse, data = test_set, response = "R")

## model complexity of each model
model_complexity <- sapply(model_list, get_complexity)

## model complexity of each model
degree_polynomial <- c(3,2,1)

##Plotting
plot(degree_polynomial, train_rmse, type = "b", 
     ylim = c(min(c(train_rmse, test_rmse)) - 0.02, 
              max(c(train_rmse, test_rmse)) + 0.02), 
     col = "dodgerblue", 
     xlab = "Degree of polynomial",
     ylab = "RMSE")
lines(degree_polynomial, test_rmse, type = "b", col = "darkorange")
#legend(2.7,0.12,legend=c("training_set","test_set"), pch=19, col=c("dodgerblue","darkorange"))
legend("topright",legend=c("training_set","test_set"), pch=19, col=c("dodgerblue","darkorange"))


## Specifically, we say that a model is overfitting if there exists a less complex model with lower Test RMSE. 
#Then a model is underfitting if there exists a more complex model with lower Test RMSE.
## Prediction
#prediction<-predict(Regression_training_fit3,test_set)
##plotting
#plot(test_set$Mya,prediction)
#lines(test_set$Mya,test_set$R, col="red")


################### 10 fold cross validation #########################
#split<-0.70
#trainIndex <- createDataPartition(data$Event, p=split, list=FALSE)
#data_train <- data[ trainIndex,]
#data_test <- data[-trainIndex,]
set.seed(1234)
library(caret)
# define training control
train_control <- trainControl(method="cv", number=10)

# train the model
modelt1 <- train(R~ Mya  + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt2 <- train(R~ Mya + I(Mya^2) + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt3 <- train(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data=data, trControl=train_control, method="lm")

# summarize results
print(modelt1)
print(modelt2)
print(modelt3)

dtfrm<-data.frame(model=NA,RMSE=NA,Rsquared=NA)
dtfrm<-rbind(dtfrm, data.frame(model="linear",RMSE=modelt1$results$RMSE,Rsquared=modelt1$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="quadratic",RMSE=modelt2$results$RMSE,Rsquared=modelt2$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="cubic",RMSE=modelt3$results$RMSE,Rsquared=modelt3$results$Rsquared))
dtfrm<-dtfrm[-1,]

dtfrm


################### repeated 10 fold cross validation #########################
#split<-0.70
#trainIndex <- createDataPartition(data$Event, p=split, list=FALSE)
#data_train <- data[ trainIndex,]
#data_test <- data[-trainIndex,]
set.seed(1234)
# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats=10)

# train the model
modelt1 <- train(R~ Mya  + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt2 <- train(R~ Mya + I(Mya^2) + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt3 <- train(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data=data, trControl=train_control, method="lm")

# summarize results
print(modelt1)
print(modelt2)
print(modelt3)

dtfrm<-data.frame(model=NA,RMSE=NA,Rsquared=NA)
dtfrm<-rbind(dtfrm, data.frame(model="linear",RMSE=modelt1$results$RMSE,Rsquared=modelt1$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="quadratic",RMSE=modelt2$results$RMSE,Rsquared=modelt2$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="cubic",RMSE=modelt3$results$RMSE,Rsquared=modelt3$results$Rsquared))
dtfrm<-dtfrm[-1,]

dtfrm ## Cubic model is the best fit model


########## This example shows the test on the rates of sequence and trait evolution model #################

## Loading previously stored data
## for pdup=0.5
load("Asym_PW_pdup0.5_sig0.02_plot_data_latest.rda")

## loading functions required
source("functions_TMM.R")
set.seed(1234)

###Checking on data with trait acceleration after 70% time since old duplication
data<-processed_corr_0.5_OC3
data$event<- factor(data$Event, levels=c("speciation", "duplication" ) )

##Loading libraries
library(gnm)
library(ggplot2)

##Plotting data with different degrees
plot(data$R ~ data$Mya, col = "gray", lwd = 2)
## linear model
model1<-lm(R~poly(Mya,1,raw=TRUE), data=data)
model2<-lm(R~poly(Mya,2,raw=TRUE), data=data)
model3<-lm(R~poly(Mya,3,raw=TRUE), data=data)
lines(model1$fitted.values ~ data$Mya, lwd = 2, col = "palegreen3")
lines(model2$fitted.values ~ data$Mya, lwd = 2, col = "orange")
lines(model3$fitted.values ~ data$Mya, lwd = 2, col = "steelblue")
legend(x = "topright", legend = c("Linear fit", "Quadratic fit", 
                                  "Cubic fit"), lwd = rep(3, 4), col = c( "palegreen3","orange", "steelblue"))

## Now to solve the overfitting issue
##Now for each dataset spliting it into "training" (70%) and "validation" (30%)
#split<-sample(2,nrow(data),replace=F,prob = c(0.7,0.3))
#training_set<-data[split==1,]
#test_set<-data[split==2,]
set.seed(1234)
split<-sample(nrow(data),size=(0.7*nrow(data)),replace=F)
training_set<-data[split,]
test_set<-data[-split,]

## R function to get complexity
get_complexity <- function(model) {
  length(coef(model)) - 1
}

## Now performing polynomial regressions with training data
Regression_training_fit1<-gnm(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data = training_set) 
summary(Regression_training_fit1)
get_complexity(Regression_training_fit1) # more complex model, complexity=5

Regression_training_fit2<-gnm(R~ Mya + I(Mya^2) + Event + Mya*Event, data = training_set)
summary(Regression_training_fit2)
get_complexity(Regression_training_fit2) #complexity=4

Regression_training_fit3<-gnm(R~ Mya + Event + Mya*Event, data = training_set)
summary(Regression_training_fit3)
get_complexity(Regression_training_fit3) #complexity=3

## R function for RMSE or more improved get_rmse
rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

get_rmse <- function(model, data, response) {
  rmse(actual = data[, response], 
       predicted = predict(model, data))
}


## model list
model_list <- list(Regression_training_fit1, Regression_training_fit2, Regression_training_fit3)

# train RMSE for each model
#sqrt(mean((training_set$R - predict(Regression_training, training_set)) ^ 2))
train_rmse <- sapply(model_list, get_rmse, data = training_set, response = "R")

# test RMSE for each model
#sqrt(mean((test_set$R - predict(Regression_training, test_set)) ^ 2))
test_rmse <- sapply(model_list, get_rmse, data = test_set, response = "R")

## model complexity of each model
model_complexity <- sapply(model_list, get_complexity)

## model complexity of each model
degree_polynomial <- c(3,2,1)

##Plotting
plot(degree_polynomial, train_rmse, type = "b", 
     ylim = c(min(c(train_rmse, test_rmse)) - 0.02, 
              max(c(train_rmse, test_rmse)) + 0.02), 
     col = "dodgerblue", 
     xlab = "Degree of polynomial",
     ylab = "RMSE")
lines(degree_polynomial, test_rmse, type = "b", col = "darkorange")
#legend(2.7,0.12,legend=c("training_set","test_set"), pch=19, col=c("dodgerblue","darkorange"))
legend("topright",legend=c("training_set","test_set"), pch=19, col=c("dodgerblue","darkorange"))


## Specifically, we say that a model is overfitting if there exists a less complex model with lower Test RMSE. 
#Then a model is underfitting if there exists a more complex model with lower Test RMSE.
## Prediction
#prediction<-predict(Regression_training_fit3,test_set)
##plotting
#plot(test_set$Mya,prediction)
#lines(test_set$Mya,test_set$R, col="red")


################### 10 fold cross validation #########################
#split<-0.70
#trainIndex <- createDataPartition(data$Event, p=split, list=FALSE)
#data_train <- data[ trainIndex,]
#data_test <- data[-trainIndex,]
set.seed(1234)
library(caret)
# define training control
train_control <- trainControl(method="cv", number=10)

# train the model
modelt1 <- train(R~ Mya  + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt2 <- train(R~ Mya + I(Mya^2) + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt3 <- train(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data=data, trControl=train_control, method="lm")

# summarize results
print(modelt1)
print(modelt2)
print(modelt3)

dtfrm<-data.frame(model=NA,RMSE=NA,Rsquared=NA)
dtfrm<-rbind(dtfrm, data.frame(model="linear",RMSE=modelt1$results$RMSE,Rsquared=modelt1$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="quadratic",RMSE=modelt2$results$RMSE,Rsquared=modelt2$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="cubic",RMSE=modelt3$results$RMSE,Rsquared=modelt3$results$Rsquared))
dtfrm<-dtfrm[-1,]

dtfrm


################### repeated 10 fold cross validation #########################
#split<-0.70
#trainIndex <- createDataPartition(data$Event, p=split, list=FALSE)
#data_train <- data[ trainIndex,]
#data_test <- data[-trainIndex,]
set.seed(1234)
# define training control
train_control <- trainControl(method="repeatedcv", number=10, repeats=10)

# train the model
modelt1 <- train(R~ Mya  + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt2 <- train(R~ Mya + I(Mya^2) + Event + Mya*Event, data=data, trControl=train_control, method="lm")
modelt3 <- train(R~ Mya + I(Mya^2) +  I(Mya^3) + Event + Mya*Event, data=data, trControl=train_control, method="lm")

# summarize results
print(modelt1)
print(modelt2)
print(modelt3)

dtfrm<-data.frame(model=NA,RMSE=NA,Rsquared=NA)
dtfrm<-rbind(dtfrm, data.frame(model="linear",RMSE=modelt1$results$RMSE,Rsquared=modelt1$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="quadratic",RMSE=modelt2$results$RMSE,Rsquared=modelt2$results$Rsquared))
dtfrm<-rbind(dtfrm, data.frame(model="cubic",RMSE=modelt3$results$RMSE,Rsquared=modelt3$results$Rsquared))
dtfrm<-dtfrm[-1,]

dtfrm ## Quadratic model is the best fit model


