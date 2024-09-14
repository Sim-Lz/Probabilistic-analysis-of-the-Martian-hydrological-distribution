rm(list=ls())

#-------------------------------------------------------------------#
#                 Import files and libraries                        #
#-------------------------------------------------------------------#

#install needed libraries 
install.packages("dplyr")
install.packages("ggplot2")
install.packages("DataExplorer")
install.packages("tidyverse")
install.packages("SmartEDA")
install.packages("dlookr")
install.packages("leaps")
install.packages("car")
install.packages("tidyr")
install.packages("heatmaply")
install.packages("ROSE")
install.packages("boot")
install.packages("glmnet")
install.packages("ROCR")
install.packages("reshape2")

#read libraries
library("dplyr")
library("ggplot2")
library("DataExplorer")
library("tidyverse")
library("SmartEDA")
library("dlookr")
library("leaps")
library("car")
library("tidyr")
library("heatmaply")
library("ROSE")
library("boot")
library("glmnet")
library("ROCR")
library("reshape2")

#Import the csv file containing our data
Variables<-read.csv("Variables data.csv")

#-------------------------------------------------------------------#
#                  Esploration data analysis                        #
#-------------------------------------------------------------------#

#First step is get a report without a responding variable
diagnose_web_report(Variables) 

#EDA with plots of variables
ExpReport(Variables, op_file = "smartEDA.html")

#EDA with plots of variables and basic information on the dataset
create_report(Variables)

Variables<-subset(Variables, select=-Altitude.above.local.surface..m.) #It is constant, so it is not usable

#-------------------------------------------------------------------#
#   Check if any single variable satisfies normality assumption     #
#                           Shapiro test                            #
#-------------------------------------------------------------------#

#Definition of Normal/Not Normal classification function
Normal_assumption<- function(shapiro_test) {
  classe <- ifelse(shapiro_test$p.value < 2.2e-16, "Not Normal", "Normal")
  
  return(classe)  
}

#Temperature
Normal_assumption(shapiro.test(Variables$Temperature..K.))#Not Normal

#Pressure
Normal_assumption(shapiro.test(Variables$Pression..Pa.)) #Not Normal

#Monthly mean surface H2O layer
Normal_assumption(shapiro.test(Variables$Monthly.mean.surface.H2O.layer..kg.m2.)) #Not Normal

#Monthly mean surface CO2 ice layer
Normal_assumption(shapiro.test(Variables$Monthly.mean.surface.CO2.ice.layer..kg.m2.)) #Not Normal

#Water vapor column
Normal_assumption(shapiro.test(Variables$Water.vapor.column..kg.m2.)) #Not Normal

#Water ice column
Normal_assumption(shapiro.test(Variables$Water.ice.column..kg.m2.)) #Not Normal

#Water ice mixing ratio
Normal_assumption(shapiro.test(Variables$Water.ice.mixing.ratio..mol.mol.)) #Not Normal

#Water ice effective radius
Normal_assumption(shapiro.test(Variables$Water.ice.effective.radius..m.)) #Not Normal

#GCM perennial surface water ice
Normal_assumption(shapiro.test(Variables$GCM.perennial.surface.water.ice..0.or.1.)) #Not Normal

#H column
Normal_assumption(shapiro.test(Variables$H.column..kg.m2.)) #Not Normal

#H2 column
Normal_assumption(shapiro.test(Variables$H2.column..kg.m2.)) #Not Normal

#Electron number density
Normal_assumption(shapiro.test(Variables$Electron.number.density..m.3.))  #Not Normal

#Total electronic content 
Normal_assumption(shapiro.test(Variables$Total.electronic.content)) #Not Normal

#All not normal distributed

#-------------------------------------------------------------------#
#Check if any single variable satisfies homoscedasticity assumption #
#                            Flinger test                           #
#-------------------------------------------------------------------#

Homoscedasticity_assumption<- function(fligner.test) {
  classe <- ifelse(fligner.test$p.value < 2.2e-16, "Not Homoscedastic", "Homoscedastic")
  
  return(classe)  
}

#Temperature
Homoscedasticity_assumption(fligner.test(Variables$Temperature..K. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# Pressure
Homoscedasticity_assumption(fligner.test(Variables$Pression..Pa. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# Monthly mean surface H2O layer
Homoscedasticity_assumption(fligner.test(Variables$Monthly.mean.surface.H2O.layer..kg.m2. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Homoscedastic

# Monthly mean surface CO2 layer
Homoscedasticity_assumption(fligner.test(Variables$Monthly.mean.surface.CO2.ice.layer..kg.m2. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Homoscedastic

# Water vapor column
Homoscedasticity_assumption(fligner.test(Variables$Water.vapor.column..kg.m2. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# Water ice column
Homoscedasticity_assumption(fligner.test(Variables$Water.ice.column..kg.m2. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Homoscedastic

# Water ice mixing ratio
Homoscedasticity_assumption(fligner.test(Variables$Water.ice.mixing.ratio..mol.mol. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Homoscedastic

# Water ice effective radius
Homoscedasticity_assumption(fligner.test(Variables$Water.ice.effective.radius..m. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Homoscedastic

# H column
Homoscedasticity_assumption(fligner.test(Variables$H.column..kg.m2. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# H2 column
Homoscedasticity_assumption(fligner.test(Variables$H2.column..kg.m2. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# Elecrton number density
Homoscedasticity_assumption(fligner.test(Variables$Electron.number.density..m.3. ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# Total electronic content
Homoscedasticity_assumption(fligner.test(Variables$Total.electronic.content ~ Variables$GCM.perennial.surface.water.ice..0.or.1., data = Variables))
# Not Homoscedastic

# Monthly mean surface H2O layer, Monthly mean surface CO2 layer, Water ice column, Water ice 
# mixing ratio, Water ice effective radius are Homoscedastic variables (p-value higher than 2.2e-16)

#-------------------------------------------------------------------#
#                  Forward stepwise selection                       #
#-------------------------------------------------------------------#

# Best Subset Selection
#check the book for explaination
regfit.full <- regsubsets(`GCM.perennial.surface.water.ice..0.or.1.` ~., data = Variables)
summary(regfit.full)

regfit.full <- regsubsets(`GCM.perennial.surface.water.ice..0.or.1.` ~., data = Variables,nvmax=12)
reg.summary<-summary(regfit.full)

reg.summary$rsq

par(mfrow=c(2,2))
plot(reg.summary$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
plot(reg.summary$adjr2 ,xlab="Number of Variables ",ylab="Adjusted RSq",type="l")

which.max(reg.summary$adjr2)
points(9,reg.summary$adjr2[9], col="red",cex=2,pch=20)

plot(reg.summary$cp, xlab="Number of Variables ", ylab="Cp", type="l")
which.min(reg.summary$cp)
points(9,reg.summary$cp [9],col="red",cex=2,pch=20)

plot(reg.summary$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
which.min(reg.summary$bic )
points(7,reg.summary$bic [7],col="red",cex=2,pch=20)

?plot.regsubsets

plot(regfit.full,scale="r2")
plot(regfit.full,scale="adjr2") 
plot(regfit.full,scale="Cp")
plot(regfit.full,scale="bic")

coef(regfit.full ,7)

Variables.full<- subset(Variables, select = c(GCM.perennial.surface.water.ice..0.or.1.,Temperature..K.,Pression..Pa.,Monthly.mean.surface.H2O.layer..kg.m2.,Water.ice.column..kg.m2.,Water.ice.mixing.ratio..mol.mol.,Water.ice.effective.radius..m.,H.column..kg.m2.))





# Forward stepwise selection
regfit.fwd<- regsubsets(`GCM.perennial.surface.water.ice..0.or.1.` ~.,data=Variables, nvmax=12, method ="forward")
summary(regfit.fwd)

par(mfrow=c(2,2))
plot(summary(regfit.fwd)$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
plot(summary(regfit.fwd)$adjr2 ,xlab="Number of Variables ",ylab="Adjusted RSq",type="l")

which.max(summary(regfit.fwd)$adjr2)
points(9,summary(regfit.fwd)$adjr2[9], col="red",cex=2,pch=20)

plot(summary(regfit.fwd)$cp, xlab="Number of Variables ", ylab="Cp", type="l")
which.min(summary(regfit.fwd)$cp)
points(9,summary(regfit.fwd)$cp [9],col="red",cex=2,pch=20)

plot(summary(regfit.fwd)$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
which.min(summary(regfit.fwd)$bic )
points(8,summary(regfit.fwd)$bic [8],col="red",cex=2,pch=20)

?plot.regsubsets

plot(regfit.fwd,scale="r2")
plot(regfit.fwd,scale="adjr2") 
plot(regfit.fwd,scale="Cp")
plot(regfit.fwd,scale="bic")

coef(regfit.fwd, 8) 

Variables.fwd<- subset(Variables, select = c(GCM.perennial.surface.water.ice..0.or.1.,Temperature..K.,Pression..Pa.,Monthly.mean.surface.H2O.layer..kg.m2.,Water.ice.column..kg.m2.,Water.ice.mixing.ratio..mol.mol.,Water.ice.effective.radius..m.,H.column..kg.m2.,H2.column..kg.m2.))



#-------------------------------------------------------------------#
#                       Logistic regression                         #
#                           first part                              #
#-------------------------------------------------------------------#

# Logistic regression model on all our data
glm.fit <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Monthly.mean.surface.CO2.ice.layer..kg.m2.+Water.vapor.column..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2.+Electron.number.density..m.3.+Total.electronic.content, data = Variables, family = binomial())
summary(glm.fit)
coef(glm.fit)


# Logistic regression model on the subset created on forward stepwise selection
glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd, family = binomial())
summary(glm.fit.fwd)
coef(glm.fit.fwd)

# Logistic regression model on the subset created on best subset selection
glm.fit.full <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2., data = Variables.full, family = binomial())
summary(glm.fit.full)
coef(glm.fit.full)


#-------------------------------------------------------------------#
#                    Warning resolution attempts                    #
#-------------------------------------------------------------------#

# Standardization of data and repetition of the forward stepwise selection           

# Define standardization function 
standardization <- function(data) {
  mean <- colMeans(data, na.rm = TRUE)
  sd <- apply(data, 2, sd, na.rm = TRUE)
  
  for (col in colnames(data)) {
    if (is.numeric(data[[col]])) {
      data[[col]] <- (data[[col]] - mean[col]) / sd[col]
    }
  }
  return(data)
}

Variables.strd <- standardization(Variables)

# We didn't manage to reverse the standardization for our dependent variable, so we decided to 
# remove it from the standardized dataframe and re-add it from the original (non-standardized) one
Variables.strd <- subset(Variables.strd, select = -GCM.perennial.surface.water.ice..0.or.1.) 
Variables.strd$GCM.perennial.surface.water.ice..0.or.1. <- Variables$GCM.perennial.surface.water.ice..0.or.1.

# Forward stepwise selection after standardization

regfit.fwd.strd<- regsubsets(`GCM.perennial.surface.water.ice..0.or.1.` ~.,data=Variables.strd, nvmax=12, method ="forward")
summary(regfit.fwd.strd)

par(mfrow=c(2,2))
plot(summary(regfit.fwd.strd)$rss ,xlab="Number of Variables ",ylab="RSS",type="l")
plot(summary(regfit.fwd.strd)$adjr2 ,xlab="Number of Variables ",ylab="Adjusted RSq",type="l")

which.max(summary(regfit.fwd.strd)$adjr2)
points(9,summary(regfit.fwd.strd)$adjr2[9], col="red",cex=2,pch=20)

plot(summary(regfit.fwd.strd)$cp, xlab="Number of Variables ", ylab="Cp", type="l")
which.min(summary(regfit.fwd.strd)$cp)
points(9,summary(regfit.fwd.strd)$cp [9],col="red",cex=2,pch=20)

plot(summary(regfit.fwd.strd)$bic ,xlab="Number of Variables ",ylab="BIC",type="l")
which.min(summary(regfit.fwd.strd)$bic )
points(8,summary(regfit.fwd.strd)$bic [8],col="red",cex=2,pch=20)

?plot.regsubsets

plot(regfit.fwd.strd,scale="r2")
plot(regfit.fwd.strd,scale="adjr2") 
plot(regfit.fwd.strd,scale="Cp")
plot(regfit.fwd.strd,scale="bic")

coef(regfit.fwd.strd, 8) #8 perché secondo il BIC l'ottavo subset è il migliore

Variables.fwd.strd<- subset(Variables.strd, select = c(GCM.perennial.surface.water.ice..0.or.1.,Temperature..K.,Pression..Pa.,Monthly.mean.surface.H2O.layer..kg.m2.,Water.ice.column..kg.m2.,Water.ice.mixing.ratio..mol.mol.,Water.ice.effective.radius..m.,H.column..kg.m2.,H2.column..kg.m2.))

# Logistic regression model on standardized data

glm.fit.fwd.strd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd.strd, family = binomial())
summary(glm.fit.fwd.strd)
coef(glm.fit.fwd.strd)

# WE ALWAYS HAVE THE SAME WARNING:
# Warning message:
#   glm.fit: fitted probabilities numerically 0 or 1 occurred 

round(predict(glm.fit.fwd, data=Variables, type = "response"),3)
table(Variables.fwd$GCM.perennial.surface.water.ice..0.or.1.)

#oversampling of class 1
Variables_balanced <- ovun.sample(GCM.perennial.surface.water.ice..0.or.1. ~ ., data = Variables, method = "over")$data #, N = N_desiderato per specificare N
round(predict(glm.fit, data = Variables_balanced, type = "response"),3)
table(Variables_balanced$GCM.perennial.surface.water.ice..0.or.1.)

glm.fit.bal <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables_balanced, family = binomial())
summary(glm.fit.bal)
coef(glm.fit.bal)


# Multicollinearity management
# Variables with indexes > 10 are multicollinear 
vif_model <- vif(glm.fit.fwd)
vif_model

# Correlation matrix creation
Variables$GCM.perennial.surface.water.ice..0.or.1.<-as.numeric(Variables$GCM.perennial.surface.water.ice..0.or.1.)
correlation_matrix <- cor(Variables.strd)
correlation_matrix

heatmaply(correlation_matrix, 
          col = colorRampPalette(c("blue", "white", "red"))(100),
          main = "Correlation matrix",
          fontsize = 10)

#Try to merge Pression..Pa. and Water.vapor.column..kg.m2. (index=0.78372909) into Press_WaVaCo_Index
# We use arithmetic mean on standardized data
Variables.strd$Press_WaVaCo_Index <- rowMeans(Variables.strd[, c("Pression..Pa.", "Water.vapor.column..kg.m2.")])

glm.fit_Press_WaVaCo <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Press_WaVaCo_Index +Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.strd, family = binomial())
summary(glm.fit_Press_WaVaCo)
coef(glm.fit_Press_WaVaCo)
#Didn't work
vif_model_Press_WaVaCo<-vif(glm.fit_Press_WaVaCo)
vif_model_Press_WaVaCo
#but we reduced multicollinearity 
correlation_matrix2 <- cor(Variables.strd)
correlation_matrix2

#Try to merge H2.column..kg.m2. and Electron.number.density..m.3. (index 0.50452656) into H2_EleNumD_Index
Variables.strd$H2_EleNumD_Index <- rowMeans(Variables.strd[, c("H2.column..kg.m2.", "Electron.number.density..m.3.")])

glm.fit_Press_WaVaCo_H2_EleNumD <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Press_WaVaCo_Index +Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2_EleNumD_Index, data = Variables.strd, family = binomial())
#Didn't work
vif_model_Press_WaVaCo_H2_EleNumD<-vif(glm.fit_Press_WaVaCo_H2_EleNumD)
vif_model_Press_WaVaCo_H2_EleNumD
#but we reduced multicollinearity 
correlation_matrix3 <- cor(Variables.strd)
correlation_matrix3

#to reset the Variables.strd data
Variables.strd<-subset(Variables.strd, select=-Press_WaVaCo_Index)
Variables.strd<-subset(Variables.strd, select=-H2_EleNumD_Index)

#-------------------------------------------------------------------#
#                        Logistic regression                        #
#                            second part                            #
#-------------------------------------------------------------------#

# Definition of estimated probabilities
summary(glm.fit.fwd)$coef
glm.probs<-predict(glm.fit.fwd,type="response")
glm.probs[1:10]

Variables$GCM.perennial.surface.water.ice..0.or.1.<-as.factor(Variables$GCM.perennial.surface.water.ice..0.or.1.)
contrasts(Variables$GCM.perennial.surface.water.ice..0.or.1.)

# Prediction Binarization
glm.pred<-rep(0,3072)
glm.pred[glm.probs >.5]=1
glm.pred

# Preliminar confusion matrix
table(glm.pred,Variables$GCM.perennial.surface.water.ice..0.or.1.)
mean(glm.pred==Variables$GCM.perennial.surface.water.ice..0.or.1.)
mean(glm.pred!=Variables$GCM.perennial.surface.water.ice..0.or.1.)


#-------------------------------------------------------------------#
#                        Cross-Validation                           #
#-------------------------------------------------------------------#

# 1st k-fold cross validation. 10 folds
set.seed(17)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 
#Errors indicate the difference between the model's predictions and the actual data during the cross-validation process;
# these appear to be relatively low, as they are all below 0.04

# 2nd k-fold cross validation. 10 folds
set.seed(1)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

# 3rd k-fold cross validation. 10 folds
set.seed(2)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

# 4th k-fold cross validation. 10 folds
set.seed(3)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

# 5th k-fold cross validation. 10 folds
set.seed(4)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

# 6th k-fold cross validation. 10 folds
set.seed(5)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

# 7th k-fold cross validation. 10 folds
set.seed(6)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

# 8th k-fold cross validation. 10 folds
set.seed(7)
cv.error.10=rep(0,10)
for (i in 1:10){
  glm.fit.fwd <- glm(`GCM.perennial.surface.water.ice..0.or.1.` ~ Temperature..K.+Pression..Pa.+Monthly.mean.surface.H2O.layer..kg.m2.+Water.ice.column..kg.m2.+Water.ice.mixing.ratio..mol.mol.+Water.ice.effective.radius..m.+H.column..kg.m2.+H2.column..kg.m2., data = Variables.fwd) #, family = binomial()
  cv.error.10[i]=cv.glm(Variables.fwd,glm.fit.fwd,K=10)$delta[1]
}
cv.error.10 

#plot single repetition
set.seed(17)
degrees <- 1:10 
for (deg in degrees) {
  formula <- reformulate(
    termlabels = c(paste("poly(Temperature..K.,", deg, ")"), 
                   "Pression..Pa.", 
                   "Monthly.mean.surface.H2O.layer..kg.m2.",
                   "Water.ice.column..kg.m2.", 
                   "Water.ice.mixing.ratio..mol.mol.", 
                   "Water.ice.effective.radius..m.", 
                   "H.column..kg.m2.", 
                   "H2.column..kg.m2."),
    response = "`GCM.perennial.surface.water.ice..0.or.1.`"
  )

  glm.fit.fwd <- glm(formula, data = Variables.fwd)

  cv.error.10[deg] <- cv.glm(Variables.fwd, glm.fit.fwd, K=10)$delta[1]
}
# Grafico dei risultati
plot(degrees, cv.error.10, type="b", xlab="Degree of Polynomial", ylab="Mean Squared Error", main="10-fold CV")



# plot all repetitions
cv_errors <- matrix(0, nrow = 10, ncol = 8)

seeds <- c(1,2,3,4,5,6,7,17)

for (j in 1:length(seeds)) {
  set.seed(seeds[j])
  cv.error.10 <- rep(0, 10)
  for (i in 1:10){
    glm.fit.fwd <- glm(formula = `GCM.perennial.surface.water.ice..0.or.1.` ~ ., 
                       data = Variables.fwd, 
                       family = binomial()) 
    cv.error.10[i] <- cv.glm(Variables.fwd, glm.fit.fwd, K=10)$delta[1]
  }
  cv_errors[, j] <- cv.error.10
}

matplot(cv_errors, type = "l", lty = 1, xlab = "Fold Number", ylab = "Cross Validation Error",
        main = "10-fold Cross Validation repetition")

legend("topright", legend = paste("Seed", seeds), col = 1:length(seeds), lty = 1)


#-------------------------------------------------------------------#
#                   Confusion matrix adjustment                     #
#-------------------------------------------------------------------#

# Threshold choise
thresholds <- seq(0, 1, by = 0.001)
false_positives <- numeric(length(thresholds))
false_negatives <- numeric(length(thresholds))
true_positives <- numeric(length(thresholds))
true_negatives <- numeric(length(thresholds))

# Actual values from the dataset
actuals <- Variables.fwd$GCM.perennial.surface.water.ice..0.or.1.

for (i in seq_along(thresholds)) {
  threshold <- thresholds[i]
  predictions <- ifelse(glm.probs >= threshold, 1, 0)
  
  # Ensure both factors are present
  confusion_matrix <- table(factor(predictions, levels = c(0, 1)), factor(actuals, levels = c(0, 1)))
  
  if (length(confusion_matrix) == 4) {
    true_positives[i] <- confusion_matrix[2, 2]
    true_negatives[i] <- confusion_matrix[1, 1]
    false_positives[i] <- confusion_matrix[2, 1]
    false_negatives[i] <- confusion_matrix[1, 2]
  } else {
    true_positives[i] <- 0
    true_negatives[i] <- 0
    false_positives[i] <- 0
    false_negatives[i] <- 0
  }
}

false_positive_rate <- false_positives / (false_positives + true_negatives)
false_negative_rate <- false_negatives / (false_negatives + true_positives)

error_data <- data.frame(Threshold = thresholds, FalsePositiveRate = false_positive_rate, FalseNegativeRate = false_negative_rate)

ggplot(error_data, aes(x = Threshold)) +
  geom_line(aes(y = FalsePositiveRate, color = "False Positive")) +
  geom_line(aes(y = FalseNegativeRate, color = "False Negative")) +
  labs(title = "Error Rates by Threshold", x = "Threshold", y = "Rate") +
  theme_minimal() +
  scale_color_manual(values = c("False Positive" = "blue", "False Negative" = "red")) +
  theme(legend.title = element_blank())

thresholds[which.min(abs(false_positive_rate - false_negative_rate))] #we chose this threshold

# Prediction Binarization with the chosen Threshold
glm.pred<-rep(0,3072)
glm.pred[glm.probs >.074]=1
glm.pred

#Confusion matrix
table(glm.pred,Variables$GCM.perennial.surface.water.ice..0.or.1.)
mean(glm.pred==Variables$GCM.perennial.surface.water.ice..0.or.1.)
mean(glm.pred!=Variables$GCM.perennial.surface.water.ice..0.or.1.)

#-------------------------------------------------------------------#
#                            ROC curve                              #
#-------------------------------------------------------------------#


truth<- Variables$GCM.perennial.surface.water.ice..0.or.1.
predictions <- glm.probs  

pred <- prediction(predictions, truth)

perf <- performance(pred, "tpr", "fpr")

# Creazione del plot
plot(perf, col = "blue", main = "ROC Curve")  
abline(a = 0, b = 1, lty = 3, col="lightgrey")


#-------------------------------------------------------------------#
#                     Probabilites dataframe                        #
#-------------------------------------------------------------------#

Probs_df <- matrix(as.numeric(glm.probs), nrow = 64, ncol = 48, byrow = TRUE)

Probs_df <- t(Probs_df)

Probs_df <- as.data.frame(Probs_df)

row_headers<- c(-180,-174,-168,-162,-157,-151,-145,-140,-134,-128,-122,-117,-111,-105,-100,-94,-88,-82,-77,-71,-65,-60,-54,-48,-42,-37,-31,-25,-20,-14,-8,-2,+2,+8,+14,+20,+25,+31,+37,+42,+48,+54,+60,+65,+71,+77,+82,+88,+94,+100,+105,+111,+117,+122,+128,+134,+140,+145,+151,+157,+162,+168,+174,+180)
colnames(Probs_df)<-row_headers

column_headers<- c(-9.00000e+01,-8.61702e+01,-8.23404e+01,-7.85106e+01,-7.46809e+01,-7.08511e+01,-6.70213e+01,-6.31915e+01,-5.93617e+01,-5.55319e+01,-5.17021e+01,-4.78723e+01,-4.40426e+01,-4.02128e+01,-3.63830e+01,-3.25532e+01,-2.87234e+01,-2.48936e+01,-2.10638e+01,-1.72340e+01,-1.34043e+01,-9.57447e+00,-5.74468e+00,-1.91489e+00,1.91489e+00,5.74468e+00,9.57447e+00,1.34043e+01,1.72340e+01,2.10638e+01,2.48936e+01,2.87234e+01,3.25532e+01,3.63830e+01,4.02128e+01,4.40426e+01,4.78723e+01,5.17021e+01,5.55319e+01,5.93617e+01,6.31915e+01,6.70213e+01,7.08511e+01,7.46809e+01,7.85106e+01,8.23404e+01,8.61702e+01,9.00000e+01)
rownames(Probs_df)<-column_headers

heatmaply(Probs_df, 
          col = colorRampPalette(c( "lightyellow","orange", "red"))(100),
          main = "Probabilities",
          fontsize = 10,
          Rowv = FALSE,
          Colv = FALSE,
          xlab="Longitude",
          ylab="Latitude")

