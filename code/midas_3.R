# 
rm(list=ls())
# setwd("/Users/carlo/Dropbox/R/Midas/data")
setwd("/home/giuliano/Dropbox/R/Midas/data")
load("midas_data.RData")
load("functions.RData")


library(Matrix)
library(ggplot2)
library(MASS)



z<-40
arg<-c("climat",8,2,"umidas","month","spread",251,4,"polynomial","day","dow_jones",251,4,"polynomial","day")
h<-60

# # data[,1]<-format(data[,1],"%d/%m/%Y")
# 
# # colnames(prov)<-c("Date","pib",paste("coco",c(1:12),sep="",collapse=NULL))
# # reg<-summary(lm(pib~coco1+coco2+coco3+coco4+coco5+coco6+coco7+coco8+coco9+coco10+coco11+coco12,prov))
# # prov2<-as.matrix(prov)
# # reg<-summary(lm(prov2[,2]~prov2[,3:ncol(prov)]))
# # 
# # "coco1":"coco3"
# # 
# 
# # init<-c(1,1,-0.1,1,-0.1,1,-0.1)
# # h<--60
# # z<-30
# 
# 
# 
# # for (i in 150:220){
# # arg<-c("climat",8,2,"umidas","month","spread",250,4,"polynomial","day","cac_40",250,4,"polynomial","day")
# # fin<-loop_rmsfe(60,30,data,arg,init)
# # print(i)
# # print(fin[[2]])
# # # fin[[3]]
# # }
# 
# 
# 
# # loop_rmsfe<-function(h,z,data,arg,init){  
#   #data : data set that contains the data we are working on
#   #arg : arguments vector. For each regressor, it contains the its name, the number of days, the degree of the constaints polynom and the method (polynomial, exponential)
#   #h : number of days of lag between the end of the estimation period and the end of the quarter we want to estimate.
#   #if h=60 for example, we will forecast the values 60 days before their release
#   #init : vector of initial values for minimization (the exponential case)
#   #This function is split in two parts: 
#   # - the first part aims at computing a regression
#   # - the second part aims at computing the daily coefficients, thanks to the regression coefficients
#   
#   
#   
#   #I put the argument vector ("arg") into a matrix
#   arguments<-matrix(arg,ncol=5,byrow=TRUE)
#   #If the method for one variable is umidas, then I replace the degree by the number of lag days
#   arguments[arguments[,4]=="umidas",3]<-(as.numeric(arguments[arguments[,4]=="umidas",2])-1)
#   
#   variables<-c(arguments[,1])
#   retards<-as.numeric(c(arguments[,2]))
#   degres<-as.numeric(c(arguments[,3]))
#   method<-c(arguments[,4])
#   hf<-c(arguments[,5])
#   nb_var<-nrow(arguments)
#   
#   
#   ############################
#   #I compute the regression###
#   ############################
#   #I create a table that contains the data for the regression
#   #that's why I call the "base" function
#   
#   #I regress.
#   
#   forecast<-c()
#   upper_bound<-c()
#   lower_bound<-c()
#   cusum<-c()
#   cusum_up<-c()
#   cusum_down<-c()
#   W<-c()
#   prov2<-base_reg(data,arg,h)
#   prov2<-data.frame(prov2[,-1])
#   prov3<-prov2
#   for (j in 1:(z-1)){  
#     nb_row<-(nrow(prov3)+j-z)
#     nb_row2<-nb_row+1
#     if (all(method==rep("exponential",nb_var))){
#       
#       
#       prov2<-prov3[1:nb_row,]
#       
#       #error message if the number of starting values and the number of parameters to estimate don't match
#       if (length(init)>(sum(degres)+1)){
#         print("There are more initial values than estimated parameters")
#       }
#       if (length(init)<(sum(degres)+1)){
#         print("There are less initial values than estimated parameters")
#       }
#       
#       #if the method specified is non-linear, I first write the formula of the function to be minimized
#       expr1<-paste("+",apply(arguments,1,expression),collapse="",sep="")
#       expr1<-paste("const",expr1,sep="")
#       #and then I compute a non linear minimization
#       reg<-nlm(fc,init,expr=expr1,variables=variables,deg=degres,data=prov)
#       #I stock the estimates in a vector with a standard name
#       estimated_coeffs<-reg$estimate  
#       
#     }else{
#       
#       #else (ie if the method is polynomial or umidas), I compute a linear model
#       reg<-summary(lm(pib~.,data=prov2))
#       #I stock the estimates in a vector with a standard name
#       estimated_coeffs<-reg$coefficient[,1]
#       sigma2<-reg$residuals%*%reg$residuals/(nrow(prov2)-ncol(prov2)+1)
#     }
#     
#     ###################################
#     #I compute the daily coefficients##
#     ###################################
#     #I create a matrix whose purpose is to receive the values of the daily coefficients
#     coefficients<-matrix(NA,nrow=max(retards),ncol=(length(retards)+1))
#     #the first column of this matrix contains the day id.
#     coefficients[,1]<-c(1:max(retards))
#     
#     #exponential case
#     if (all(method==rep("exponential",nb_var))){
#       #I first have to assign the estimated values to the matching coefficients
#       #that's why I compute a loop that gives their values to each subset of "estimated_coeffs"
#       #I first assign a value to the intercept, named "const"
#       const<-reg$estimate[1]  
#       #I assign an initial value to a variable, very useful for this loop
#       s<-2
#       for (i in 1:length(variables)){
#         #I select a subset of the estimated vector, which matches to the corresponding coefficients
#         values<-estimated_coeffs[s:(degres[i]+1)]
#         #I increase the value of s, in order to select the proprer subset of estimates during the next iteration
#         s<-(degres[i]+2)
#         #I compute the daily coefficients by computing a matrix product and then by dividing each coefficient by the sum of coefficients
#         M<-Mat(retards[i],degres[i])
#         prov2<-sapply(M[,-1]%*%values,exp)/sum(sapply(M[,-1]%*%values,exp))
#         #I put the daily coefficients (prov2) into a matrix containing the whole set of coefficients. 
#         #the first column of the matrix contains the day id
#         coefficients[,i+1]<-c(prov2,rep(NA,(nrow(coefficients)-length(prov2))))
#       }  
#     }else{
#       #polynomial case
#       const<-estimated_coeffs[1]
#       #n is a useful variable for the following loop. I have to initialise it.
#       n<-1
#       for (i in 1:length(variables)){
#         #At each iteration I need to select the m-th to the n-th coefficients stocked into reg$coefficients
#         #the value of m is 1 plus the value of the last coefficient selected during the previous iteration
#         m<-(n+1)
#         #the value of n is m plus the number of coefficients i need to select, ie degres[i]
#         n<-(m+degres[i])
#         #I prepare the transformation matrix
#         M<-Mat(retards[i],degres[i])
#         #I put the computed daily coefficients into one column of the "coefficients" matrix.
#         if (method[i]=="polynomial"){
#           #if the constraint is polynomial, I have to compute a matrix product with the estimated coefficients and the polynomial matrix
#           coefficients[,(i+1)]<-c(M%*%reg$coefficients[m:n,1],rep(NA,(length(coefficients[,(i+1)])-length(M%*%reg$coefficients[m:n,1]))))
#         }
#         if (method[i]=="umidas"){
#           #if the model is unconstrained, I just put the coefficents into the matrix, without any matrix product.
#           coefficients[,(i+1)]<-c(reg$coefficients[m:n,1],rep(NA,(length(coefficients[,(i+1)])-length(reg$coefficients[m:n,1]))))
#         }
#       }
#     }
#     
# #     arg2<-c("climat",8,2,"umidas","month","spread",251,4,"umidas","day","dow_jones",251,4,"umidas","day")
# #     test<-base_reg(data,arg2,h)
# #     test<-test[,-1]
# #     forecast2[j]<-const
# #     f<-2
# #     for (k in 1:nb_var){
# #       g<-retards[k]+f-1
# #       #foo is the produt of the coefficients and the regressors related to one variable
# #       foo<-as.matrix(test[nb_row2,f:g])%*%as.matrix(na.omit(coefficients[,(k+1)]))
# # #       print(na.omit(coefficients[,(k+1)]))
# # #       print(test[nb_row2,f:g])
# #       forecast2[j]<-forecast2[j]+foo
# #       f<-g+1
# #     }
#        forecast[j]<-reg$coefficients[1,1]+as.matrix(prov3[nb_row2,2:ncol(prov3)])%*%(as.matrix(reg$coefficients[2:nrow(reg$coefficients),1]))
#     
#     #block that computes the confidence intervals of the forecasting values
#       X<-as.matrix(prov2[,-1])
#       #forecasting error variance
#       fev<-sqrt(sigma2*(1+as.matrix(prov3[nb_row2,2:ncol(prov3)])%*%ginv(t(X)%*%X)%*%t(as.matrix(prov3[nb_row2,2:ncol(prov3)]))))
#       kk<-qt(0.90,nrow(X)-ncol(X))*fev
#       upper_bound[j]<-forecast[j]+kk
#       lower_bound[j]<-forecast[j]-kk
#     
#     #block thaht computes the cusum test
#       #computation of the recursive residual
#       a<-1.143
#       k<-ncol(X)
#       cusum[j]<-(prov3[nb_row2,1]-forecast[j])/sqrt(1+as.matrix(prov3[nb_row2,2:ncol(prov3)])%*%ginv(t(X)%*%X)%*%t(as.matrix(prov3[nb_row2,2:ncol(prov3)])))
#       cusum_up[j]<-a*sqrt(nrow(prov3)-ncol(X))+2*a*(nb_row2-k)/sqrt(nrow(prov3)-ncol(X))
#       cusum_down[j]<-a*sqrt(nrow(prov3)-ncol(X))-2*a*(nb_row2-k)/sqrt(nrow(prov3)-ncol(X))
#       print(cusum_down[i])
#       W[j]<-sum(cusum)/sqrt(sigma2)
#   }
# 
#   abcisses<-1:length(forecast)
#   realite<-c(prov3[(nrow(prov3)+1-length(forecast)):nrow(prov3),1])
#   graph_rmsfe<-data.frame(abcisses,realite,forecast,upper_bound,lower_bound,cusum)
#   rmsfe<-sqrt(mean((graph_rmsfe[,3]-graph_rmsfe[,2])^2))
#   
#   graph_cusum<-data.frame(abcisses,cusum,cusum_up,cusum_down)
# 
# 
#   courbe<-ggplot()
#   courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=realite,colour="realite"))
# #   courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=lower_bound,colour="bh"))
# #   courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=cusum,colour="cusum"))  
#   courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=forecast,colour="prevision"))
# courbe  
# 
# 
# courbe_cusum<-ggplot()
# courbe_cusum<-courbe_cusum + geom_path(data=graph_cusum,aes(x=abcisses,y=cusum,colour="cusum"))
# courbe_cusum<-courbe_cusum + geom_path(data=graph_cusum,aes(x=abcisses,y=cusum_up,colour="cusum_up"))
# courbe_cusum<-courbe_cusum + geom_path(data=graph_cusum,aes(x=abcisses,y=cusum_down,colour="cusum_down"))  
# # courbe_cusum<-courbe_cusum + geom_path(data=graph_cusum,aes(x=abcisses,y=W,colour="W"))
# courbe_cusum  
# 
# return(list(graph_rmsfe,rmsfe,courbe,sigma2))
# }

start<-Sys.time()
z<-30
arg<-c("climat",8,2,"umidas","month","spread",251,4,"polynomial","day","dow_jones",251,4,"polynomial","day")
h<-60


fin<-loop_rmsfe(60,50,data,arg,init)
end<-Sys.time()

end 
start
end-start
# fin[[1]]
fin[[2]]
fin[[3]]
test<-base_reg(data,arg,87)
coco2<-regression_midas(data,arg,-30)
coco2