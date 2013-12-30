rm(list=ls())
setwd("/Users/carlo/Dropbox/R/Midas/donnees")

library(Matrix)
library(ggplot2)

Mat<-function(a,deg){
  #mat : matrice vide
  #a : nombre de jours de retards du regresseurs pris en compte ; cela correspond au nombre de lignes de mat
  mat<-matrix(NA,nrow=a,ncol=(deg+1))
  t<-deg
  for (i in 0:t){
    mat[,i+1]<-c(0:(a-1))^i
  }
  return(mat)
}



################################
#Diff?renciateur de variables###
################################
di<-function(x)
{
  return(c(x[2:length(x)]-x[1:(length(x)-1)]))
}



#function that finds the first row that does not contain any NA
first_row<-function(data){
  #this function finds the first non NA row of every variable and then returns the maximum of these values
  first_rows<-c()
  for (j in 3:ncol(data)){
    foo<-which(!is.na(data[,j]))[1]  
    first_rows<-c(first_rows,foo)
  }
  res<-max(first_rows)
}


#this function finds the first (most remote) observation of the regressand that can be put into the midas table. 
#for example : if I want to regress on the latests 60 days of one regressor, with a forecasting delay (h) of 10 days, 
#and if there are only 30 observations before the first non NA observation of the regressand, I will have to delete these observation.
find_first_row<-function(data,variables,retards,j,h){
  
  #I look for the date of the first non NA value of the j-st variable + the number of lags of this variable I want to consider + the forecasting delay (h)
  table_prov<-data[!is.na(data[,variables[j]]),c("Date",variables[j])]
  scope_first_obs<-table_prov[retards[j],"Date"]
  scope_first_obs<-scope_first_obs+h
  #I look for the number of the first non NA obs of the regressand, whose date is greater than scope_first_obs.
  first_non_NA_dependent<-min(which(data[!is.na(data[,2]),"Date"]>scope_first_obs))
  return(first_non_NA_dependent)
}




polynomial_add_on<-function(data_set,nb_var,variables,retards,degres,method){
  mat_list<-list()
  for (i in 1:nb_var){
    n_row<-retards[i]  
    if(method[i]=="polynomial"){
      n_col<-degres[i]
      ma<-Mat(retards[i],degres[i])
      mat_list[[i]]<-ma
    }else{
      ma<-diag(1,nrow=n_row,ncol=n_row)
      mat_list[[i]]<-ma
    }
  }
  return(bdiag(mat_list))
}



base_reg<-function(data,arg,h){
  
  arguments<-matrix(arg,ncol=5,byrow=TRUE)
  variables<-c(arguments[,1])
  retards<-as.numeric(c(arguments[,2]))
  degres<-as.numeric(c(arguments[,3]))
  method<-c(arguments[,4])
  hf<-c(arguments[,5])
  nb_var<-nrow(arguments)
  # print(arguments)
  
  #first I delete the first rows which are not filled with values, for some variable.
  fr<-first_row(data)
  data<-data[-c(1:fr),]
  #I need the first column of the "data" dataset to be exhaustive, ie to have every day between the first and the last date.
  full_dates<-data.frame(Date=seq(data[1,1], data[nrow(data),1], by="days"))
  data<-merge(data,full_dates,all.y=TRUE)
  
  #I look for the indicia of the rows of the regressand that are not NA.
  indicia<-which(!is.na(data[,2]))  
  #I create a table with the dates and the non-NA regressand.
  q<-data.frame(data[indicia,c(1,2)])
  #I first have to find the date of the first forecastable value 
  #given the number of lags and the forecasting delay (h) I have chosen
  start<-max(sapply(1:length(variables),find_first_row,data=data,variables=variables,retards=retards,h=h))
  #the the regressands whose values arer < start are useless. So I delete them
  q<-q[-c(1:start),]
  
  #I take in consideration the forecasting delay (h)
  #d'abord, je d?cale vers le haut (je recule) toutes les dates d'un nombre de jours correspondant au jours mis en h.
  delayed_dates<-data[indicia,1]-h
  #now I have to find the number of the row of each date in the vector (delayed_dates)
  indicia2<-(sapply(delayed_dates,function(x)which(data[,1]==x)))
  
  #I create a table that will contain the regressors for the midas regression
  stack<-as.data.frame(matrix(NA,ncol=sum(retards)))
  
  for(i in start:length(na.omit(data[,2]))){  
    #prov corresponds to one row of the midas table
    #the number of columns of prov is equal to the sum of the lags of each regressor
    prov<-c()
    
    for (j in 1:nb_var){
      #for the j-st vaiable, I first select every value between the first one and the most close to the forecasted value (whose id is given by indicia2)
      foo1<-t(na.omit(data[1:indicia2[[i]],variables[j]]))
      lg<-length(foo1)
      
      #then I only keep the good number of lags, ie the "retards[j]" last values
      foo1<-rev(foo1[-c(1:(max(lg-retards[j],0)))])
      
      #then I paste these hash into prov
      prov<-c(prov,foo1)
    }     
    
    #I stack the rows (prov)
    stack[(i-start),]<-prov   
    
  }
  #I have to compute a matrix transformation if method = polynomial
  #that's why I compute a matricial product, with a block diagonal matrix.
  #each block correspond to the set of lags of a given varibale.
  #if the method for this variable is not "polynomial", the corresponding block is an identity matrix.
  matrix_multiplicator<-polynomial_add_on(data_set,nb_var,variables,retards,degres,method)
  midas_regressors<-as.matrix(stack)%*%matrix_multiplicator
  midas_regressors<-as.matrix(midas_regressors)
  
  #I create the midas data set by binding the regressor and the regressand.  
  #now, for each regressands observation I have the matching regressors.
  midas_table<-as.data.frame(cbind(q,midas_regressors))
  
  #I change the names of the columns of the table
  column_names<-col_names(base_midas,variables,retards,degres,method,nb_var)
  colnames(midas_table)<-column_names
  
  return(midas_table)
  
}

# column_names<-col_names(base_midas,variables,retards,degres,method,nb_var)

col_names<-function(table,variables,retards,degres,method,nb_var){
  a<-3
  nb_col<-c()
  names<-c("Date","pib")
  for (i in 1:nb_var){
    if(method[i]=="polynomial"){
      nb_col[i]<-degres[i]+1
    }else{
      nb_col[i]<-retards[i]
    }
    
    b<-a+nb_col[i]-1
    
    names[a:b]<-c(paste(variables[i],c(1:nb_col[i]),sep="",collapse=NULL))
    a<-a+nb_col[i]
    
  }
  return(names)  
}





regression_midas<-function(data,arg,h,init){
  #data : data set that contains the data we are working on
  #arg : arguments vector. For each regressor, it contains the its name, the number of days, the degree of the constaints polynom and the method (polynomial, exponential)
  #h : number of days of lag between the end of the estimation period and the end of the quarter we want to estimate.
  #if h=60 for example, we will forecast the values 60 days before their release
  #init : vector of initial values for minimization (the exponential case)
  #This function is split in two parts: 
  # - the first part aims at computing a regression
  # - the second part aims at computing the daily coefficients, thanks to the regression coefficients
  
  #   arg<-c("cac_40",250,2,"exponential","change",200,2,"exponential","dow_jones",150,2,"exponential")
  #init<-c(c(1),c(1,-0.01),c(1,-0.01),c(1,-0.01))
  #h<-60
  #I put the argument vector ("arg") into a matrix
  arguments<-matrix(arg,ncol=5,byrow=TRUE)
  #If the method for one variable is umidas, then I replace the degree by the number of lag days
  arguments[arguments[,4]=="umidas",3]<-(as.numeric(arguments[arguments[,4]=="umidas",2])-1)
  
  variables<-c(arguments[,1])
  retards<-as.numeric(c(arguments[,2]))
  degres<-as.numeric(c(arguments[,3]))
  method<-c(arguments[,4])
  hf<-c(arguments[,5])
  nb_var<-nrow(arguments)
  
  ############################
  #I compute the regression###
  ############################
  #I create a table that contains the data for the regression
  #that's why I call the "base" function
  prov<-base_reg(data,arg,h)
  prov<-data.frame(prov[,-1])
  
  #I regress.
  if (all(method==rep("exponential",nb_var))){
    prov<-base_reg(data,arg,h)
    prov<-data.frame(prov[,-1])
    
    #error message if the number of starting values and the number of parameters to estimate don't match
    if (length(init)>(sum(degres)+1)){
      print("There are more initial values than estimated parameters")
    }
    if (length(init)<(sum(degres)+1)){
      print("There are less initial values than estimated parameters")
    }
    
    #if the method specified is non-linear, I first write the formula of the function to be minimized
    expr1<-paste("+",apply(arguments,1,expression),collapse="",sep="")
    expr1<-paste("const",expr1,sep="")
    #and then I compute a non linear minimization
    reg<-nlm(fc,init,expr=expr1,variables=variables,deg=degres,data=prov)
    #I stock the estimates in a vector with a standard name
    estimated_coeffs<-reg$estimate  
    
  }else{
    prov<-base_reg(data,arg,h)
    prov<-data.frame(prov[,-1])
    #else (ie if the method is polynomial or umidas), I compute a linear model
    reg<-summary(lm(pib~.,data=prov))
    #I stock the estimates in a vector with a standard name
    estimated_coeffs<-reg$coefficient[,1]
    
  }
  ###################################
  #I compute the daily coefficients##
  ###################################
  #I create a matrix whose purpose is to receive the values of the daily coefficients
  coefficients<-matrix(NA,nrow=max(retards),ncol=(length(retards)+1))
  #the first column of this matrix contains the day id.
  coefficients[,1]<-c(1:max(retards))
  
  #exponential case
  if (all(method==rep("exponential",nb_var))){
    #I first have to assign the estimated values to the matching coefficients
    #that's why I compute a loop that gives their values to each subset of "estimated_coeffs"
    #I first assign a value to the intercept, named "const"
    const<-reg$estimate[1]  
    #I assign an initial value to a variable, very useful for this loop
    s<-2
    for (i in 1:length(variables)){
      #I select a subset of the estimated vector, which matches to the corresponding coefficients
      values<-estimated_coeffs[s:(degres[i]+1)]
      #I increase the value of s, in order to select the proprer subset of estimates during the next iteration
      s<-(degres[i]+2)
      #I compute the daily coefficients by computing a matrix product and then by dividing each coefficient by the sum of coefficients
      M<-Mat(retards[i],degres[i])
      prov2<-sapply(M[,-1]%*%values,exp)/sum(sapply(M[,-1]%*%values,exp))
      #I put the daily coefficients (prov2) into a matrix containing the whole set of coefficients. 
      #the first column of the matrix contains the day id
      coefficients[,i+1]<-c(prov2,rep(NA,(nrow(coefficients)-length(prov2))))
    }  
  }else{
    #polynomial case
    const<-estimated_coeffs[1]
    #n is a useful variable for the following loop. I have to initialise it.
    n<-1
    for (i in 1:length(variables)){
      #At each iteration I need to select the m-th to the n-th coefficients stocked into reg$coefficients
      #the value of m is 1 plus the value of the last coefficient selected during the previous iteration
      m<-(n+1)
      #the value of n is m plus the number of coefficients i need to select, ie degres[i]
      n<-(m+degres[i])
      #I prepare the transformation matrix
      M<-Mat(retards[i],degres[i])
      #I put the computed daily coefficients into one column of the "coefficients" matrix.
      if (method[i]=="polynomial"){
        #if the constraint is polynomial, I have to compute a matrix product with the estimated coefficients and the polynomial matrix
        coefficients[,(i+1)]<-c(M%*%reg$coefficients[m:n,1],rep(NA,(length(coefficients[,(i+1)])-length(M%*%reg$coefficients[m:n,1]))))
      }
      if (method[i]=="umidas"){
        #if the model is unconstrained, I just put the coefficents into the matrix, without any matrix product.
        coefficients[,(i+1)]<-c(reg$coefficients[m:n,1],rep(NA,(length(coefficients[,(i+1)])-length(reg$coefficients[m:n,1]))))
      }
    }
  }
  
  return(coefficients)
}







##############################################################################
#fonctions nécessaires à l'écriture du programme de minimisation non linéaire#
##############################################################################

#fonction qui écrit un monôme de degré "i" associé au jème regresseur
monom<-function(i,j,variable){
  #i : degré du monôme
  #j : numéro du regresseur auquel sera associé le coefficient contenant le monôme en question
  #co : vecteur contenant les coefficients à estimer
  res<-paste("+co_",variable,"[",i,"]*",j,"^",i,sep="")
  return(res)
}


#fonction qui concatène les monômes de degré < deg, relatifs au k-ème coefficient, 
#et qui insère dans une exponentielle
coefficient<-function(k,deg,variable){
  #k : numéro du regresseur auquel est associé le coefficient en question
  #deg: degré du polynôme situé dans l'exponentielle
  
  res<-paste(sapply(c(1:deg),monom,j=k,variable=variable)[1:deg],collapse="",sep="")
  #j'élimine le premier signe "+"
  res<-substring(res,2)
  #je passe à l'exponentielle
  res<-paste("exp(",res,")",sep="")
  return(res)
}

#fonction qui concatème les monômes de degré < deg, relatifs au k-ème coefficient, 
#et qui insère dans une exponentielle
#et qui multiplie par le regresseur idoine
coefficient_prod<-function(k,deg,variable){
  #k : numéro du regresseur auquel est associé le coefficient en question
  #deg: degré du polynôme situé dans l'exponentielle 
  
  res<-paste(sapply(c(1:deg),monom,j=k,variable=variable)[1:deg],collapse="",sep="")
  #j'élimine le premier signe "+"
  res<-substring(res,2)
  #je passe à l'exponentielle et je multiplie par le k-ième regresseur
  res<-paste("exp(",res,")*","prov[,'",paste(variable,k,sep=""),"']",sep="")
  return(res)
}


#fonction qui écrit l'expression pour la minimisation non linéaire
expression<-function(x){
  #x: argument vector made of three components : 
  #a : nombre de coefficients (ie nombre de regresseurs)
  #deg : degré du polynôme situé dans l'exponentielle (ie nombre de coefficients à estimer)
  #variable: name of the variable on which it will be regressed
  
  #I split the argument vectors into more explicit variables
  variable<-x[1]
  a<-x[2]
  deg<-x[3]
  
  #j'écris un vecteur contenant les coefficients non normalisés multipliés par leur regresseur 
  coeff_prod<-paste(sapply(c(1:a),coefficient_prod,deg=deg,variable=variable),sep="")
  
  #je crée le dénominateur : il suffit de sommer les coefficients
  denom<-paste("+",sapply(c(1:a),coefficient,deg=deg,variable=variable),collapse="",sep="")
  #je retire le premier signe "+"
  denom<-substring(denom,2)
  
  #j'ecris l'expression finale :
  #je somme les "norm_coeffs"
  res<-paste("+",coeff_prod,collapse="",sep="")
  res<-substring(res,2)
  #je divise par norm
  res<-paste("(",res,")/(",denom,")",sep="")
  #j'ajoute une constante
  #res<-paste("const+",res,sep="")
  return(res)
}

#fonction à minimiser
#l'argument de cette fonction est un vecteur constitué des paramètres que nous cherchons à estimer
#la longueur de ce vecteur est donc de deg+1 (il ne faut pas oublier qu'il y a aussi une constante à estimer)
fc<-function(x,expr,variables,deg,data){  
  #I assign the values that are given in the argument to the coefficients to be estimated
  #that's why I compute a loop that gives their values to each subset of "coefficients"
  #I first assign a value to the intercept, named "const"
  const<-x[1]
  #I assign an initial value to a variable, very useful for this loop
  s<-2
  for (i in 1:length(variables)){
    #I select a subset of the argument vector, which matches to the corresponding coefficients
    values<-x[s:(deg[i]+1)]
    s<-(deg[i]+2)
    eval(parse(text=paste("co_",variables[i],"<-values",sep="")))
  }
  res<-sum((data[,"pib"]-eval(parse(text=expr)))^2)
  return(res)
}

#fonction qui permet de dessiner la courbe des coefficents estimés 
grapheur <- function(d,vect)
{
  vect<-regression_midas(data,arg,60)
  d<-2
  #vect<-res[,c(1,2)]
  #vect est un tableau à deux colonnes : 
  # - la première indique le numéro du coefficient
  # - la seconde indique sa valeur
  tableau<-data.frame(vect)
  p <- ggplot(tableau, aes(y=vect[,d], x = vect[,1]))
  p + geom_line(aes(y=vect[,d], colour="400 jours"))+
    xlab("Jours avant la fin de l'observation") +
    ylab("Valeur des coefficients") +
    theme(legend.title = element_blank())
}



loop_rmsfe<-function(h,z,data,arg,init){  
  #data : data set that contains the data we are working on
  #arg : arguments vector. For each regressor, it contains the its name, the number of days, the degree of the constaints polynom and the method (polynomial, exponential)
  #h : number of days of lag between the end of the estimation period and the end of the quarter we want to estimate.
  #if h=60 for example, we will forecast the values 60 days before their release
  #init : vector of initial values for minimization (the exponential case)
  #This function is split in two parts: 
  # - the first part aims at computing a regression
  # - the second part aims at computing the daily coefficients, thanks to the regression coefficients
  
  
  
  #I put the argument vector ("arg") into a matrix
  arguments<-matrix(arg,ncol=5,byrow=TRUE)
  #If the method for one variable is umidas, then I replace the degree by the number of lag days
  arguments[arguments[,4]=="umidas",3]<-(as.numeric(arguments[arguments[,4]=="umidas",2])-1)
  
  variables<-c(arguments[,1])
  retards<-as.numeric(c(arguments[,2]))
  degres<-as.numeric(c(arguments[,3]))
  method<-c(arguments[,4])
  hf<-c(arguments[,5])
  nb_var<-nrow(arguments)
  
  
  ############################
  #I compute the regression###
  ############################
  #I create a table that contains the data for the regression
  #that's why I call the "base" function
  
  #I regress.
  
  forecast<-c()
  upper_bound<-c()
  lower_bound<-c()
  prov2<-base_reg(data,arg,h)
  prov2<-data.frame(prov2[,-1])
  prov3<-prov2
  for (j in 1:(z-1)){  
    nb_row<-(nrow(prov3)+j-z)
    nb_row2<-nb_row+1
    if (all(method==rep("exponential",nb_var))){
      
      
      prov2<-prov3[1:nb_row,]
      
      #error message if the number of starting values and the number of parameters to estimate don't match
      if (length(init)>(sum(degres)+1)){
        print("There are more initial values than estimated parameters")
      }
      if (length(init)<(sum(degres)+1)){
        print("There are less initial values than estimated parameters")
      }
      
      #if the method specified is non-linear, I first write the formula of the function to be minimized
      expr1<-paste("+",apply(arguments,1,expression),collapse="",sep="")
      expr1<-paste("const",expr1,sep="")
      #and then I compute a non linear minimization
      reg<-nlm(fc,init,expr=expr1,variables=variables,deg=degres,data=prov)
      #I stock the estimates in a vector with a standard name
      estimated_coeffs<-reg$estimate  
      
    }else{
      
      #else (ie if the method is polynomial or umidas), I compute a linear model
      reg<-summary(lm(pib~.,data=prov2))
      #I stock the estimates in a vector with a standard name
      estimated_coeffs<-reg$coefficient[,1]
      sigma2<-reg$residuals%*%reg$residuals/(nrow(prov2)-ncol(prov2)+1)
    }
    
    ###################################
    #I compute the daily coefficients##
    ###################################
    #I create a matrix whose purpose is to receive the values of the daily coefficients
    coefficients<-matrix(NA,nrow=max(retards),ncol=(length(retards)+1))
    #the first column of this matrix contains the day id.
    coefficients[,1]<-c(1:max(retards))
    
    #exponential case
    if (all(method==rep("exponential",nb_var))){
      #I first have to assign the estimated values to the matching coefficients
      #that's why I compute a loop that gives their values to each subset of "estimated_coeffs"
      #I first assign a value to the intercept, named "const"
      const<-reg$estimate[1]  
      #I assign an initial value to a variable, very useful for this loop
      s<-2
      for (i in 1:length(variables)){
        #I select a subset of the estimated vector, which matches to the corresponding coefficients
        values<-estimated_coeffs[s:(degres[i]+1)]
        #I increase the value of s, in order to select the proprer subset of estimates during the next iteration
        s<-(degres[i]+2)
        #I compute the daily coefficients by computing a matrix product and then by dividing each coefficient by the sum of coefficients
        M<-Mat(retards[i],degres[i])
        prov2<-sapply(M[,-1]%*%values,exp)/sum(sapply(M[,-1]%*%values,exp))
        #I put the daily coefficients (prov2) into a matrix containing the whole set of coefficients. 
        #the first column of the matrix contains the day id
        coefficients[,i+1]<-c(prov2,rep(NA,(nrow(coefficients)-length(prov2))))
      }  
    }else{
      #polynomial case
      const<-estimated_coeffs[1]
      #n is a useful variable for the following loop. I have to initialise it.
      n<-1
      for (i in 1:length(variables)){
        #At each iteration I need to select the m-th to the n-th coefficients stocked into reg$coefficients
        #the value of m is 1 plus the value of the last coefficient selected during the previous iteration
        m<-(n+1)
        #the value of n is m plus the number of coefficients i need to select, ie degres[i]
        n<-(m+degres[i])
        #I prepare the transformation matrix
        M<-Mat(retards[i],degres[i])
        #I put the computed daily coefficients into one column of the "coefficients" matrix.
        if (method[i]=="polynomial"){
          #if the constraint is polynomial, I have to compute a matrix product with the estimated coefficients and the polynomial matrix
          coefficients[,(i+1)]<-c(M%*%reg$coefficients[m:n,1],rep(NA,(length(coefficients[,(i+1)])-length(M%*%reg$coefficients[m:n,1]))))
        }
        if (method[i]=="umidas"){
          #if the model is unconstrained, I just put the coefficents into the matrix, without any matrix product.
          coefficients[,(i+1)]<-c(reg$coefficients[m:n,1],rep(NA,(length(coefficients[,(i+1)])-length(reg$coefficients[m:n,1]))))
        }
      }
    }
    
    forecast[j]<-reg$coefficients[1,1]+as.matrix(prov3[nb_row2,2:ncol(prov3)])%*%(as.matrix(reg$coefficients[2:nrow(reg$coefficients),1]))
    
    #block that computes the confidence intervals of the forecasting values
    X<-as.matrix(prov2[,-1])
    #forecasting error variance
    fev<-sqrt(sigma2*(1+as.matrix(prov3[nb_row2,2:ncol(prov3)])%*%ginv(t(X)%*%X)%*%t(as.matrix(prov3[nb_row2,2:ncol(prov3)]))))
    kk<-qt(0.90,nrow(X)-ncol(X))*fev
    upper_bound[j]<-forecast[j]+kk
    lower_bound[j]<-forecast[j]-kk
    
    
  }
  
  abcisses<-1:length(forecast)
  realite<-c(prov3[(nrow(prov3)+1-length(forecast)):nrow(prov3),1])
  graph_rmsfe<-data.frame(abcisses,realite,forecast,upper_bound,lower_bound)
  rmsfe<-sqrt(mean((graph_rmsfe[,3]-graph_rmsfe[,2])^2))
  
  
  courbe<-ggplot()
  courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=realite,colour="realite"))
#   courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=lower_bound,colour="bb"))
#   courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=upper_bound,colour="bh"))  
  courbe<-courbe + geom_path(data=graph_rmsfe,aes(x=abcisses,y=forecast,colour="prevision"))
  courbe  
  
  
  
  return(list(graph_rmsfe,rmsfe,courbe,sigma2))
}



save(list = ls(all=TRUE), file = "functions.RData")


