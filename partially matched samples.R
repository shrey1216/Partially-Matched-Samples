#function partially_matched() which implements both Kim's and Looney and Jones's statistical procedures.
partially_matched <- function(a,b, method = c("Kim", "Looney"), 
                              alternative = c("two.sided", "less", "greater"),
                              conf.level = 0.95) {
  #putting vectors a and b into one matrix, and declaring data1, data2, data3, and data4
  data <- matrix(c(a,b), nrow = 2, byrow = T) 
  data1 <- c(); data2 <- c(); data3 <- c(); data4 <- c()

  
  

  #sorting the data by which values are missing into data1, data2, data3, and data4
  for(i in 1:length(a)) {                                           
    if(identical((is.na(data[,i]) == c(FALSE, FALSE)), c(TRUE, TRUE))){
      data1 <- c(data1,data[,i])
    } else if (identical((is.na(data[,i]) == c(FALSE, TRUE)), c(TRUE, TRUE))){
      data2 <- c(data2,data[,i])
    }  else if (identical((is.na(data[,i]) == c(TRUE, FALSE)), c(TRUE, TRUE))){
      data3 <- c(data3,data[,i])
    }  else if (identical((is.na(data[,i]) == c(TRUE, TRUE)), c(TRUE, TRUE))){
      data4 <- c(data4,data[,i])
    }                                         
  }
  
  data1 <- matrix(c(data1), nrow = 2, byrow = F)
  data2 <- matrix(c(data2), nrow = 2, byrow = F)
  data3 <- matrix(c(data3), nrow = 2, byrow = F)
  data4 <- matrix(c(data4), nrow = 2, byrow = F)
  
  
  
  
  #mean, standard deviation, n, error, lower and upper bounds, and df for the output
  sample_mean <- mean(data, na.rm = TRUE)
  sample_sd <- sd(data, na.rm = TRUE)
  sample_n <- length(data1) + length(data2[1,]) + length(data3[2,])
  error <- qnorm(1 - ((1 - conf.level) / 2)) * (sample_sd/sqrt(sample_n))
  lower_bound = sample_mean - error
  upper_bound = sample_mean + error
  percent <- conf.level * 100
  df = sample_n - 1
  
  
  #based on the method in the input
  if (method == "Kim"){
    #Kim et al.’s Modified t-Statistic
    
    #n1, n2, n3, and n4 are the number of paired samples in each set 
    n1 <- ncol(data1); n2 <- ncol(data2); n3 <- ncol(data3); n4 <- ncol(data4)
    
    
    #Let D and SD be the mean and the standard deviation of the difference of Sample 1 and Sample
    #2 in data1, respectively

    D <-  mean(data1[1,] - data1[2,])
    SD <- sd(data1[1,] - data1[2,])
    
    #Let T and ST be the mean and the standard deviation of Sample 1 in data2, respectively
    
    T_ <- mean(data2[1,])
    ST <- sd(data2[1,])
    
    #Define N and SN for Sample 2 in data3, respectively
    
    N <- mean(data3[2,])
    SN <- sd(data3[2,])
    
    #Let nH be the harmonic mean of n2 and n3
    
    nH <- harmonic.mean(n2, n3)

    
    #Plugging them into the formula
    numerator <- n1 * D + nH * (T_ - N)
    denominator <- sqrt(n1*(SD**2) + (nH**2)*((ST**2)/n2 + (SN**2)/n3))
    t3 <- numerator/denominator
    
    
    #calculating the p-value
    p_value <- 0
    if (alternative == "two.sided") {
      p_value <- 2*(1- pnorm(abs(t3)))
    } else if(alternative == "greater") {
      p_value <- 1 - pnorm(t3)
    } else if(alternative == "less") {
      p_value <- pnorm(t3)
    } else {
      stop("argument should be one of the following: \"two.sided\", \"greater\", \"less\"")
    }
    
    #printing results
    cat("\n     Kim et al Modified t-Statistic\n\n      df:", df,"\n      p-value:", p_value, "\n      t-statistic:", t3, "\n     ", percent,"percent confidence interval:\n      ", lower_bound,"", upper_bound,
        "\n      sample estimates:\n      mean of x \n             ", sample_mean)
    
    
  } else if (method == "Looney"){
    #Looney and Jones’s Corrected Z-Test
    
    
    #n1, n2, n3, and n4 are the number of paired samples in each set 
    n1 <- ncol(data1); n2 <- ncol(data2); n3 <- ncol(data3); n4 <- ncol(data4)
    
    
    #Let T∗ and ST∗ be the mean and the standard deviation of Sample 1 in the combined 
    #data of data1 and data2, respectively.

    T1 <- mean(c(data1[1,],data2[1,]))
    ST1 <- sd(c(data1[1,],data2[1,]))
    
    #Define N∗ and SN∗ for Sample 2 in the combined data of data1 and data3, respectively.

    N1 <- mean(c(data1[2,],data3[2,]))
    SN1 <- sd(c(data1[2,],data3[2,]))
    
    #Let STN1 be the sample covariance of Sample 1 and Sample 2 in data1.
    STN1 <- cov(data1[1,], data1[2,])
    
    #All the values have been defined.

    numerator <- T1 - N1
    denominator <- sqrt((ST1**2)/(n1 + n2) + (SN1**2)/(n1 + n3) - (2 * n1 * STN1)/((n1 + n2) * (n1 + n3)))
    z_corr <- numerator/denominator
    
    #calculating the p-value
    p_value <- 0
    if (alternative == "two.sided") {
      p_value <- 2 * (1 - pnorm(abs(z_corr)))
    } else if(alternative == "greater") {
      p_value <- 1 - pnorm(z_corr)
    } else if(alternative == "less") {
      p_value <- pnorm(z_corr)
    } else {
      stop("argument should be one of the following: \"two.sided\", \"greater\", \"less\"")
    }
    
    #printing results
    cat("\n     Looney and Jones's Corrected Z-Test\n\n      df:", df,"\n      p-value:", p_value, "\n      Z-corrected:", z_corr, "\n     ", percent,"percent confidence interval:\n      ", lower_bound,"", upper_bound,
        "\n      sample estimates:\n      mean of x \n             ", sample_mean)
    
  }else{
    stop("argument should be one of the following: \"Kim\", \"Looney\"")
  }
  
}




#sample data with missing values
a <- c(1, NA, 3, 4, 5, 6, NA, NA, 9, NA, NA, NA)
b <- c(12, 11, 10, NA, NA, NA, 6, 5, 4, NA, NA, 20)

#test function
#methods = "Kim", "Looney"
#alternative = "two.sided", "less", "greater
#conf.level = 0.95 by default (any decimal percentage)
partially_matched(a,b, method = "Kim", alternative = "less", conf.level = 0.95)
partially_matched(a,b, method = "Looney", alternative = "less", conf.level = 0.95)
