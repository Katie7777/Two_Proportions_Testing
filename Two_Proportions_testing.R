################# 2*2 Contingency Table P value #################
gen_contigTab_pval <-
  function(x1, x2, n1, n2, test.type = 'two_sample_prop') {
    #test.type: chisq, cmh_unstrat, fisher_exact, two_sample_prop
    #xi: number of events in group i
    #ni: total number of pats in group i
    non_res_1 = n1 - x1
    non_res_2 = n2 - x2
    
    E1 = n1 * (x1 + x2) / (n1 + n2)
    E2 = n2 * (x1 + x2) / (n1 + n2)
    E3 = (n1 + n2 - x1 - x2) * n1 / (n1 + n2)
    E4 = (n1 + n2 - x1 - x2) * n2 / (n1 + n2)
    
    if (test.type == "chisq") {
      #chi square statistic formula
      chisq_ = (x1 - E1) ^ 2 / E1 + (x2 - E2) ^ 2 / E2 +
        (non_res_1 - E3) ^ 2 / E3 + (non_res_2 - E4) ^ 2 / E4
      
      pval = 1 - pchisq(chisq_, 1)
      
      #testing using function in R-chisq; confirmed
      data <- matrix(c(x1, n1-x1, x2, n2-x2), ncol=2, byrow = TRUE)
      data <-as.table(data)
      # by default, chisq.test() will do continuity correction.
      # correction will be more conservative: larger pval with fewer rejections
      pval_ = chisq.test(data, correct = F)$p.value
      
      # #simplified formula (identical to chisq_)
      # chisq_2 = 2 * (x1^2 + x2^2)/(x1+x2) +
      #   2 * (non_res_1^2 + non_res_2^2)/(non_res_1 + non_res_2) - (n1+n2)
    } else if (test.type == "cmh_unstrat") {
      #non-stratefied cmh test (very minor differece from chi square test)
      cmh_ = (x1 - n1 * (x1 + x2) / (n1 + n2)) ^ 2 /
        (n1 * n2 * (x1 + x2) * (non_res_1 + non_res_2) / (n1 + n2) ^ 2 /
           (n1 + n2 - 1))
      
      pval = 1 - pchisq(cmh_, 1)
    } else if (test.type == "fisher_exact") {
      #Fisher exact test
      pval = fisher.test(rbind(c(x1, x2), c(non_res_1, non_res_2)),
                         alternative = "greater")$p.value
    } else if (test.type == "two_sample_prop") {
      #two sample proportion test
      p1_hat = x1 / n1
      p2_hat = x2 / n2
      
      #compare arm1 and arm3
      z_stat = (p1_hat - p2_hat) / (sqrt(p1_hat * (1 - p1_hat) / n1 + p2_hat * (1 - p2_hat) /
                                           n2))
      pval = 1 - pnorm(z_stat)
    }
    
    return(pval)
  }


################# Simulation #########################
one_sim_contigTab <-
  function(n, p1, p2, p3, alpha, test.type = 'two_sample_prop') {
    #test.type: chisq, cmh_unstrat, fisher_exact, two_sample_prop
    n1 = floor(n / 3)
    n2 = floor(n / 3)
    n3 = n - n1 - n2
    
    x1 = rbinom(1, n1, p1)
    # print(c("x1",x1))
    x2 = rbinom(1, n2, p2)
    # print(c("x2",x2))
    x3 = rbinom(1, n3, p3)
    # print(c("x3",x3))
    

    #compare arm1 and arm3
    p_value_13 = gen_contigTab_pval(x1, x3, n1, n3, test.type)
    
    
    #compare arm2 and arm3
    p_value_23 = gen_contigTab_pval(x2, x3, n2, n3, test.type)
    
    #family wise error rate (Hochberg)
    result = 0
    rank_pvalue = sort(c(p_value_13, p_value_23), decreasing = T)
    if (rank_pvalue[1] < alpha) {
      result = 1
    } else if (rank_pvalue[2] < alpha / 2) {
      result = 1
    }
    
    return(result)
  }


run_simulation <-
  function(n_sim, n, p1, p2, p3, alpha, test.type = 'two_sample_prop') {
    sum = 0
    for (i in 1:n_sim) {
      res = one_sim_contigTab(
        n = n,
        p1 = p1,
        p2 = p2,
        p3 = p3,
        alpha = alpha,
        test.type
      )
      sum = sum + res
    }
    power = sum / n_sim
    return(power)
  }


######################## RUN SIMULATION #######################
#test.type: chisq, cmh_unstrat, fisher_exact, two_sample_prop

################# TEST ALPHA #########################
alpha = 0.1
p1 = 0.5
p2 = 0.5
p3 = 0.5
n_sim = 10000
n = 120
test.type = 'chisq'

run_simulation(n_sim, n, p1, p2, p3, alpha, test.type)

################# TEST POWER #########################
alpha = 0.1
p1 = 0.7
p2 = 0.6
p3 = 0.5
n_sim = 10000
n = 120
test.type = 'chisq'

run_simulation(n_sim, n, p1, p2, p3, alpha, test.type)
