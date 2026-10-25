# print shows the three tests and every edge

    Code
      print(nct_fit(1))
    Output
      Network Comparison Test: group 1 vs group 2
      15 permutations, independent groups, n = 120 and 120
      
        Test                  Statistic       p
        M  network structure      0.370   0.062 
        S  global strength        0.658   0.062 
             group 1 2.170 vs group 2 1.513
        E  edges differing      0 of 10
      
        Edges, largest difference first
        Edge     group 1  group 2  |diff|       p
        V1 - V2    0.238    0.608   0.370   0.062 
        V1 - V3    0.321    0.000   0.321   0.062 
        V2 - V4    0.314    0.000   0.314   0.062 
        V2 - V3    0.060    0.298   0.237   0.062 
        V3 - V5    0.315    0.130   0.185   0.188 
        V3 - V4    0.188    0.021   0.167   0.188 
        V1 - V5    0.254    0.108   0.145   0.312 
        V2 - V5    0.185    0.062   0.123   0.438 
        V1 - V4    0.206    0.161   0.045   0.875 
        V4 - V5    0.090    0.125   0.035   0.875 
      
        * p < 0.05

