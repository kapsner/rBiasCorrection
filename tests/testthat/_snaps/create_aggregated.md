# test functioning of aggregated function

    Code
      round(table_prep(df_agg), 1)
    Output
         true_methylation   CpG    sd CpG_true_diff relative_error
                    <num> <num> <num>         <num>          <num>
      1:              0.0   3.0   0.8           3.0             NA
      2:             12.5  10.3   3.0           2.2           17.8
      3:             25.0  17.3   1.8           7.7           30.7
      4:             37.5  26.2   3.6          11.3           30.1
      5:             50.0  36.8   2.7          13.2           26.3
      6:             62.5  48.3   2.7          14.2           22.7
      7:             75.0  57.8   4.5          17.2           22.9
      8:             87.5  65.0   2.0          22.5           25.7
      9:            100.0  93.0   0.7           7.0            7.0

---

    Code
      round(table_prep(df_agg), 1)
    Output
            CpG    sd
          <num> <num>
       1:  63.1  51.0
       2:  20.6   6.9
       3:  28.7  17.5
       4:  41.1  26.4
       5:   8.7    NA
       6:  15.7  10.4
       7:  16.5    NA
       8:  31.5   7.7
       9:  70.0  13.8
      10:   2.3    NA

---

    Code
      round(table_prep(df_agg), 2)
    Output
            CpG    sd
          <num> <num>
       1: 20.63  6.88
       2: 69.98 13.77
       3: 16.48    NA
       4: 63.14 50.99
       5: 29.61 16.16
       6:  2.29    NA
       7: 31.51  7.67
       8: 41.10 26.39
       9:  8.69    NA
      10: 15.72 10.42

