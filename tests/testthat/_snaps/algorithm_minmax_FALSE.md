# algorithm test, type 1, minmax = FALSE

    Code
      round(table_prep(rv$reg_stats), 1)
    Output
          relative_error SSE_hyperbolic R2_hyperbolic a_hyperbolic b_hyperbolic
                   <num>          <num>         <num>        <num>        <num>
       1:           22.9           77.3             1       -108.6       -937.7
       2:           10.4           45.8             1       -141.9      -2393.1
       3:           23.2           67.1             1       -170.4      -1477.5
       4:           14.1           58.3             1       -162.2      -2115.9
       5:           25.3            7.5             1       -151.4      -1383.3
       6:           13.6           11.8             1       -125.0      -1550.6
       7:           43.6           72.1             1        -75.3       -772.8
       8:            9.4           75.3             1       -278.0      -3872.7
       9:           25.8           33.1             1        -70.7       -843.9
      10:           19.0           35.4             1       -127.7      -1525.8
          d_hyperbolic b1_hyperbolic s_hyperbolic SSE_cubic R2_cubic a_cubic b_cubic
                 <num>         <num>        <num>     <num>    <num>   <num>   <num>
       1:       -232.1           0.6          4.1      71.0        1       0       0
       2:       -268.3           0.6          8.9      46.0        1       0       0
       3:       -327.2           0.7          4.5      66.8        1       0       0
       4:       -299.7           0.7          7.1      57.2        1       0       0
       5:       -304.2           0.7          4.6       6.4        1       0       0
       6:       -247.3           0.6          6.3      12.3        1       0       0
       7:       -226.7           0.6          3.5      69.0        1       0       0
       8:       -438.9           0.8          8.8      73.3        1       0       0
       9:       -184.5           0.5          4.6      34.5        1       0       0
      10:       -263.0           0.6          5.8      36.0        1       0       0
          c_cubic d_cubic
            <num>   <num>
       1:     0.8     1.9
       2:     0.5     9.3
       3:     0.6     3.9
       4:     0.7     6.2
       5:     0.4     5.2
       6:     0.6     6.1
       7:     0.6     1.7
       8:     0.5     9.6
       9:     0.7     2.7
      10:     0.6     5.2

---

    Code
      round(table_prep(rv$final_results), 1)
    Output
          CpG#1_c CpG#2_h CpG#3_c CpG#4_c CpG#5_c CpG#6_h CpG#7_c CpG#8_c CpG#9_h
            <num>   <num>   <num>   <num>   <num>   <num>   <num>   <num>   <num>
       1:      NA    76.7    99.5      NA    69.3    79.3      NA    69.5    81.0
       2:    31.4    31.0    28.0    29.0    26.8    30.2    24.1    28.7    26.6
       3:      NA    42.7    45.1    49.4    43.5    41.8      NA    35.5    45.7
       4:      NA    58.9    65.1    77.4    51.9    56.8      NA    57.5    57.0
       5:     9.3     5.0     8.1    10.2    10.2     8.9    11.2     9.9     6.8
       6:    20.9    20.2    24.1    22.1    23.4    18.7    24.0    20.3    22.9
       7:    22.4    17.8    26.4    27.6    26.4    29.9    28.9    38.9    38.9
       8:      NA    40.9    48.1    47.0    42.2    42.8    68.1    44.6    43.7
       9:      NA    86.0      NA      NA    84.1    86.8      NA    79.7    88.0
      10:     0.5     3.1     0.4     3.3     1.4     1.5     0.5     0.0     0.0
          row_means_h
                <num>
       1:        77.1
       2:        28.4
       3:        42.5
       4:        57.3
       5:         8.8
       6:        21.8
       7:        28.7
       8:        43.4
       9:        86.4
      10:         0.0

---

    Code
      round(table_prep(rv$substitutions), 1)
    Output
      Null data.table (0 rows and 0 cols)

---

    Code
      round(table_prep(solved_eq2[["results"]]), 1)
    Output
         CpG#1_c CpG#2_h CpG#3_c CpG#4_c CpG#5_c CpG#6_h CpG#7_c CpG#8_c CpG#9_h
           <num>   <num>   <num>   <num>   <num>   <num>   <num>   <num>   <num>
      1:     1.3     1.1     1.4     0.9     1.5     0.1     2.1     2.0     0.0
      2:    11.7    11.4    10.6    12.3    10.3    11.9    10.0     9.5    10.5
      3:    24.1    26.2    25.6    24.4    24.5    26.5    23.1    26.4    26.9
      4:    50.4    35.1    38.3    41.6    37.1    35.3    42.6    34.7    36.8
      5:      NA    47.7    57.4    57.9    49.5    50.1      NA    50.9    52.0
      6:      NA    67.1    80.1   100.0    59.9    65.0      NA    62.3    64.9
      7:      NA    75.8   100.0      NA    72.6    73.7      NA    74.2    74.6
      8:      NA    84.4      NA      NA    81.5    87.1      NA    76.7    84.5
      9:      NA   100.0      NA      NA    94.9   100.0      NA    94.3   100.0
         row_means_h
               <num>
      1:         0.3
      2:        11.0
      3:        25.4
      4:        36.5
      5:        50.7
      6:        65.3
      7:        75.5
      8:        83.2
      9:       100.0

---

    Code
      round(table_prep(solved_eq2[["substitutions"]]), 1)
    Output
      Null data.table (0 rows and 0 cols)

---

    Code
      round(table_prep(rv$fileimport_cal_corrected), 1)
    Output
         CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
         <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
      1:   1.3   1.1   1.4   0.9   1.5   0.1   2.1   2.0   0.0       0.3
      2:  11.7  11.4  10.6  12.3  10.3  11.9  10.0   9.5  10.5      11.0
      3:  24.1  26.2  25.6  24.4  24.5  26.5  23.1  26.4  26.9      25.4
      4:  50.4  35.1  38.3  41.6  37.1  35.3  42.6  34.7  36.8      36.5
      5:    NA  47.7  57.4  57.9  49.5  50.1    NA  50.9  52.0      50.7
      6:    NA  67.1  80.1 100.0  59.9  65.0    NA  62.3  64.9      65.3
      7:    NA  75.8 100.0    NA  72.6  73.7    NA  74.2  74.6      75.5
      8:    NA  84.4    NA    NA  81.5  87.1    NA  76.7  84.5      83.2
      9:    NA 100.0    NA    NA  94.9 100.0    NA  94.3 100.0     100.0

---

    Code
      round(table_prep(rv$fileimport_cal_corrected_h), 1)
    Output
         CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
         <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
      1:   0.0   1.1   0.5   0.0   2.4   0.1   0.0   2.8   0.0       0.3
      2:  12.2  11.4  10.8  12.5  10.2  11.9  10.2   9.3  10.5      11.0
      3:  24.5  26.2  25.5  24.3  24.0  26.5  24.6  25.5  26.9      25.4
      4:  38.2  35.1  36.5  38.1  37.3  35.3  37.8  34.0  36.8      36.5
      5:  52.3  47.7  50.8  48.6  50.9  50.1  53.6  51.8  52.0      50.7
      6:  65.5  67.1  64.9  67.7  62.4  65.0  65.9  64.7  64.9      65.3
      7:  75.0  75.8  77.6  74.2  76.4  73.7  75.7  78.4  74.6      75.5
      8:  81.5  84.4  80.4  82.9  86.2  87.1  79.4  81.3  84.5      83.2
      9: 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0 100.0     100.0

---

    Code
      round(table_prep(rv$substitutions_corrected_h), 1)
    Output
      Null data.table (0 rows and 0 cols)

---

    Code
      round(table_prep(rv$reg_stats_corrected_h), 1)
    Output
          relative_error SSE_hyperbolic R2_hyperbolic a_hyperbolic b_hyperbolic
                   <num>          <num>         <num>        <num>        <num>
       1:            2.8           39.6             1        985.3       -657.9
       2:            4.5           46.2             1      30221.7       3393.1
       3:            4.4           59.9             1       1333.7       -945.8
       4:            2.8           46.9             1       1637.4       -754.1
       5:            3.6           16.6             1      -5621.1      -1917.3
       6:            2.9           15.4             1      13767.3      -1255.0
       7:            5.5           76.6             1        722.0      -1065.0
       8:            7.0           86.2             1       4784.8         68.5
       9:            4.6           21.6             1       1273.9      -1022.8
      10:            3.4           26.7             1       1740.6      -1166.7
          d_hyperbolic b1_hyperbolic s_hyperbolic SSE_cubic R2_cubic a_cubic b_cubic
                 <num>         <num>        <num>     <num>    <num>   <num>   <num>
       1:        903.5           1.1          0.7      39.5        1       0       0
       2:      30297.8           1.0          0.1      42.6        1       0       0
       3:       1253.1           1.1          0.8      59.1        1       0       0
       4:       1558.4           1.1          0.5      46.9        1       0       0
       5:      -5734.1           1.0          0.3      11.1        1       0       0
       6:      13682.2           1.0          0.1      15.3        1       0       0
       7:        634.4           1.2          1.7      76.5        1       0       0
       8:       4739.0           1.0          0.0      74.4        1       0       0
       9:       1182.0           1.1          0.9      21.5        1       0       0
      10:       1653.9           1.1          0.7      26.3        1       0       0
          c_cubic d_cubic
            <num>   <num>
       1:     1.1    -0.9
       2:     0.9     1.0
       3:     1.0    -0.4
       4:     1.0    -0.4
       5:     0.8     1.4
       6:     1.0     0.0
       7:     1.1    -1.6
       8:     0.7     1.5
       9:     1.1    -1.0
      10:     1.0    -0.4

---

    Code
      round(table_prep(rv$fileimport_cal_corrected_c), 1)
    Output
         CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
         <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
      1:   1.3   0.5   1.4   0.9   1.5   0.3   2.1   2.0   1.5       1.3
      2:  11.7  11.5  10.6  12.3  10.3  11.8  10.0   9.5  10.2      11.0
      3:  24.1  26.8  25.6  24.4  24.5  27.0  23.1  26.4  27.0      26.0
      4:  50.4  35.9  38.3  41.6  37.1  36.8  42.6  34.7  45.9      39.5
      5:    NA  48.7  57.4  57.9  49.5  54.4    NA  50.9    NA      60.3
      6:    NA  69.2  80.1 100.0  59.9  73.9    NA  62.3    NA      87.2
      7:    NA  78.6 100.0    NA  72.6  86.3    NA  74.2    NA        NA
      8:    NA  88.3    NA    NA  81.5 100.0    NA  76.7    NA        NA
      9:    NA 100.0    NA    NA  94.9    NA    NA  94.3    NA        NA

---

    Code
      round(table_prep(rv$substitutions_corrected_c), 1)
    Output
      Null data.table (0 rows and 0 cols)

---

    Code
      round(table_prep(rv$reg_stats_corrected_c), 1)
    Output
          relative_error SSE_hyperbolic R2_hyperbolic a_hyperbolic b_hyperbolic
                   <num>          <num>         <num>        <num>        <num>
       1:           14.9            1.7             1        -36.9       -124.7
       2:            4.7           49.8             1       1860.6      -1362.0
       3:           16.0           14.5             1       -189.8        -51.5
       4:           18.2           31.4             1        -60.7       -294.7
       5:            5.1           11.3             1       1781.8        210.3
       6:           10.3           33.6             1       -866.6        641.7
       7:           13.8            0.3             1        -43.4       -152.6
       8:            7.3           71.5             1        801.2        -54.7
       9:           16.3            3.6             1        -74.5        -84.7
      10:           16.3            3.4             1       -105.6       -180.0
          d_hyperbolic b1_hyperbolic s_hyperbolic SSE_cubic R2_cubic a_cubic b_cubic
                 <num>         <num>        <num>     <num>    <num>   <num>   <num>
       1:        -67.5          -0.5          1.9       0.0        1       0       0
       2:       1724.8           1.1          0.8      35.7        1       0       0
       3:       -216.2           0.5          0.3       8.3        1       0       0
       4:       -103.8           0.0          2.9      25.3        1       0       0
       5:       1791.0           1.1          0.1      10.8        1       0       0
       6:       -827.1           0.9          0.8      21.2        1       0       0
       7:        -79.2          -0.3          2.0       0.0        1       0       0
       8:        765.3           1.1          0.1      68.9        1       0       0
       9:        -99.9           0.0          0.9       0.0        1       0       0
      10:       -140.4           0.3          1.3       3.6        1       0       0
          c_cubic d_cubic
            <num>   <num>
       1:     1.1     1.3
       2:     0.8     0.8
       3:     0.6     1.5
       4:     1.4     0.2
       5:     0.9     0.4
       6:     0.7     0.9
       7:     0.4     2.1
       8:     0.9     0.6
       9:     0.2     1.5
      10:     0.9     1.0

