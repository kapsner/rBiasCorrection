# test normal function of file import of type 1

    Code
      round(table_prep(exp_type_1[["dat"]]), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
       1:   4.3   7.9   3.9   4.7   5.4   6.0   3.0  12.6   4.9       5.8
       2:  92.5 100.0  80.8  90.0  73.6  98.7  66.3  85.4  88.7      86.2
       3:  92.6 100.0  84.7 100.0  81.9 100.0  65.4 100.0  97.0      91.3
       4:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6
       5:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7
       6:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5
       7:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7
       8:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2
       9:   2.5  11.6   6.7   8.5   7.4   7.0   4.8  10.2   4.7       7.0
      10:  16.8  31.5  22.2  23.4  20.5  23.0  14.1  21.7  18.6      21.3
      11:  22.8  20.7  19.3  31.2  23.7  20.2  17.4  22.7  19.0      21.9
      12:    NA    NA    NA    NA  35.2  38.4  31.0  35.3  29.1      33.8
      13:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6
      14:  45.7  62.0  47.5  57.6  46.6  50.7  35.1  48.3  46.0      48.8
      15:  45.3  51.2  47.5  51.5  44.9  48.0  35.0  54.6  46.6      47.2
      16:   2.6  10.1   4.3   7.7   4.4   7.1   2.3  10.5   2.6       5.7
      17:  49.6  61.1  45.2  53.8  49.0  61.1  28.0  79.2  30.2      50.8
      18:  49.5  56.9  49.6  57.1  34.4  44.2  39.4  51.6  49.5      48.0
      19:  62.7  77.4  63.0  73.9  63.6  69.4  53.0  72.8  70.5      67.4
      20:   8.7  11.8   8.9  13.1   9.9  11.2   7.7  15.3   7.5      10.4
      21:   8.3  16.6  14.0  16.0  12.2   9.5   8.7  13.9   8.3      11.9
      22:  23.1  25.9  22.9  25.2  22.2  24.6  18.0  30.2  22.2      23.8
      23:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5  24.7      22.2
      24:  19.2  25.9  19.0  23.2  18.2  21.9  10.8  28.3  10.7      19.7
      25:  30.0  36.0  31.0  32.3  28.3  32.2  23.4  32.4  33.9      31.1
      26:  33.7  35.3  34.7  36.8  31.4  35.6  27.2  41.0  34.8      34.5
      27:  39.3  47.0  39.0  41.9  41.6  43.7  22.8  66.7  23.8      40.7
      28:  35.4  36.1  35.2  38.4  31.0  35.3  29.1  38.7  36.4      35.1
      29:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2
      30:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9
      31:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3
      32:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7
      33:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6
      34:   2.3  10.7   4.2   8.4   5.8   7.1   2.0   7.5   3.2       5.7
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means

---

    Code
      round(table_prep(exp_type_1), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
       1:   4.3   7.9   3.9   4.7   5.4   6.0   3.0  12.6   4.9       5.8
       2:  92.5 100.0  80.8  90.0  73.6  98.7  66.3  85.4  88.7      86.2
       3:  92.6 100.0  84.7 100.0  81.9 100.0  65.4 100.0  97.0      91.3
       4:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6
       5:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7
       6:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5
       7:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7
       8:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2
       9:   2.5  11.6   6.7   8.5   7.4   7.0   4.8  10.2   4.7       7.0
      10:  16.8  31.5  22.2  23.4  20.5  23.0  14.1  21.7  18.6      21.3
      11:  22.8  20.7  19.3  31.2  23.7  20.2  17.4  22.7  19.0      21.9
      12:    NA    NA    NA    NA  35.2  38.4  31.0  35.3  29.1      33.8
      13:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6
      14:  45.7  62.0  47.5  57.6  46.6  50.7  35.1  48.3  46.0      48.8
      15:  45.3  51.2  47.5  51.5  44.9  48.0  35.0  54.6  46.6      47.2
      16:   2.6  10.1   4.3   7.7   4.4   7.1   2.3  10.5   2.6       5.7
      17:  49.6  61.1  45.2  53.8  49.0  61.1  28.0  79.2  30.2      50.8
      18:  49.5  56.9  49.6  57.1  34.4  44.2  39.4  51.6  49.5      48.0
      19:  62.7  77.4  63.0  73.9  63.6  69.4  53.0  72.8  70.5      67.4
      20:   8.7  11.8   8.9  13.1   9.9  11.2   7.7  15.3   7.5      10.4
      21:   8.3  16.6  14.0  16.0  12.2   9.5   8.7  13.9   8.3      11.9
      22:  23.1  25.9  22.9  25.2  22.2  24.6  18.0  30.2  22.2      23.8
      23:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5  24.7      22.2
      24:  19.2  25.9  19.0  23.2  18.2  21.9  10.8  28.3  10.7      19.7
      25:  30.0  36.0  31.0  32.3  28.3  32.2  23.4  32.4  33.9      31.1
      26:  33.7  35.3  34.7  36.8  31.4  35.6  27.2  41.0  34.8      34.5
      27:  39.3  47.0  39.0  41.9  41.6  43.7  22.8  66.7  23.8      40.7
      28:  35.4  36.1  35.2  38.4  31.0  35.3  29.1  38.7  36.4      35.1
      29:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2
      30:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9
      31:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3
      32:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7
      33:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6
      34:   2.3  10.7   4.2   8.4   5.8   7.1   2.0   7.5   3.2       5.7
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means

---

    Code
      round(table_prep(exp_type_1), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
       1:   4.3   7.9   3.9   4.7   5.4   6.0   3.0  12.6   4.9       5.8
       2:  92.5 100.0  80.8  90.0  73.6  98.7  66.3  85.4  88.7      86.2
       3:  92.6 100.0  84.7 100.0  81.9 100.0  65.4 100.0  97.0      91.3
       4:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6
       5:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7
       6:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5
       7:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7
       8:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2
       9:   2.5  11.6   6.7   8.5   7.4   7.0   4.8  10.2   4.7       7.0
      10:  16.8  31.5  22.2  23.4  20.5  23.0  14.1  21.7  18.6      21.3
      11:  22.8  20.7  19.3  31.2  23.7  20.2  17.4  22.7  19.0      21.9
      12:    NA    NA    NA    NA  35.2  38.4  31.0  35.3  29.1      33.8
      13:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6
      14:  45.7  62.0  47.5  57.6  46.6  50.7  35.1  48.3  46.0      48.8
      15:  45.3  51.2  47.5  51.5  44.9  48.0  35.0  54.6  46.6      47.2
      16:   2.6  10.1   4.3   7.7   4.4   7.1   2.3  10.5   2.6       5.7
      17:  49.6  61.1  45.2  53.8  49.0  61.1  28.0  79.2  30.2      50.8
      18:  49.5  56.9  49.6  57.1  34.4  44.2  39.4  51.6  49.5      48.0
      19:  62.7  77.4  63.0  73.9  63.6  69.4  53.0  72.8  70.5      67.4
      20:   8.7  11.8   8.9  13.1   9.9  11.2   7.7  15.3   7.5      10.4
      21:   8.3  16.6  14.0  16.0  12.2   9.5   8.7  13.9   8.3      11.9
      22:  23.1  25.9  22.9  25.2  22.2  24.6  18.0  30.2  22.2      23.8
      23:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5  24.7      22.2
      24:  19.2  25.9  19.0  23.2  18.2  21.9  10.8  28.3  10.7      19.7
      25:  30.0  36.0  31.0  32.3  28.3  32.2  23.4  32.4  33.9      31.1
      26:  33.7  35.3  34.7  36.8  31.4  35.6  27.2  41.0  34.8      34.5
      27:  39.3  47.0  39.0  41.9  41.6  43.7  22.8  66.7  23.8      40.7
      28:  35.4  36.1  35.2  38.4  31.0  35.3  29.1  38.7  36.4      35.1
      29:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2
      30:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9
      31:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3
      32:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7
      33:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6
      34:   2.3  10.7   4.2   8.4   5.8   7.1   2.0   7.5   3.2       5.7
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means

---

    Code
      round(table_prep(cal_type_1), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>
       1:   4.3   7.9   3.9   4.7   5.4   6.0   3.0  12.6   4.9       5.8
       2:   2.5  11.6   6.7   8.5   7.4   7.0   4.8  10.2   4.7       7.0
       3:   2.6  10.1   4.3   7.7   4.4   7.1   2.3  10.5   2.6       5.7
       4:   2.3  10.7   4.2   8.4   5.8   7.1   2.0   7.5   3.2       5.7
       5:   3.1   7.6   4.8   4.5   6.0   4.6   2.6  12.5   3.2       5.4
       6:   8.7  11.8   8.9  13.1   9.9  11.2   7.7  15.3   7.5      10.4
       7:   8.3  16.6  14.0  16.0  12.2   9.5   8.7  13.9   8.3      11.9
       8:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6
       9:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5  24.7      22.2
      10:  19.2  25.9  19.0  23.2  18.2  21.9  10.8  28.3  10.7      19.7
      11:  16.8  31.5  22.2  23.4  20.5  23.0  14.1  21.7  18.6      21.3
      12:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7
      13:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5
      14:  22.8  20.7  19.3  31.2  23.7  20.2  17.4  22.7  19.0      21.9
      15:  23.1  25.9  22.9  25.2  22.2  24.6  18.0  30.2  22.2      23.8
      16:  30.0  36.0  31.0  32.3  28.3  32.2  23.4  32.4  33.9      31.1
      17:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7
      18:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2
      19:  33.7  35.3  34.7  36.8  31.4  35.6  27.2  41.0  34.8      34.5
      20:  39.3  47.0  39.0  41.9  41.6  43.7  22.8  66.7  23.8      40.7
      21:  35.4  36.1  35.2  38.4  31.0  35.3  29.1  38.7  36.4      35.1
      22:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6
      23:  49.6  61.1  45.2  53.8  49.0  61.1  28.0  79.2  30.2      50.8
      24:  45.7  62.0  47.5  57.6  46.6  50.7  35.1  48.3  46.0      48.8
      25:  49.5  56.9  49.6  57.1  34.4  44.2  39.4  51.6  49.5      48.0
      26:  45.3  51.2  47.5  51.5  44.9  48.0  35.0  54.6  46.6      47.2
      27:  51.3  65.0  49.1  62.0  49.2  61.1  40.9  58.3  54.7      54.6
      28:  62.7  77.4  63.0  73.9  63.6  69.4  53.0  72.8  70.5      67.4
      29:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2
      30:  53.1  64.6  55.0  52.4  52.5  57.7  43.6  63.7  55.4      55.3
      31:  55.0  64.9  56.8  60.0  53.5  60.4  30.6  82.8  35.1      55.5
      32:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9
      33:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3
      34:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7
      35:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6
      36:  92.5 100.0  80.8  90.0  73.6  98.7  66.3  85.4  88.7      86.2
      37:  92.6 100.0  84.7 100.0  81.9 100.0  65.4 100.0  97.0      91.3
      38:  92.9  99.3  86.2  93.0  84.6  88.4  71.6 100.0 100.0      90.7
      39:  92.6 100.0  82.8  95.0  77.8  99.3  65.9  92.7  92.9      88.8
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means

# test normal function of file import of type 2

    Code
      round(table_prep(exp_type_2), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>     <num>
       1:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6         9
       2:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7         9
       3:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5         9
       4:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7         9
       5:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2         9
       6:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2         9
       7:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9         9
       8:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3         9
       9:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7         9
      10:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6         9
      11:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5    NA      21.9         8
      12:   4.3   7.9   3.9   4.7   5.4   6.0   3.0    NA    NA       5.0         7
      13:  92.5 100.0  80.8  90.0  73.6  98.7  66.3    NA    NA      86.0         7
      14:  92.6 100.0  84.7 100.0  81.9 100.0  65.4    NA    NA      89.2         7
      15:   2.5  11.6   6.7   8.5   7.4    NA    NA    NA    NA       7.3         5
      16:  16.8  31.5  22.2  23.4  20.5    NA    NA    NA    NA      22.9         5
      17:  22.8  20.7  19.3  31.2  23.7    NA    NA    NA    NA      23.5         5
      18:  35.2  38.4  31.0  35.3  29.1    NA    NA    NA    NA      33.8         5
      19:  38.9  47.7  37.8  42.1  39.5    NA    NA    NA    NA      41.2         5
      20:  45.7  62.0  47.5  57.6  46.6    NA    NA    NA    NA      51.9         5
      21:  45.3  51.2  47.5  51.5  44.9    NA    NA    NA    NA      48.1         5
      22:   2.3  10.7   4.2   8.4    NA    NA    NA    NA    NA       6.4         4
      23:  19.2  25.9  19.0    NA    NA    NA    NA    NA    NA      21.3         3
      24:  30.0  36.0  31.0    NA    NA    NA    NA    NA    NA      32.3         3
      25:  33.7  35.3  34.7    NA    NA    NA    NA    NA    NA      34.5         3
      26:  39.3  47.0  39.0    NA    NA    NA    NA    NA    NA      41.8         3
      27:  35.4  36.1  35.2    NA    NA    NA    NA    NA    NA      35.5         3
      28:   2.6  10.1    NA    NA    NA    NA    NA    NA    NA       6.3         2
      29:  49.6  61.1    NA    NA    NA    NA    NA    NA    NA      55.4         2
      30:  49.5  56.9    NA    NA    NA    NA    NA    NA    NA      53.2         2
      31:  62.7  77.4    NA    NA    NA    NA    NA    NA    NA      70.0         2
      32:   8.7    NA    NA    NA    NA    NA    NA    NA    NA       8.7         1
      33:   8.3    NA    NA    NA    NA    NA    NA    NA    NA       8.3         1
      34:  23.1    NA    NA    NA    NA    NA    NA    NA    NA      23.1         1
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count

---

    Code
      round(table_prep(exp_type_2), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>     <num>
       1:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6         9
       2:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7         9
       3:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5         9
       4:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7         9
       5:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2         9
       6:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2         9
       7:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9         9
       8:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3         9
       9:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7         9
      10:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6         9
      11:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5    NA      21.9         8
      12:   4.3   7.9   3.9   4.7   5.4   6.0   3.0    NA    NA       5.0         7
      13:  92.5 100.0  80.8  90.0  73.6  98.7  66.3    NA    NA      86.0         7
      14:  92.6 100.0  84.7 100.0  81.9 100.0  65.4    NA    NA      89.2         7
      15:   2.5  11.6   6.7   8.5   7.4    NA    NA    NA    NA       7.3         5
      16:  16.8  31.5  22.2  23.4  20.5    NA    NA    NA    NA      22.9         5
      17:  22.8  20.7  19.3  31.2  23.7    NA    NA    NA    NA      23.5         5
      18:  35.2  38.4  31.0  35.3  29.1    NA    NA    NA    NA      33.8         5
      19:  38.9  47.7  37.8  42.1  39.5    NA    NA    NA    NA      41.2         5
      20:  45.7  62.0  47.5  57.6  46.6    NA    NA    NA    NA      51.9         5
      21:  45.3  51.2  47.5  51.5  44.9    NA    NA    NA    NA      48.1         5
      22:   2.3  10.7   4.2   8.4    NA    NA    NA    NA    NA       6.4         4
      23:  19.2  25.9  19.0    NA    NA    NA    NA    NA    NA      21.3         3
      24:  30.0  36.0  31.0    NA    NA    NA    NA    NA    NA      32.3         3
      25:  33.7  35.3  34.7    NA    NA    NA    NA    NA    NA      34.5         3
      26:  39.3  47.0  39.0    NA    NA    NA    NA    NA    NA      41.8         3
      27:  35.4  36.1  35.2    NA    NA    NA    NA    NA    NA      35.5         3
      28:   2.6  10.1    NA    NA    NA    NA    NA    NA    NA       6.3         2
      29:  49.6  61.1    NA    NA    NA    NA    NA    NA    NA      55.4         2
      30:  49.5  56.9    NA    NA    NA    NA    NA    NA    NA      53.2         2
      31:  62.7  77.4    NA    NA    NA    NA    NA    NA    NA      70.0         2
      32:   8.7    NA    NA    NA    NA    NA    NA    NA    NA       8.7         1
      33:   8.3    NA    NA    NA    NA    NA    NA    NA    NA       8.3         1
      34:  23.1    NA    NA    NA    NA    NA    NA    NA    NA      23.1         1
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count

---

    Code
      round(table_prep(exp_type_2), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>     <num>
       1:  13.8  18.5   8.5  14.2   7.7  18.0   4.9  15.9  11.6      12.6         9
       2:  15.1  20.5  19.2  17.8  18.1  19.0  11.7  19.6  18.6      17.7         9
       3:  19.1  28.6  16.4  21.6  13.6  21.9  13.1  26.4  14.6      19.5         9
       4:  30.1  36.8  30.1  31.3  31.2  32.2  17.5  48.2  18.8      30.7         9
       5:  25.1  38.7  29.2  38.5  26.2  31.5  19.6  31.2  22.9      29.2         9
       6:  60.6  66.1  60.7  64.5  57.8  60.4  44.3  65.7  61.5      60.2         9
       7:  65.6  74.2  60.4  68.9  67.4  81.1  44.6  76.1  63.9      66.9         9
       8:  62.8  81.4  63.0  69.8  63.4  74.4  48.4  69.3  72.6      67.3         9
       9:  66.7  78.7  61.2  76.6  67.7  77.5  44.5  76.8  68.1      68.7         9
      10:  94.2 100.0  88.3  94.4  88.3  92.6  75.5 100.0 100.0      92.6         9
      11:  16.5  19.6  19.8  23.9  19.1  24.4  15.1  36.5    NA      21.9         8
      12:   4.3   7.9   3.9   4.7   5.4   6.0   3.0    NA    NA       5.0         7
      13:  92.5 100.0  80.8  90.0  73.6  98.7  66.3    NA    NA      86.0         7
      14:  92.6 100.0  84.7 100.0  81.9 100.0  65.4    NA    NA      89.2         7
      15:   2.5  11.6   6.7   8.5   7.4    NA    NA    NA    NA       7.3         5
      16:  16.8  31.5  22.2  23.4  20.5    NA    NA    NA    NA      22.9         5
      17:  22.8  20.7  19.3  31.2  23.7    NA    NA    NA    NA      23.5         5
      18:  35.2  38.4  31.0  35.3  29.1    NA    NA    NA    NA      33.8         5
      19:  38.9  47.7  37.8  42.1  39.5    NA    NA    NA    NA      41.2         5
      20:  45.7  62.0  47.5  57.6  46.6    NA    NA    NA    NA      51.9         5
      21:  45.3  51.2  47.5  51.5  44.9    NA    NA    NA    NA      48.1         5
      22:   2.3  10.7   4.2   8.4    NA    NA    NA    NA    NA       6.4         4
      23:  19.2  25.9  19.0    NA    NA    NA    NA    NA    NA      21.3         3
      24:  30.0  36.0  31.0    NA    NA    NA    NA    NA    NA      32.3         3
      25:  33.7  35.3  34.7    NA    NA    NA    NA    NA    NA      34.5         3
      26:  39.3  47.0  39.0    NA    NA    NA    NA    NA    NA      41.8         3
      27:  35.4  36.1  35.2    NA    NA    NA    NA    NA    NA      35.5         3
      28:   2.6  10.1    NA    NA    NA    NA    NA    NA    NA       6.3         2
      29:  49.6  61.1    NA    NA    NA    NA    NA    NA    NA      55.4         2
      30:  49.5  56.9    NA    NA    NA    NA    NA    NA    NA      53.2         2
      31:  62.7  77.4    NA    NA    NA    NA    NA    NA    NA      70.0         2
      32:   8.7    NA    NA    NA    NA    NA    NA    NA    NA       8.7         1
      33:   8.3    NA    NA    NA    NA    NA    NA    NA    NA       8.3         1
      34:  23.1    NA    NA    NA    NA    NA    NA    NA    NA      23.1         1
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count

---

    Code
      round(table_prep(cal_type_2), 1)
    Output
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count
          <num> <num> <num> <num> <num> <num> <num> <num> <num>     <num>     <num>
       1:  38.9  47.7  37.8  42.1  39.5  43.7  31.9    NA    NA      40.2         7
       2:  33.7  35.3  34.7  36.8  31.4  35.6  27.2    NA    NA      33.5         7
       3:  39.3  47.0  39.0  41.9  41.6  43.7  22.8    NA    NA      39.3         7
       4:  39.3  47.0  39.0  41.9  41.6  43.7  22.8    NA    NA      39.3         7
       5:  35.4  36.1  35.2  38.4  31.0  35.3  29.1    NA    NA      34.3         7
       6:  35.4  36.1  35.2  38.4  31.0  35.3  29.1  38.7  36.4      35.1         9
       7:  39.3  47.0  39.0  41.9  41.6  43.7  22.8  66.7  23.8      40.7         9
       8:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6         9
       9:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6         9
      10:  35.4  36.1  35.2  38.4  31.0    NA    NA    NA    NA      35.2         5
      11:  39.3  47.0  39.0  41.9  41.6    NA    NA    NA    NA      41.8         5
      12:  33.7  35.3  34.7  36.8  31.4    NA    NA    NA    NA      34.4         5
      13:  35.4  36.1  35.2  38.4  31.0    NA    NA    NA    NA      35.2         5
      14:  38.9  47.7    NA    NA    NA    NA    NA    NA    NA      43.3         2
      15:  39.3  47.0    NA    NA    NA    NA    NA    NA    NA      43.2         2
      16:  35.4  36.1    NA    NA    NA    NA    NA    NA    NA      35.7         2
      17:  38.9    NA    NA    NA    NA    NA    NA    NA    NA      38.9         1
      18:  42.1    NA    NA    NA    NA    NA    NA    NA    NA      42.1         1
      19:  38.4    NA    NA    NA    NA    NA    NA    NA    NA      38.4         1
      20:  42.1    NA    NA    NA    NA    NA    NA    NA    NA      42.1         1
      21:  38.4    NA    NA    NA    NA    NA    NA    NA    NA      38.4         1
      22:  41.9    NA    NA    NA    NA    NA    NA    NA    NA      41.9         1
      23:  39.3  47.0  39.0  41.9  41.6  43.7  22.8  66.7    NA      42.8         8
      24:  33.7  35.3  34.7  36.8  31.4  35.6  27.2  41.0    NA      34.5         8
      25:  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5    NA      40.8         8
      26:  37.8  42.1  39.5    NA    NA    NA    NA    NA    NA      39.8         3
      27:  35.2  38.4  31.0    NA    NA    NA    NA    NA    NA      34.8         3
      28:  35.2  38.4  31.0    NA    NA    NA    NA    NA    NA      34.8         3
      29:  39.0  41.9  41.6    NA    NA    NA    NA    NA    NA      40.8         3
      30:  34.7  36.8  31.4    NA    NA    NA    NA    NA    NA      34.3         3
      31:  38.9  47.7  37.8  42.1  39.5  43.7  31.9  42.4  41.5      40.6         9
      32:  33.7  35.3  34.7  36.8  37.8  42.1  39.5  43.7  31.9      37.3         9
      33:  35.4  36.1  35.2  38.4    NA    NA    NA    NA    NA      36.3         4
          CpG#1 CpG#2 CpG#3 CpG#4 CpG#5 CpG#6 CpG#7 CpG#8 CpG#9 row_means CpG_count

