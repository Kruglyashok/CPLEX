\ENCODING=ISO-8859-1
\Problem name: ilog.cplex

Minimize
 obj: 40 x1 + 22 x2 + 6 x3 + 14 x4 + 180 x5 + 45 x6 + 20 x7 + 7 x8 + x9 + x10
      + x11 + x12 + x13 + x14 + x15
Subject To
 c1:  x1 + x5  = 1
 c2:  x5 + x6 + x7  = 1
 c3:  x1 - x2 - x3 - x6  = 0
 c4:  x2 - x4  = 0
 c5:  x3 + x4 - x7  = 0
 c6:  - x1 + x8 + x9 >= 0
 c7:  - x2 + x8 + x10 >= 0
 c8:  - x3 + x8 + x11 >= 0
 c9:  - x4 + x8 + x12 >= 0
 c10: - x5 + x8 + x13 >= 0
 c11: - x6 + x8 + x14 >= 0
 c12: - x7 + x8 + x15 >= 0
Bounds
 0 <= x1 <= 1
 0 <= x2 <= 1
 0 <= x3 <= 1
 0 <= x4 <= 1
 0 <= x5 <= 1
 0 <= x6 <= 1
 0 <= x7 <= 1
      x8 >= 0
      x9 >= 0
      x10 >= 0
      x11 >= 0
      x12 >= 0
      x13 >= 0
      x14 >= 0
      x15 >= 0
Generals
 x1  x2  x3  x4  x5  x6  x7  x8  x9  x10  x11  x12  x13  x14  x15 
End
