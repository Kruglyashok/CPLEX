\ENCODING=ISO-8859-1
\Problem name: ilog.cplex

Minimize
 obj: 40 x1 + 22 x2 + 6 x3 + 14 x4 + 180 x5 + 45 x6 + 20 x7
Subject To
 c1: x1 + x5  = 1
 c2: x5 + x6 + x7  = 1
 c3: x1 - x2 - x3 - x6  = 0
 c4: x2 - x4  = 0
 c5: x3 + x4 - x7  = 0
Bounds
      x1 >= 0
      x2 >= 0
      x3 >= 0
      x4 >= 0
      x5 >= 0
      x6 >= 0
      x7 >= 0
Generals
 x1  x2  x3  x4  x5  x6  x7 
End
