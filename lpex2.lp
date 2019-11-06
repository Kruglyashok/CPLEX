\ENCODING=ISO-8859-1
\Problem name: ilog.cplex

Maximize
 obj: - x1 + x2
Subject To
 c1: - x1 + x3 >= -40
 c2: - x3 + x4 >= -22
 c3: - x3 + x5 >= -6
 c4: - x4 + x5 >= -14
 c5: - x1 + x6 >= -180
 c6: - x3 + x6 >= -45
 c7: - x5 + x6 >= -20
Bounds
      x1 >= 0
      x2 >= 0
      x3 >= 0
      x4 >= 0
      x5 >= 0
      x6 >= 0
Generals
 x1  x2  x3  x4  x5  x6 
End
