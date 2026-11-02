NAME objsense_maximize
OBJSENSE
    MAXIMIZE
ROWS
 N  profit
 L  c1
 L  c2
COLUMNS
    x  profit  3  c1  1
    x  c2  1
    y  profit  2  c1  1
    y  c2  3
RHS
    rhs  c1  4  c2  6
BOUNDS
 UP bnd  x  3
ENDATA
