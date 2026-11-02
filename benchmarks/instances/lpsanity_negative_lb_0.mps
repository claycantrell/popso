NAME negative_lb
ROWS
 N  obj
 G  c1
 L  c2
COLUMNS
    x  obj  1  c1  1
    x  c2  1
    y  obj  2  c1  1
    y  c2  -1
RHS
    rhs  c1  -8  c2  4
BOUNDS
 LO bnd  x  -3
 UP bnd  x  5
 LO bnd  y  -10
ENDATA
