* x <= 1, y <= 1 (bounds) but x + y >= 3
NAME infeasible_bounds
ROWS
 N  obj
 G  c1
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
RHS
    rhs  c1  3
BOUNDS
 UP bnd  x  1
 UP bnd  y  1
ENDATA
