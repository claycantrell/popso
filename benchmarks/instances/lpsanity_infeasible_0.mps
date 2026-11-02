* x + y <= 1 and x + y >= 2
NAME infeasible
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    x  obj  1  c1  1
    x  c2  1
    y  obj  1  c1  1
    y  c2  1
RHS
    rhs  c1  1  c2  2
ENDATA
