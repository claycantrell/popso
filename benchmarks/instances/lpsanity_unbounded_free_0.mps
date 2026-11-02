* min x with x free, x - y <= 3, y in [0, 10]
NAME unbounded_free
ROWS
 N  obj
 L  c1
COLUMNS
    x  obj  1  c1  1
    y  c1  -1
RHS
    rhs  c1  3
BOUNDS
 FR bnd  x
 UP bnd  y  10
ENDATA
