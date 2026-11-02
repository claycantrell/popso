NAME upper_only
OBJSENSE
    MAX
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    x  obj  1  c1  1
    y  obj  1  c1  1
    z  obj  -1  c2  1
RHS
    rhs  c1  100  c2  -4
BOUNDS
 UP bnd  x  5
 MI bnd  y
 UP bnd  y  7
 MI bnd  z
ENDATA
