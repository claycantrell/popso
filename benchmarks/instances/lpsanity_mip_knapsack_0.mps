NAME mip_knapsack
OBJSENSE MAX
ROWS
 N  obj
 L  cap
COLUMNS
    a  obj  10  cap  4
    b  obj  13  cap  6
    c  obj  7  cap  3
    d  obj  8  cap  4
    e  obj  4  cap  2
RHS
    rhs  cap  11
BOUNDS
 BV bnd  a
 BV bnd  b
 BV bnd  c
 BV bnd  d
 BV bnd  e
ENDATA
