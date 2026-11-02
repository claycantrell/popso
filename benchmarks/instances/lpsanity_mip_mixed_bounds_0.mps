NAME mip_mixed_bounds
ROWS
 N  obj
 G  c1
 G  c2
 L  c3
 L  c4
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    x  c2  1  c3  1
    x  c4  1
    y  obj  -1  c1  1
    y  c2  -1
    v  c4  1
    MARKER  'MARKER'  'INTEND'
    w  obj  2  c2  1
RHS
    rhs  obj  -100
    rhs  c1  -3.5  c2  0.5
    rhs  c3  4  c4  6
BOUNDS
 FR bnd  x
 LO bnd  y  -5
 UP bnd  y  5
 FX bnd  v  3
ENDATA
