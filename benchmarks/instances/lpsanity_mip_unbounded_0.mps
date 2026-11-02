NAME mip_unbounded
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 L  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    x  c2  -1
    y  obj  1  c1  -1
    y  c2  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  2.5  c2  2.5
BOUNDS
 PL bnd  x
 PL bnd  y
ENDATA
