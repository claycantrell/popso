NAME mip_marker_defaults
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    x  obj  1  c1  1
    y  obj  1  c1  1
    z  obj  -1  c2  1
    MARKER  'MARKER'  'INTEND'
RHS
    rhs  c1  7.5  c2  -3.5
BOUNDS
 LO bnd  y  2
 MI bnd  z
ENDATA
