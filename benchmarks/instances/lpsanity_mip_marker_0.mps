NAME mip_marker
OBJSENSE MAX
ROWS
 N  obj
 L  c1
 G  c2
COLUMNS
    MARKER  'MARKER'  'INTORG'
    a  obj  5  c1  2
    a  c2  1
    b  obj  4  c1  3
    b  c2  1
    MARKER  'MARKER'  'INTEND'
    c  obj  3  c1  1
    c  c2  -1
    d  obj  1  c1  1
RHS
    rhs  c1  7.5  c2  1
BOUNDS
 BV bnd  a
 UP bnd  b  3
 LI bnd  c  -2
 UI bnd  c  2
 UP bnd  d  1.5
ENDATA
