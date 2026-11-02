NAME jobshop_tiny_1
ROWS
 N  OBJ
 G  prec0_0 
 G  prec0_1 
 G  mk0     
 G  prec1_0 
 G  prec1_1 
 G  mk1     
 G  prec2_0 
 G  prec2_1 
 G  mk2     
 L  d0_1_0a 
 L  d0_1_0b 
 L  d0_2_0a 
 L  d0_2_0b 
 L  d1_2_0a 
 L  d1_2_0b 
 L  d0_1_1a 
 L  d0_1_1b 
 L  d0_2_1a 
 L  d0_2_1b 
 L  d1_2_1a 
 L  d1_2_1b 
 L  d0_1_2a 
 L  d0_1_2b 
 L  d0_2_2a 
 L  d0_2_2b 
 L  d1_2_2a 
 L  d1_2_2b 
COLUMNS
    s0_0      prec0_0   1
    s0_0      prec0_1   -1
    s0_0      d0_1_0a   1
    s0_0      d0_1_0b   -1
    s0_0      d0_2_0a   1
    s0_0      d0_2_0b   -1
    s0_1      prec0_1   1
    s0_1      mk0       -1
    s0_1      d0_1_1a   1
    s0_1      d0_1_1b   -1
    s0_1      d0_2_1a   1
    s0_1      d0_2_1b   -1
    s0_2      prec0_0   -1
    s0_2      d0_1_2a   1
    s0_2      d0_1_2b   -1
    s0_2      d0_2_2a   1
    s0_2      d0_2_2b   -1
    s1_0      prec1_0   -1
    s1_0      d0_1_0a   -1
    s1_0      d0_1_0b   1
    s1_0      d1_2_0a   1
    s1_0      d1_2_0b   -1
    s1_1      prec1_1   1
    s1_1      mk1       -1
    s1_1      d0_1_1a   -1
    s1_1      d0_1_1b   1
    s1_1      d1_2_1a   1
    s1_1      d1_2_1b   -1
    s1_2      prec1_0   1
    s1_2      prec1_1   -1
    s1_2      d0_1_2a   -1
    s1_2      d0_1_2b   1
    s1_2      d1_2_2a   1
    s1_2      d1_2_2b   -1
    s2_0      prec2_0   -1
    s2_0      d0_2_0a   -1
    s2_0      d0_2_0b   1
    s2_0      d1_2_0a   -1
    s2_0      d1_2_0b   1
    s2_1      prec2_0   1
    s2_1      prec2_1   -1
    s2_1      d0_2_1a   -1
    s2_1      d0_2_1b   1
    s2_1      d1_2_1a   -1
    s2_1      d1_2_1b   1
    s2_2      prec2_1   1
    s2_2      mk2       -1
    s2_2      d0_2_2a   -1
    s2_2      d0_2_2b   1
    s2_2      d1_2_2a   -1
    s2_2      d1_2_2b   1
    cmax      OBJ       1
    cmax      mk0       1
    cmax      mk1       1
    cmax      mk2       1
    MARKER    'MARKER'                 'INTORG'
    z0_1_0    d0_1_0a   565
    z0_1_0    d0_1_0b   -565
    z0_2_0    d0_2_0a   565
    z0_2_0    d0_2_0b   -565
    z1_2_0    d1_2_0a   565
    z1_2_0    d1_2_0b   -565
    z0_1_1    d0_1_1a   565
    z0_1_1    d0_1_1b   -565
    z0_2_1    d0_2_1a   565
    z0_2_1    d0_2_1b   -565
    z1_2_1    d1_2_1a   565
    z1_2_1    d1_2_1b   -565
    z0_1_2    d0_1_2a   565
    z0_1_2    d0_1_2b   -565
    z0_2_2    d0_2_2a   565
    z0_2_2    d0_2_2b   -565
    z1_2_2    d1_2_2a   565
    z1_2_2    d1_2_2b   -565
    MARKER    'MARKER'                 'INTEND'
RHS
    RHS1      prec0_0   48
    RHS1      prec0_1   69
    RHS1      mk0       51
    RHS1      prec1_0   29
    RHS1      prec1_1   85
    RHS1      mk1       20
    RHS1      prec2_0   86
    RHS1      prec2_1   82
    RHS1      mk2       95
    RHS1      d0_1_0a   496
    RHS1      d0_1_0b   -29
    RHS1      d0_2_0a   496
    RHS1      d0_2_0b   -86
    RHS1      d1_2_0a   536
    RHS1      d1_2_0b   -86
    RHS1      d0_1_1a   514
    RHS1      d0_1_1b   -20
    RHS1      d0_2_1a   514
    RHS1      d0_2_1b   -82
    RHS1      d1_2_1a   545
    RHS1      d1_2_1b   -82
    RHS1      d0_1_2a   517
    RHS1      d0_1_2b   -85
    RHS1      d0_2_2a   517
    RHS1      d0_2_2b   -95
    RHS1      d1_2_2a   480
    RHS1      d1_2_2b   -95
BOUNDS
 UP BND1      s0_0      565
 UP BND1      s0_1      565
 UP BND1      s0_2      565
 UP BND1      s1_0      565
 UP BND1      s1_1      565
 UP BND1      s1_2      565
 UP BND1      s2_0      565
 UP BND1      s2_1      565
 UP BND1      s2_2      565
 UP BND1      cmax      565
 BV BND1      z0_1_0  
 BV BND1      z0_2_0  
 BV BND1      z1_2_0  
 BV BND1      z0_1_1  
 BV BND1      z0_2_1  
 BV BND1      z1_2_1  
 BV BND1      z0_1_2  
 BV BND1      z0_2_2  
 BV BND1      z1_2_2  
ENDATA
