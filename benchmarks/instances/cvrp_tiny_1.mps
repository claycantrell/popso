NAME cvrp_tiny_1
ROWS
 N  OBJ
 E  out1    
 E  in1     
 E  out2    
 E  in2     
 E  out3    
 E  in3     
 E  out4    
 E  in4     
 E  out5    
 E  in5     
 E  depot   
 L  load1_2 
 L  load1_3 
 L  load1_4 
 L  load1_5 
 L  load2_1 
 L  load2_3 
 L  load2_4 
 L  load2_5 
 L  load3_1 
 L  load3_2 
 L  load3_4 
 L  load3_5 
 L  load4_1 
 L  load4_2 
 L  load4_3 
 L  load4_5 
 L  load5_1 
 L  load5_2 
 L  load5_3 
 L  load5_4 
COLUMNS
    MARKER    'MARKER'                 'INTORG'
    x0_1      OBJ       60
    x0_1      in1       1
    x0_1      depot     1
    x0_2      OBJ       64
    x0_2      in2       1
    x0_2      depot     1
    x0_3      OBJ       42
    x0_3      in3       1
    x0_3      depot     1
    x0_4      OBJ       72
    x0_4      in4       1
    x0_4      depot     1
    x0_5      OBJ       23
    x0_5      in5       1
    x0_5      depot     1
    x1_0      OBJ       60
    x1_0      out1      1
    x1_0      depot     -1
    x1_2      OBJ       87
    x1_2      out1      1
    x1_2      in2       1
    x1_2      load1_2   16
    x1_3      OBJ       46
    x1_3      out1      1
    x1_3      in3       1
    x1_3      load1_3   16
    x1_4      OBJ       83
    x1_4      out1      1
    x1_4      in4       1
    x1_4      load1_4   16
    x1_5      OBJ       65
    x1_5      out1      1
    x1_5      in5       1
    x1_5      load1_5   16
    x2_0      OBJ       64
    x2_0      out2      1
    x2_0      depot     -1
    x2_1      OBJ       87
    x2_1      in1       1
    x2_1      out2      1
    x2_1      load2_1   16
    x2_3      OBJ       99
    x2_3      out2      1
    x2_3      in3       1
    x2_3      load2_3   16
    x2_4      OBJ       16
    x2_4      out2      1
    x2_4      in4       1
    x2_4      load2_4   16
    x2_5      OBJ       41
    x2_5      out2      1
    x2_5      in5       1
    x2_5      load2_5   16
    x3_0      OBJ       42
    x3_0      out3      1
    x3_0      depot     -1
    x3_1      OBJ       46
    x3_1      in1       1
    x3_1      out3      1
    x3_1      load3_1   16
    x3_2      OBJ       99
    x3_2      in2       1
    x3_2      out3      1
    x3_2      load3_2   16
    x3_4      OBJ       103
    x3_4      out3      1
    x3_4      in4       1
    x3_4      load3_4   16
    x3_5      OBJ       62
    x3_5      out3      1
    x3_5      in5       1
    x3_5      load3_5   16
    x4_0      OBJ       72
    x4_0      out4      1
    x4_0      depot     -1
    x4_1      OBJ       83
    x4_1      in1       1
    x4_1      out4      1
    x4_1      load4_1   16
    x4_2      OBJ       16
    x4_2      in2       1
    x4_2      out4      1
    x4_2      load4_2   16
    x4_3      OBJ       103
    x4_3      in3       1
    x4_3      out4      1
    x4_3      load4_3   16
    x4_5      OBJ       50
    x4_5      out4      1
    x4_5      in5       1
    x4_5      load4_5   16
    x5_0      OBJ       23
    x5_0      out5      1
    x5_0      depot     -1
    x5_1      OBJ       65
    x5_1      in1       1
    x5_1      out5      1
    x5_1      load5_1   16
    x5_2      OBJ       41
    x5_2      in2       1
    x5_2      out5      1
    x5_2      load5_2   16
    x5_3      OBJ       62
    x5_3      in3       1
    x5_3      out5      1
    x5_3      load5_3   16
    x5_4      OBJ       50
    x5_4      in4       1
    x5_4      out5      1
    x5_4      load5_4   16
    MARKER    'MARKER'                 'INTEND'
    u1        load1_2   1
    u1        load1_3   1
    u1        load1_4   1
    u1        load1_5   1
    u1        load2_1   -1
    u1        load3_1   -1
    u1        load4_1   -1
    u1        load5_1   -1
    u2        load1_2   -1
    u2        load2_1   1
    u2        load2_3   1
    u2        load2_4   1
    u2        load2_5   1
    u2        load3_2   -1
    u2        load4_2   -1
    u2        load5_2   -1
    u3        load1_3   -1
    u3        load2_3   -1
    u3        load3_1   1
    u3        load3_2   1
    u3        load3_4   1
    u3        load3_5   1
    u3        load4_3   -1
    u3        load5_3   -1
    u4        load1_4   -1
    u4        load2_4   -1
    u4        load3_4   -1
    u4        load4_1   1
    u4        load4_2   1
    u4        load4_3   1
    u4        load4_5   1
    u4        load5_4   -1
    u5        load1_5   -1
    u5        load2_5   -1
    u5        load3_5   -1
    u5        load4_5   -1
    u5        load5_1   1
    u5        load5_2   1
    u5        load5_3   1
    u5        load5_4   1
RHS
    RHS1      out1      1
    RHS1      in1       1
    RHS1      out2      1
    RHS1      in2       1
    RHS1      out3      1
    RHS1      in3       1
    RHS1      out4      1
    RHS1      in4       1
    RHS1      out5      1
    RHS1      in5       1
    RHS1      load1_2   12
    RHS1      load1_3   7
    RHS1      load1_4   9
    RHS1      load1_5   13
    RHS1      load2_1   8
    RHS1      load2_3   7
    RHS1      load2_4   9
    RHS1      load2_5   13
    RHS1      load3_1   8
    RHS1      load3_2   12
    RHS1      load3_4   9
    RHS1      load3_5   13
    RHS1      load4_1   8
    RHS1      load4_2   12
    RHS1      load4_3   7
    RHS1      load4_5   13
    RHS1      load5_1   8
    RHS1      load5_2   12
    RHS1      load5_3   7
    RHS1      load5_4   9
BOUNDS
 BV BND1      x0_1    
 BV BND1      x0_2    
 BV BND1      x0_3    
 BV BND1      x0_4    
 BV BND1      x0_5    
 BV BND1      x1_0    
 BV BND1      x1_2    
 BV BND1      x1_3    
 BV BND1      x1_4    
 BV BND1      x1_5    
 BV BND1      x2_0    
 BV BND1      x2_1    
 BV BND1      x2_3    
 BV BND1      x2_4    
 BV BND1      x2_5    
 BV BND1      x3_0    
 BV BND1      x3_1    
 BV BND1      x3_2    
 BV BND1      x3_4    
 BV BND1      x3_5    
 BV BND1      x4_0    
 BV BND1      x4_1    
 BV BND1      x4_2    
 BV BND1      x4_3    
 BV BND1      x4_5    
 BV BND1      x5_0    
 BV BND1      x5_1    
 BV BND1      x5_2    
 BV BND1      x5_3    
 BV BND1      x5_4    
 LO BND1      u1        8
 UP BND1      u1        16
 LO BND1      u2        4
 UP BND1      u2        16
 LO BND1      u3        9
 UP BND1      u3        16
 LO BND1      u4        7
 UP BND1      u4        16
 LO BND1      u5        3
 UP BND1      u5        16
ENDATA
