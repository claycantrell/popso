import highspy, numpy as np, sys
inf = float('inf')
toks = open(sys.argv[1]).read().split('\n')
k = 0; mism = 0; tot = 0
names = {0:'optimal',1:'infeasible',2:'unbounded'}
while k < len(toks) and toks[k].startswith('LP'):
    _, n, m, st, obj = toks[k].split(); n, m, st, obj = int(n), int(m), int(st), float(obj); k += 1
    c=[];cl=[];cu=[]
    for j in range(n): a,b,d = map(float, toks[k].split()); c.append(a); cl.append(b); cu.append(d); k+=1
    rl=[];ru=[]
    for i in range(m): a,b = map(float, toks[k].split()); rl.append(a); ru.append(b); k+=1
    nnz = int(toks[k]); k+=1
    h = highspy.Highs(); h.setOptionValue('output_flag', False); h.setOptionValue('presolve','off')
    lp = highspy.HighsLp(); lp.num_col_=n; lp.num_row_=m
    lp.col_cost_=np.array(c); lp.col_lower_=np.array(cl); lp.col_upper_=np.array(cu); lp.row_lower_=np.array(rl); lp.row_upper_=np.array(ru)
    cols=[[] for _ in range(n)]
    for _ in range(nnz): j,i,v = toks[k].split(); cols[int(j)].append((int(i),float(v))); k+=1
    start=[0]; idx=[]; val=[]
    for j in range(n):
        for i,v in cols[j]: idx.append(i); val.append(v)
        start.append(len(idx))
    lp.a_matrix_.format_=highspy.MatrixFormat.kColwise; lp.a_matrix_.start_=np.array(start); lp.a_matrix_.index_=np.array(idx,dtype=np.int32); lp.a_matrix_.value_=np.array(val); lp.a_matrix_.num_col_=n; lp.a_matrix_.num_row_=m
    h.passModel(lp); h.run(); ms = h.modelStatusToString(h.getModelStatus()).lower()
    hs = 'optimal' if ms=='optimal' else 'infeasible' if 'infeasible' in ms and 'unbounded' not in ms else 'unbounded' if ms=='unbounded' else ms
    tot += 1
    ours = names.get(st, str(st))
    ok = ours == hs or (hs.startswith('primal infeasible or unbounded') and ours in ('infeasible','unbounded'))
    if ok and ours == 'optimal' and abs(obj - h.getInfo().objective_function_value) > 1e-6*(1+abs(obj)): ok = False
    if not ok:
        mism += 1
        if mism < 10: print('mismatch LP', tot-1, 'ours', ours, obj, 'highs', ms, h.getInfo().objective_function_value)
print('checked', tot, 'mismatches', mism)
