"""Phase 2 deterministic ancillary summaries; no change to effect estimation.

Original author-owned implementation under the package's documented MIT scope.
Policy: decimal half-even at 1e-12, average ranks, count descending / gene ID
ascending, exactly ceiling(G/10) selected genes. No stochastic tie-breaking.
"""
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
import math

def average_ranks(values):
    order=sorted(range(len(values)),key=lambda i:values[i])
    ranks=[0.0]*len(values);i=0
    while i<len(order):
        j=i+1
        while j<len(order) and values[order[j]]==values[order[i]]:j+=1
        for k in order[i:j]:ranks[k]=(i+1+j)/2
        i=j
    return ranks

def influence_summary(records):
    records=sorted(records,key=lambda r:str(r['gene_id']))
    if not records or len({str(r['gene_id']) for r in records})!=len(records):
        raise ValueError('Expected nonempty, uniquely identified gene records')
    with localcontext() as ctx:
        ctx.prec=50
        values=[abs(Decimal(str(r['signed_change'])).quantize(Decimal('0.000000000001'),rounding=ROUND_HALF_EVEN)) for r in records]
        if not all(v.is_finite() for v in values):raise ValueError('Nonfinite influence')
        counts=[int(r['total_units']) for r in records]
        x,y=average_ranks(counts),average_ranks(values)
        mx,my=math.fsum(x)/len(x),math.fsum(y)/len(y)
        den=math.sqrt(math.fsum((v-mx)**2 for v in x)*math.fsum((v-my)**2 for v in y))
        rho=math.fsum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else float('nan')
        k=(len(records)+9)//10
        chosen=sorted(range(len(records)),key=lambda i:(-counts[i],str(records[i]['gene_id'])))[:k]
        cutoff=counts[chosen[-1]]
        total=sum(values,Decimal(0))
        share=float(sum((values[i] for i in chosen),Decimal(0))/total) if total else float('nan')
        above=[i for i,n in enumerate(counts) if n>cutoff]
        tied=sorted(values[i] for i,n in enumerate(counts) if n==cutoff)
        slots=k-len(above);fixed=sum((values[i] for i in above),Decimal(0))
        lo=float((fixed+sum(tied[:slots],Decimal(0)))/total) if total else float('nan')
        hi=float((fixed+sum(tied[-slots:],Decimal(0)))/total) if total else float('nan')
    return {'spearman_total_units_vs_absolute_change':rho,
            'top_TMD_count_decile_share_of_total_absolute_change':share,
            'selected_gene_ids':[str(records[i]['gene_id']) for i in chosen],
            'K':k,'count_cutoff':cutoff,'cutoff_tied_genes':len(tied),'slots_from_cutoff_tie':slots,
            'tie_share_min':lo,'tie_share_max':hi,'rounded_zero_n':sum(v==0 for v in values)}
