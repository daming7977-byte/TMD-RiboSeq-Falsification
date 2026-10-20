"""Recalculate only ancillary summaries from saved gene-deletion effects.

Writes to a required separate output directory. No bootstrap, refit or raw reads.
Original author-owned implementation under the package's documented MIT scope.
"""
from pathlib import Path
import argparse,csv,json,hashlib,math
from nargab_region_composite.influence_summary_policy import influence_summary
ROOT=Path(__file__).resolve().parent
FIELDS=['spearman_total_units_vs_absolute_change','top_TMD_count_decile_share_of_total_absolute_change']
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output-dir',type=Path,required=True)
    output=parser.parse_args().output_dir.resolve()
    if output==ROOT or output.is_relative_to(ROOT/'reproduction_output'):
        raise ValueError('Use a separate output directory, not saved reproduction outputs')
    output.mkdir(parents=True,exist_ok=True);record={}
    for mode in ['available_case','complete_case']:
        source=ROOT/'reproduction_output'/mode;dest=output/mode;dest.mkdir(exist_ok=True)
        effects=read(source/'REGION_COMPOSITE_GENE_INFLUENCE_v1.tsv')
        summaries=read(source/'REGION_COMPOSITE_INFLUENCE_SUMMARY_v1.tsv');details=[];members=[]
        for row in summaries:
            subset=[r for r in effects if r['region']==row['region']];result=influence_summary(subset)
            for key in FIELDS:row[key]='NA' if not math.isfinite(result[key]) else repr(result[key])
            details.append({'region':row['region'],**{k:v for k,v in result.items() if k!='selected_gene_ids'}})
            lookup={r['gene_id']:r for r in subset}
            members.extend({'region':row['region'],'gene_id':g,'total_units':lookup[g]['total_units']} for g in result['selected_gene_ids'])
        write(dest/'REGION_COMPOSITE_INFLUENCE_SUMMARY_v1.tsv',summaries)
        write(dest/'INFLUENCE_TIE_DETAILS_PHASE2.tsv',details);write(dest/'TOP_COUNT_DECILE_MEMBERS_PHASE2.tsv',members)
        record[mode]={'source_gene_effect_sha256':hashlib.sha256((source/'REGION_COMPOSITE_GENE_INFLUENCE_v1.tsv').read_bytes()).hexdigest(),'summary_sha256':hashlib.sha256((dest/'REGION_COMPOSITE_INFLUENCE_SUMMARY_v1.tsv').read_bytes()).hexdigest()}
    (output/'SUMMARY_ONLY_REGENERATION.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))
if __name__=='__main__':main()
