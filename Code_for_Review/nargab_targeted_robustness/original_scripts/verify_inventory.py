from pathlib import Path
import csv,hashlib,json
P=Path(__file__).resolve().parent
with (P/'P3_OVERLAP_INVENTORY.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
with (P/'pair_overlap_flags.tsv').open() as f:flags=list(csv.DictReader(f,delimiter='\t'))
seen={};cases=[]
for r in rows:
 assert int(r['unmasked_nt'])==r['mask_bits'].count('0')
 assert 0<=int(r['unmasked_other_gene_overlap_nt'])<=int(r['full_other_gene_overlap_nt'])<=31
 assert int(r['unmasked_other_gene_overlap_nt'])<=int(r['unmasked_nt'])
 assert not (set(r['owner_gene_ids'].split(','))&set(r['other_gene_ids_full'].split(',')))
 assert r['continuous_target_cds_stop_path_valid']=='YES'
 key=(r['paired_real_unit_id'],r['side'],r['pseudo_index'])
 val=(r['geometry'],r['mask_bits'],r['full_other_gene_overlap_nt'],r['unmasked_other_gene_overlap_nt'])
 assert key not in seen or seen[key]==val
 seen[key]=val
 if r['owner_gene_ids']=='ENSG00000251357':cases.append(r)
flagged={k[0] for k,v in seen.items() if int(v[3])>0}
assert flagged=={r['paired_real_unit_id'] for r in flags if r['exclude_pair_unmasked_other_gene']=='YES'}
sh=hashlib.sha256((P/'pair_overlap_flags.tsv').read_bytes()).hexdigest()
assert sh=='e1c050928519bff91a3ad1c449dd3a334641f3e0fa8e280a8d4fc086a53c19a0'
with (P/'ENSG00000251357_exact_window_evidence.tsv').open('w') as f:
 w=csv.DictWriter(f,delimiter='\t',fieldnames=list(cases[0]));w.writeheader();w.writerows(cases)
out={'all_gates':'PASS','analysis_window_rows':len(rows),'unique_pair_window_components':len(seen),'pair_rows':len(flags),'positive_pairs_union':len(flagged),'flags_unchanged_from_pre_descriptive_label_refinement_sha256':sh,'invalid_target_paths':0,'target_case_records':len(cases)}
(P/'inventory_self_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
