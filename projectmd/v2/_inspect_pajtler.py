"""Quick inspection: GSE64415 series matrix subtype × primary/relapse breakdown."""
import gzip, re
from collections import Counter

path = 'output/external_cache/GSE64415_series_matrix.txt.gz'
fields = {}
with gzip.open(path, 'rt') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!Sample_geo_accession'):
            fields['gsm'] = line.strip().split('\t')[1:]
        elif line.startswith('!Sample_title'):
            fields['title'] = line.strip().split('\t')[1:]
        elif line.startswith('!Sample_characteristics_ch1'):
            parts = line.strip().split('\t')[1:]
            sample = next((p for p in parts if p and p != '""'), '')
            m = re.match(r'"([^:]+):', sample)
            if m:
                key = m.group(1).strip()
                fields[key] = [re.sub(r'^"[^:]+:\s*', '', p).rstrip('"') for p in parts]

n = len(fields.get('gsm', []))
print(f'Total samples: {n}')
print(f'Fields: {list(fields.keys())}')
print()

sub = fields.get('molecular subgroup', [])
pr = fields.get('primary/relapse', [])
print('Subgroup counts:')
for k, v in sorted(Counter(sub).items(), key=lambda x: -x[1]):
    print(f'  {k}: {v}')
print('\nPrimary/relapse counts:')
for k, v in sorted(Counter(pr).items(), key=lambda x: -x[1]):
    print(f'  {k}: {v}')
print('\nSubgroup × primary/relapse:')
for (s, p), v in sorted(Counter(zip(sub, pr)).items()):
    print(f'  {s:20s} | {p:10s} : {v}')
