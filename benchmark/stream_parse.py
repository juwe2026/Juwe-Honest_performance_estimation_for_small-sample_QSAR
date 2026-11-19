import sys, re, csv, sqlite3, io, time
KEEP = {
 'activities': ['activity_id','assay_id','doc_id','molregno','standard_relation','standard_value','standard_units','standard_type','pchembl_value','data_validity_comment','potential_duplicate','activity_comment','src_id','text_value','standard_text_value'],
 'assays': ['assay_id','doc_id','description','assay_type','assay_organism','tid','confidence_score','chembl_id','src_id','assay_cell_type','bao_format','assay_category'],
 'docs': ['doc_id','journal','year','volume','first_page','pubmed_id','doi','chembl_id','title','doc_type','src_id','authors','patent_id'],
 'target_dictionary': ['tid','target_type','pref_name','organism','chembl_id'],
 'target_components': ['tid','component_id'],
 'component_class': ['component_id','protein_class_id'],
 'protein_classification': ['protein_class_id','parent_id','pref_name','short_name','protein_class_desc','class_level'],
 'molecule_hierarchy': ['molregno','parent_molregno','active_molregno'],
 'compound_structures': ['molregno','canonical_smiles','standard_inchi_key'],
 'molecule_dictionary': ['molregno','pref_name','chembl_id','structure_type','inorganic_flag','polymer_flag','molecule_type'],
 'source': ['src_id','src_description','src_short_name'],
 'version': None,
}
ACT_TYPES = {'MIC','IC50','Ki','GI50','CC50'}
act_quick=[]
act_q_s=[f"'{t}'" for t in ACT_TYPES]
TUP=re.compile(r"\((?:[^'()]+|'(?:[^'\\]+|\\.)*')*\)")
db = sqlite3.connect(sys.argv[1]); db.execute('PRAGMA journal_mode=OFF'); db.execute('PRAGMA synchronous=OFF')
cols = {}; cur_table=None; in_create=False; colbuf=[]
counts = {}; bad = {}
buf = {}
def flush(t):
    if buf.get(t):
        k = KEEP[t]
        db.executemany(f'INSERT INTO {t} VALUES ({",".join("?"*len(k))})', buf[t]); buf[t]=[]
t0=time.time(); n=0
for raw in sys.stdin.buffer:
    n+=1
    if n % 5000000 == 0:
        print(n, round(time.time()-t0), {k:v for k,v in counts.items()}, flush=True)
    if raw.startswith(b'CREATE TABLE `'):
        cur_table = raw.split(b'`')[1].decode(); in_create=True; colbuf=[]; continue
    if in_create:
        s=raw.strip()
        if s.startswith(b'`'):
            colbuf.append(s.split(b'`')[1].decode()); continue
        if s.startswith(b')'):
            in_create=False; cols[cur_table]=colbuf
            if cur_table in KEEP:
                k = KEEP[cur_table] or colbuf
                KEEP[cur_table]=k
                db.execute(f'CREATE TABLE {cur_table} ({",".join(k)})')
                counts[cur_table]=0; bad[cur_table]=0; buf[cur_table]=[]
        continue
    if not raw.startswith(b'INSERT INTO `'): continue
    t = raw[13:raw.index(b'`',13)].decode()
    if t not in KEEP or t not in cols: continue
    line = raw.decode('utf-8', errors='replace')
    i = line.index(' VALUES (')+8
    idx = [cols[t].index(c) for c in KEEP[t]]
    for m in TUP.finditer(line, i):
        tup = m.group()[1:-1]
        if t=='activities' and not any(q in tup for q in act_q_s): continue
        try:
            row = next(csv.reader([tup], quotechar="'", escapechar='\\', doublequote=False))
        except Exception:
            bad[t]+=1; continue
        if len(row)!=len(cols[t]):
            bad[t]+=1; continue
        if t=='activities' and row[cols[t].index('standard_type')] not in ACT_TYPES: continue
        buf[t].append([None if row[j]=='NULL' else row[j] for j in idx]); counts[t]+=1
        if len(buf[t])>=50000: flush(t)
for t in list(buf): flush(t)
db.commit()
print('DONE', n, round(time.time()-t0), counts, 'bad', bad, flush=True)
