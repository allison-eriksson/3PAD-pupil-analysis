from pathlib import Path
from zipfile import ZipFile
import re,json,hashlib
from openpyxl import load_workbook
import argparse
parser = argparse.ArgumentParser(description="Patch verified values into a new workbook without changing its format.")
parser.add_argument("--original-workbook", type=Path, required=True)
parser.add_argument("--audit-dir", type=Path, required=True)
parser.add_argument("--output-workbook", type=Path, required=True)
args = parser.parse_args()
source, audit, out = args.original_workbook, args.audit_dir, args.output_workbook
updates = json.loads((audit/'workbook_updates.json').read_text())
with ZipFile(source) as a,ZipFile(out,'x') as b:
 for info in a.infolist():
  data=a.read(info.filename)
  if info.filename=='xl/worksheets/sheet1.xml':
   text=data.decode('utf-8')
   for u in updates:
    value=u['corrected']
    assert isinstance(value,(float,int)) and abs(value-u['corrected'])<1e-12
    pattern=r'(<c\b[^>]*\br="'+u['cell']+r'"[^>]*>)(.*?)(</c>)'
    def replace(m):
     middle,n=re.subn(r'(<v>).*?(</v>)',lambda v:v[1]+format(value,'.17g')+v[2],m[2],count=1)
     assert n==1
     return m[1]+middle+m[3]
    text,n=re.subn(pattern,replace,text,count=1);assert n==1,u['cell']
   data=text.encode('utf-8')
  b.writestr(info,data)
a=load_workbook(source);b=load_workbook(out);assert a.sheetnames==b.sheetnames
allowed={u['cell']:u for u in updates};changes=[]
for s1,s2 in zip(a,b):
 assert (s1.max_row,s1.max_column)==(s2.max_row,s2.max_column)==(100,830)
 for row1,row2 in zip(s1,s2):
  for c1,c2 in zip(row1,row2):
   assert c1._style==c2._style and c1.data_type==c2.data_type,(c1.coordinate,'style/type changed')
   assert (c1.value is None)==(c2.value is None),(c1.coordinate,'missingness changed')
   if c1.value!=c2.value:
    assert c1.coordinate in allowed,c1.coordinate
    assert abs(c2.value-allowed[c1.coordinate]['corrected'])<1e-12
    changes.append(c1.coordinate)
with ZipFile(source) as za,ZipFile(out) as zb:
 assert za.namelist()==zb.namelist()
 unchanged=[n for n in za.namelist() if n!='xl/worksheets/sheet1.xml'];assert all(za.read(n)==zb.read(n) for n in unchanged)
report=dict(rows=99,columns=830,sheets=a.sheetnames,updated_columns=sorted(set(u['variable'] for u in updates)),changed_cells=len(changes),other_cells_unchanged=True,types_and_styles_preserved=True,missingness_preserved=True,all_other_zip_parts_identical=True,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),corrected_sha256=hashlib.sha256(out.read_bytes()).hexdigest())
(audit/'workbook_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
