"""Presentation-only terminology patch to frozen SVG; no data/statistics computation."""
import argparse,copy,re
from pathlib import Path
from lxml import etree as E

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--input',type=Path,required=True)
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--figure',choices=['Figure_2','Figure_3','Supplementary_Figure_2'],required=True)
 a=p.parse_args()
 if a.output.exists():p.error('Output exists; refusing overwrite')
 before=E.parse(str(a.input));root=copy.deepcopy(before.getroot());edits=0
 for el in root.xpath('//*[local-name()="text"]'):
  for node in el.iter():
   if node.text:
    old=node.text;new=re.sub(r'\bPreregistered\b','Prespecified',old)
    if a.figure=='Figure_2':new=new.replace('were frozen before biological reveal.','were frozen before TMD-signal reveal.')
    if new!=old:node.text=new;edits+=1
 assert edits=={'Figure_2':3,'Figure_3':1,'Supplementary_Figure_2':1}[a.figure], 'Unexpected source text'
 x=copy.deepcopy(before.getroot());y=copy.deepcopy(root)
 for tree in (x,y):
  for node in tree.iter():node.text=None;node.tail=None
 assert E.tostring(x,method='c14n')==E.tostring(y,method='c14n'), 'Non-text change'
 def numbers(n):return re.findall(r'[-+−]?\d+(?:\.\d+)?',' '.join(n.itertext()))
 assert numbers(before.getroot())==numbers(root), 'Numerical text change'
 E.ElementTree(root).write(str(a.output),encoding='utf-8',xml_declaration=True)
 print('SVG terminology updated; numerical text and all non-text elements unchanged.')
if __name__=='__main__':main()
