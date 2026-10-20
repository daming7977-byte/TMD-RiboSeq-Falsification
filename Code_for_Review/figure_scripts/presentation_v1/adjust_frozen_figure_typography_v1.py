#!/usr/bin/env python3
"""Adjust only selected SVG text layout; never load data or execute analysis."""
import argparse, copy, hashlib, io, json, re, subprocess
from pathlib import Path
from lxml import etree as ET
from pypdf import PdfReader, PdfWriter
from pypdf.generic import RectangleObject

NS='http://www.w3.org/2000/svg'
def gid(root, ident):
    found=root.xpath('//*[@id=$x]',x=ident)
    assert len(found)==1,ident
    return found[0]
def words(group):return ' '.join(' '.join(group.itertext()).split())
def text(root, ident):
    elements=gid(root,ident).findall('{%s}text'%NS)
    assert len(elements)==1
    return elements[0]
def move(element,x=None,y=None):
    element.attrib.pop('transform',None)
    if x is not None:element.set('x',str(x))
    if y is not None:element.set('y',str(y))
def reflow(root,ident,lines,x,ys,style):
    group=gid(root,ident)
    old=words(group)
    assert old==' '.join(lines),(ident,old,lines)
    for child in list(group):group.remove(child)
    for line,y in zip(lines,ys):
        node=ET.SubElement(group,'{%s}text'%NS,style=style,x=str(x),y=str(y))
        node.text=line
    assert words(group)==old

def explicit_glyph_positions(root):
    """Expand SVG per-character x lists without moving any glyph (librsvg compatibility)."""
    expanded=0
    for node in list(root.xpath('//*[local-name()="tspan"]')):
        positions=node.get('x','').split()
        if len(positions)<=1:continue
        assert len(positions)==len(node.text)
        parent=node.getparent();index=parent.index(node);tail=node.tail
        expected=[(c,x,node.get('y'),node.get('style')) for c,x in zip(node.text,positions)]
        made=[]
        for char,x in zip(node.text,positions):
            new=copy.deepcopy(node);new.set('x',x);new.text=char;new.tail=None
            parent.insert(index,new);index+=1;made.append(new)
        made[-1].tail=tail;parent.remove(node)
        assert [(n.text,n.get('x'),n.get('y'),n.get('style')) for n in made]==expected
        expanded+=1
    return expanded

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--figure1-svg',required=True,type=Path)
    ap.add_argument('--figure4-svg',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--renderer',default='rsvg-convert')
    # Recorded frozen PDF media boxes; SVG and PDF historical tight bounding boxes differ slightly.
    ap.add_argument('--figure1-page-size',default='516.65703125,465.2315')
    ap.add_argument('--figure4-page-size',default='529.9962606383,519.5752066116')
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    audit=[];compatibility=[]
    for number,path in [(1,a.figure1_svg),(4,a.figure4_svg)]:
        original=ET.parse(str(path))
        compatibility.append(dict(figure=number,expanded_tspan_count=explicit_glyph_positions(original.getroot()),glyph_coordinates_and_characters='IDENTICAL_TO_FROZEN_SVG',reason='Explicit per-character SVG coordinates for librsvg compatibility'))
        root=copy.deepcopy(original.getroot())
        if number==1:
            changed=['text_24','text_26','text_28','text_48','text_69','text_70','text_71','text_72','text_73','text_74']
            # Panel B: retain all words/values and move the negative-bar label above its bar.
            move(text(root,'text_24'),y=163.5)
            t=text(root,'text_28');t.set('style',t.get('style').replace('8.5px','7.4px'))
            reflow(root,'text_26',['pseudo','> real'],307.8,[56,64],"fill: #d55e00; font: 700 6.7px 'Arial'; text-anchor: middle")
            # Panel D: separate the unchanged effect caption from the FAIL caption.
            move(text(root,'text_48'),y=305.5)
            # Panel F: reflow text inside its existing, unchanged boxes.
            move(text(root,'text_69'),y=295)
            move(text(root,'text_71'),y=295)
            reflow(root,'text_70',['strong broad','signal ≠','TMD-specific'],404.641688,[303,310,317],"font: 6.2px 'Arial'; text-anchor: middle")
            reflow(root,'text_72',['coordinate-null','PASS ≠ TMD','validation'],468.428625,[303,310,317],"font: 6.2px 'Arial'; text-anchor: middle")
            reflow(root,'text_73',['Credible TMD-specific','inference requires'],436.535156,[350,359],"font: 700 7.5px 'Arial'; text-anchor: middle")
            reflow(root,'text_74',['matched controls +','technical robustness','+ positional specificity','+ preserved inference hierarchy'],436.535156,[367,374,381,388],"font: 6.2px 'Arial'; text-anchor: middle")
        else:
            changed=['text_123']
            reflow(root,'text_123',['Direct topology × compartment','interaction'],443.558448,[263.5,273],"font: 8.5px 'Arial'; text-anchor: middle")
        # Remove only whitelisted text groups in paired copies; everything else must match.
        left=copy.deepcopy(original.getroot());right=copy.deepcopy(root)
        for ident in changed:
            before=gid(left,ident);after=gid(right,ident)
            assert words(before)==words(after),(number,ident)
            audit.append(dict(figure=number,element=ident,text=words(before),before=ET.tostring(before,encoding='unicode'),after=ET.tostring(after,encoding='unicode')))
            before.getparent().remove(before);after.getparent().remove(after)
        assert ET.tostring(left,method='c14n')==ET.tostring(right,method='c14n'),number
        out=a.output/f'Figure_{number}.svg'
        ET.ElementTree(root).write(str(out),encoding='utf-8',xml_declaration=True)
        subprocess.run([a.renderer,'--format','pdf','--background-color','white','--output',str(a.output/f'Figure_{number}.pdf'),str(out)],check=True)
        pdf=a.output/f'Figure_{number}.pdf'
        reader=PdfReader(io.BytesIO(pdf.read_bytes()));page=reader.pages[0]
        content=page.get_contents().get_data()
        width,height=map(float,getattr(a,f'figure{number}_page_size').split(','))
        page.mediabox=RectangleObject([0,0,width,height])
        writer=PdfWriter();writer.add_page(page)
        with pdf.open('wb') as handle:writer.write(handle)
        assert PdfReader(pdf).pages[0].get_contents().get_data()==content
        subprocess.run([a.renderer,'--format','png','--background-color','white','--dpi-x','400','--dpi-y','400','--output',str(a.output/f'Figure_{number}.png'),str(out)],check=True)
    (a.output/'FIGURE_TEXT_LAYOUT_DIFF_v1.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False))
    (a.output/'SVG_GLYPH_COMPATIBILITY_v1.json').write_text(json.dumps(compatibility,indent=2))
    print(json.dumps(dict(text_groups_changed=len(audit),nontext_svg_elements='UNCHANGED',text_content='UNCHANGED',scientific_execution='NOT_RUN')))

if __name__=='__main__':main()
