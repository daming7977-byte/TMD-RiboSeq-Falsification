import fs from 'node:fs/promises';
import {Workbook} from '@oai/artifact-tool';
const tables=JSON.parse(await fs.readFile(process.argv[2],'utf8'));
function col(n){let s='';while(n){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26)}return s}
function field(x){const s=String(x??'UNKNOWN / NOT DOCUMENTED');return /[\t\r\n"]/.test(s)?'"'+s.replaceAll('"','""')+'"':s}
for(const [name,rows] of Object.entries(tables)){
 const w=Workbook.create(),s=w.worksheets.add('Audit');
 if(!rows.every(r=>r.length===rows[0].length))throw Error('Ragged '+name);
 // Escape ISO date text so the workbook does not reduce timestamp precision.
 const isDateText=x=>typeof x==='string'&&/^\d{4}-\d{2}-\d{2}(?:T.*)?$/.test(x);
 const encoded=rows.map(r=>r.map(x=>isDateText(x)?"'"+x:x));
 const range=s.getRange(`A1:${col(rows[0].length)}${rows.length}`);range.values=encoded;
 const actual=range.values.map((r,ri)=>r.map((x,ci)=>isDateText(rows[ri][ci])?String(x).slice(1):x));
 for(let r=0;r<rows.length;r++)for(let c=0;c<rows[0].length;c++)if(JSON.stringify(actual[r]?.[c])!==JSON.stringify(rows[r][c]))throw Error('Value roundtrip mismatch '+name+' row '+r+' col '+c+' expected '+JSON.stringify(rows[r][c])+' actual '+JSON.stringify(actual[r]?.[c]));
 const check=await w.inspect({kind:'table',range:'Audit!A1:D3',include:'values',tableMaxRows:3,tableMaxCols:4});
 await fs.writeFile(name,actual.map(r=>r.map(field).join('\t')).join('\n')+'\n');
 console.log(name,rows.length-1,'rows; roundtrip verified');
}
