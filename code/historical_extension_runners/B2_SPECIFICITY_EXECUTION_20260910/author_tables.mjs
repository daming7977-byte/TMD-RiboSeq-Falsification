import fs from 'node:fs/promises';
import crypto from 'node:crypto';
import {Workbook} from '@oai/artifact-tool';
const input=process.argv[2]??'tables.json';
const tables=JSON.parse(await fs.readFile(input,'utf8'));
function col(n){let s='';while(n){n--;s=String.fromCharCode(65+n%26)+s;n=Math.floor(n/26)}return s}
function field(x){const s=String(x??'UNKNOWN / NOT DOCUMENTED');return /[\t\r\n"]/.test(s)?'"'+s.replaceAll('"','""')+'"':s}
const hashes=[];
for(const [name,rows] of Object.entries(tables)){
 const w=Workbook.create(),s=w.worksheets.add('Audit');
 if(!rows.every(r=>r.length===rows[0].length))throw Error('Ragged '+name);
 const range=s.getRange(`A1:${col(rows[0].length)}${rows.length}`);range.values=rows;
 const text=range.values.map(r=>r.map(field).join('\t')).join('\n')+'\n';
 await fs.writeFile(name,text);await fs.chmod(name,0o444);
 hashes.push(crypto.createHash('sha256').update(text).digest('hex')+'  '+name);
 console.log(name,rows.length-1);
}
if(input==='tables.json'){
 await fs.writeFile('NUMERICAL_OUTPUTS_SHA256.txt',hashes.join('\n')+'\n');
 await fs.chmod('NUMERICAL_OUTPUTS_SHA256.txt',0o444);
 await fs.writeFile('NUMERICAL_OUTPUT_FREEZE_TIME.txt',new Date().toISOString()+'; prior to aggregate visualization\n');
 await fs.chmod('NUMERICAL_OUTPUT_FREEZE_TIME.txt',0o444);
}
