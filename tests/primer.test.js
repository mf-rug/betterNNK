import {test} from 'node:test';
import assert from 'node:assert/strict';
import {CODE,IUPAC,METHODS,expand,reverseComplement,normalize,frames,guessFrame,generate,tsv,resolveMethod,distribution,libraryStats,mixedLibraryStats} from '../web/primer.js';
test('standard genetic code includes 64 codons and exactly three stops',()=>{
 assert.equal(Object.keys(CODE).length,64);assert.deepEqual(Object.keys(CODE).filter(c=>CODE[c]==='*').sort(),['TAA','TAG','TGA']);
 assert.equal(CODE.ATG,'M');assert.equal(CODE.TGG,'W');assert.equal(CODE.AGT,'S');
});
for(const method of ['22c','20c'])test(`${method}: exact coverage, no stops, correct molar weights`,()=>{
 const counts={},codons=[];for(const [c,weight] of METHODS[method].parts){const expanded=expand(c);assert.equal(expanded.length,weight);for(const triplet of expanded){codons.push(triplet);counts[CODE[triplet]]=(counts[CODE[triplet]]||0)+1;}}
 assert.equal(new Set(codons).size,METHODS[method].total);assert.equal(Object.keys(counts).sort().join(''),'ACDEFGHIKLMNPQRSTVWY');
 for(const [aa,count] of Object.entries(counts))assert.equal(count,method==='22c'&&'LV'.includes(aa)?2:1);
});
test('reverse codes match Tang main text and Kille SI Table S4',()=>{
 assert.deepEqual(['NDT','VHG','VMA','ATG','TGG'].map(reverseComplement),['AHN','CDB','TKB','CAT','CCA']);
 for(const code of Object.keys(IUPAC))assert.equal(reverseComplement(reverseComplement(code)),code);
});
test('published CHMO-426 primer pair reproduced from Kille SI',()=>{
 const seq='CATGTTCATGGTGGCTGGACCGAATGGC';
 const rows=generate(seq,'+2',13,'22c','CHMO426',true);
 assert.equal(rows[0].sequence,'CATGTTCATGGTGNDTGGACCGAATGGC');
 assert.equal(rows[1].sequence,'GCCATTCGGTCCAHNCACCATGAACATG');
 assert.equal(rows[3].sequence,'GCCATTCGGTCCCDBCACCATGAACATG');
});
test('Tang published forward-primer flanks preserved',()=>{
 const rows=generate('TCTTCAACGCTCCGACATA','+3',8,'20c');
 assert.deepEqual(rows.map(r=>r.sequence),['TCTTCAACNDTCCGACATA','TCTTCAACVMACCGACATA','TCTTCAACATGCCGACATA','TCTTCAACTGGCCGACATA']);
});
test('all six frames: coordinates, coding replacements, both orientations and partial ends',()=>{
 for(const seq of ['ATGACGTTAGCTGACGTTGCA','AATGACGTTAGCTGACGTTGCA','AAATGACGTTAGCTGACGTTGCA']){
 for(const f of frames(seq))for(const c of f.codons)for(const method of ['22c','20c']){
 const rows=generate(seq,f.id,c.start,method,'x',true);
 for(let i=0;i<rows.length;i+=2){const r=rows[i];assert.equal(r.sequence.length,seq.length);assert.equal(r.sequence.slice(0,c.start),seq.slice(0,c.start));assert.equal(r.sequence.slice(c.end),seq.slice(c.end));assert.equal(rows[i+1].sequence,reverseComplement(r.sequence));const coding=f.strand===1?r.replacement:reverseComplement(r.replacement);assert.equal(coding,r.codon);}
 assert.equal(c.triplet,f.strand===1?seq.slice(c.start,c.end):reverseComplement(seq.slice(c.start,c.end)));
 }
 }
});
test('input validation never silently drops non-DNA symbols',()=>{
 assert.equal(normalize(' acg\nt T ').sequence,'ACGTT');
 for(const input of ['ATN','ATU','123ATG','>name\nATG','AT','A'.repeat(3001)])assert.ok(normalize(input).error);
 assert.equal(normalize('').error,'');
});
test('frame suggestion is deterministic and excludes empty frames',()=>{
 assert.equal(guessFrame(frames('ATG')),'+1');
 assert.equal(guessFrame(frames('TAACCC')),'+2');
 assert.equal(frames('ATG').filter(f=>f.codons.length).length,2);
});
test('invalid selections rejected; names cannot inject TSV fields or formulas',()=>{
 assert.throws(()=>generate('ATGCCC','+1',1,'22c'));
 assert.throws(()=>generate('ATGCCC','+1',0,'bogus'));
 const rows=generate('ATGCCC','+1',0,'22c','=evil\t\n<script>');
 assert.ok(rows.every(r=>/^[A-Za-z0-9_.-]+$/.test(r.name)));
 const plain=tsv(rows,false,false).trim().split('\n');assert.equal(plain.length,3);assert.ok(plain.every(r=>r.split('\t').length===2));
 assert.equal(tsv(rows,true).trim().split('\n')[0].split('\t').length,10);
});

test('NNK and NNS distribution includes stops without renormalizing sense codons',()=>{
 for(const recipe of [resolveMethod('NNK'),resolveMethod('other','NNS')]){
 const d=distribution(recipe);assert.equal(d.codons,32);assert.equal(d.aminoAcids,20);assert.equal(d.stopCodons,1);
 assert.equal(d.probabilities['*'],1/32);assert.equal(d.probabilities.L,3/32);assert.equal(d.probabilities.W,1/32);
 }
});
test('restricted and custom codons give exact normalized distributions',()=>{
 for(const [codon,count,aas,stops] of [['NDT',12,12,0],['VMA',6,6,0],['VRK',12,9,0],['NNN',64,20,3],['ATG',1,1,0],['TAA',1,0,1]]){
 const d=distribution(resolveMethod('other',codon));assert.equal(d.codons,count);assert.equal(d.aminoAcids,aas);assert.equal(d.stopCodons,stops);
 assert.ok(Math.abs(Object.values(d.probabilities).reduce((a,b)=>a+b,0)-1)<1e-12);
 }
 assert.ok(Math.abs(distribution(resolveMethod('20c')).probabilities.A-0.05)<1e-12);
 assert.ok(Math.abs(distribution(resolveMethod('22c')).probabilities.L-2/22)<1e-12);
 for(const c of ['', 'NN', 'NNNN', 'N-U', 'XYZ'])assert.equal(resolveMethod('other',c),null);
 assert.equal(resolveMethod('other',' ndt ').parts[0][0],'NDT');
});
test('shared names have no input suffix; custom reverse triplets and TSV agree',()=>{
 const rows=generate('ATGCCC','-1',3,'other','MyPrimer',true,'NNS');
 assert.equal(rows[0].name,'MyPrimer_nt4_other_NNS');assert.equal(rows[1].name,'MyPrimer_nt4_other_NNS_rc');
 assert.equal(rows[0].sequence,'ATGSNN');assert.equal(rows[1].sequence,'NNSCAT');
 assert.ok(tsv(rows).includes('MyPrimer_nt4_other_NNS\tATGSNN'));
 assert.throws(()=>generate('ATGCCC','+1',0,'other','p',false,'XYZ'));
});

test('two-site 22c combinatorial weights match Kille SI Table S1',()=>{
 const rows=generate('ATGGCTGAA','+1',[0,6],'22c','test',false,'NNS','combinatorial');
 assert.equal(rows.length,9);
 assert.deepEqual(rows.map(r=>r.ratio).sort((a,b)=>a-b),[1,9,9,12,12,81,108,108,144]);
 assert.equal(rows.reduce((sum,r)=>sum+r.ratio,0),484);
 assert.ok(Math.abs(rows.reduce((sum,r)=>sum+r.fraction,0)-1)<1e-12);
 assert.equal(new Set(rows.map(r=>r.pool)).size,1);
 assert.equal(new Set(rows.map(r=>r.name)).size,9);
 assert.ok(rows.every(r=>r.sequence.slice(3,6)==='GCT'));
 assert.equal(rows.find(r=>r.codon==='NDT/VHG').sequence,'NDTGCTVHG');
});
test('separate libraries preserve all other selected sites and normalize each strand pool',()=>{
 const rows=generate('ATGGCTGAA','+1',[0,6],'22c','test',true,'NNS','separate');
 assert.equal(rows.length,12);assert.equal(new Set(rows.map(r=>r.pool)).size,4);
 for(const pool of new Set(rows.map(r=>r.pool)))assert.ok(Math.abs(rows.filter(r=>r.pool===pool).reduce((sum,r)=>sum+r.fraction,0)-1)<1e-12);
 for(const row of rows.filter(r=>r.strand==='input')){
 assert.equal(row.starts.length,1);
 assert.equal(row.sequence.slice(row.starts[0]===0?6:0,row.starts[0]===0?9:3),row.starts[0]===0?'GAA':'ATG');
 }
 assert.ok(tsv(rows,true).includes('site_nt7_complement'));
});
test('reverse combinatorial selection assigns components to ascending input positions',()=>{
 const rows=generate('ATGGCTGAA','-1',[6,0],'20c','p',true,'NNS','combinatorial');
 assert.equal(rows.length,32);
 const row=rows.find(r=>r.strand==='input'&&r.codon==='NDT/VMA');
 assert.equal(row.sequence,'AHNGCTTKB');assert.equal(row.ratio,72);
 assert.equal(rows.find(r=>r.name===row.name+'_rc').sequence,'VMAAGCNDT');
});
test('multi-site preview counts, sense diversity and stop probability',()=>{
 const a=libraryStats(resolveMethod('20c'),3,'combinatorial',true);
 assert.equal(a.oligos,128);assert.equal(a.variants,8000n);assert.equal(a.libraries,1);
 const b=libraryStats(resolveMethod('NNK'),2,'combinatorial');
 assert.equal(b.oligos,1);assert.ok(Math.abs(b.stopProbability-(1-(31/32)**2))<1e-12);
 const c=libraryStats(resolveMethod('NNK'),2,'separate');
 assert.equal(c.libraries,2);assert.equal(c.variants,20n);assert.equal(c.stopProbability,1/32);
});
test('selection sorting, deduplication, invalid selections and explosion guard',()=>{
 assert.equal(generate('ATGGCTGAA','+1',[6,0,6],'22c').length,6);
 assert.throws(()=>generate('ATGGCTGAA','+1',[],'22c'));
 assert.throws(()=>generate('ATGGCTGAA','+1',[0,1],'22c'));
 assert.throws(()=>generate('ATG'.repeat(6),'+1',[0,3,6,9,12,15],'20c','p',false,'NNS','combinatorial'),/exceeds/);
 assert.equal(generate('ATG'.repeat(6),'+1',[0,3,6,9,12,15],'20c').length,24);
});

test('mixed 22c/Tang components, weights and fractions normalize in each pool',()=>{
 const overrides={6:{method:'20c'}};
 const rows=generate('ATGGCTGAA','+1',[6,0],'22c','p',true,'NNS','combinatorial',overrides);
 assert.equal(rows.length,24);
 const input=rows.filter(r=>r.strand==='input');
 assert.equal(input.reduce((n,r)=>n+r.ratio,0),440);
 assert.ok(Math.abs(input.reduce((n,r)=>n+r.fraction,0)-1)<1e-12);
 const row=input.find(r=>r.codon==='VHG/VMA');assert.equal(row.sequence,'VHGGCTVMA');assert.equal(row.ratio,54);assert.equal(row.method,'22c/20c');
 assert.ok(tsv(rows,true).includes('22c/20c'));
});
test('mixed restricted libraries report per-site diversity and product stop probabilities',()=>{
 const methods=[resolveMethod('NNK'),resolveMethod('other','NDT'),resolveMethod('other','NNN')];
 const stats=mixedLibraryStats(methods,'combinatorial');
 assert.equal(stats.variants,4800n);assert.equal(stats.oligos,1);
 assert.ok(Math.abs(stats.stopProbability-(1-(31/32)*(61/64)))<1e-12);
 const separate=mixedLibraryStats(methods,'separate');
 assert.deepEqual(separate.perSite.map(s=>s.variants),[20n,12n,20n]);
 assert.deepEqual(separate.perSite.map(s=>s.stopProbability),[1/32,0,3/64]);
});
test('mixed reverse-frame custom triplets and separate library counts',()=>{
 const overrides={0:{method:'other',custom:'VMA'},6:{method:'NNK'}};
 const rows=generate('ATGGCTGAA','-1',[0,6],'22c','p',true,'NNS','combinatorial',overrides);
 assert.equal(rows.length,2);assert.equal(rows[0].sequence,'TKBGCTMNN');assert.equal(rows[1].sequence,'NNKAGCVMA');
 assert.equal(rows[0].method,'VMA/NNK');
 const separate=generate('ATGGCTGAA','+1',[0,6],'22c','p',false,'NNS','separate',{6:{method:'other',custom:'NDT'}});
 assert.equal(separate.length,4);assert.equal(separate.find(r=>r.start==='7').sequence,'ATGGCTNDT');
 assert.throws(()=>generate('ATGGCTGAA','+1',[0,6],'22c','p',false,'NNS','combinatorial',{6:{method:'other',custom:'XY'}}));
});
