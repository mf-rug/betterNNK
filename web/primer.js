// Standard genetic code. Codons ordered by T,C,A,G at each position.
const bases = 'TCAG';
const aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
export const CODE = Object.fromEntries([...bases].flatMap((a,i)=>[...bases].flatMap((b,j)=>[...bases].map((c,k)=>[a+b+c,aas[i*16+j*4+k]]))));
export const IUPAC = {A:'A',C:'C',G:'G',T:'T',R:'AG',Y:'CT',S:'CG',W:'AT',K:'GT',M:'AC',B:'CGT',D:'AGT',H:'ACT',V:'ACG',N:'ACGT'};
const complements = {A:'T',T:'A',C:'G',G:'C',R:'Y',Y:'R',S:'S',W:'W',K:'M',M:'K',B:'V',V:'B',D:'H',H:'D',N:'N'};
export const METHODS = {
 '22c':{name:'22c-trick',total:22,parts:[['NDT',12],['VHG',9],['TGG',1]]},
 'NNK':{name:'NNK',total:1,parts:[['NNK',1]]},
 '20c':{name:'20c · Tang',total:20,parts:[['NDT',12],['VMA',6],['ATG',1],['TGG',1]]}
};
export function resolveMethod(id, custom='NNS'){
 if(id!=='other')return METHODS[id]??null;
 const codon=custom.trim().toUpperCase();
 if(!/^[ACGTRYSWKMBDHVN]{3}$/.test(codon))return null;
 return {name:codon,total:1,parts:[[codon,1]]};
}
export function distribution(method){
 const probabilities=Object.fromEntries([...'ACDEFGHIKLMNPQRSTVWY*'].map(a=>[a,0]));
 const unique=new Set();
 for(const [codon,weight] of method.parts){
 const variants=expand(codon);
 for(const variant of variants){unique.add(variant);probabilities[CODE[variant]]+=weight/method.total/variants.length;}
 }
 return {probabilities,codons:unique.size,aminoAcids:Object.entries(probabilities).filter(([a,p])=>a!=='*'&&p>0).length,stopCodons:[...unique].filter(c=>CODE[c]==='*').length};
}
export function normalize(raw){
 const sequence = raw.replace(/\s/g,'').toUpperCase();
 if (!sequence) return {sequence:'',error:''};
 if (/[^ACGT]/.test(sequence)) return {sequence:'',error:'Use only A, C, G and T. Spaces and line breaks are accepted; remove FASTA headers, numbers and existing degeneracy codes.'};
 if(sequence.length<3) return {sequence:'',error:'Enter at least 3 bases to select a complete codon.'};
 if(sequence.length>3000) return {sequence:'',error:'Please enter a primer of at most 3,000 bases.'};
 return {sequence,error:''};
}
export function reverseComplement(seq){return [...seq].reverse().map(b=>{if(!complements[b])throw Error('Invalid IUPAC base');return complements[b];}).join('');}
export function expand(codon){return [...codon].reduce((all,b)=>all.flatMap(s=>[...IUPAC[b]].map(x=>s+x)),['']);}
export function frames(seq){
 return [1,-1].flatMap(strand=>[0,1,2].map(offset=>{
 const coding = strand===1?seq:reverseComplement(seq);
 const codons=[];
 for(let pos=offset;pos+3<=seq.length;pos+=3){
 const triplet=coding.slice(pos,pos+3);
 const start=strand===1?pos:seq.length-pos-3;
 codons.push({triplet,aa:CODE[triplet],start,end:start+3,index:codons.length+1});
 }
 return {id:`${strand===1?'+':'-'}${offset+1}`,strand,offset,codons,stops:codons.filter(c=>c.aa==='*').length};
 }));
}
export function guessFrame(all){return all.find(f=>f.codons.length&&f.stops===0)?.id??all.find(f=>f.codons.length)?.id??'+1';}
export const MAX_OLIGOS=2048;
export function mixedLibraryStats(methods,mode='separate',paired=false){
 if(!['separate','combinatorial'].includes(mode))throw Error('Unknown library mode.');
 const distributions=methods.map(distribution);
 const perSite=distributions.map(d=>({variants:BigInt(d.aminoAcids),stopProbability:d.probabilities['*']}));
 return {oligos:methods.length?(mode==='combinatorial'?methods.reduce((n,m)=>n*m.parts.length,1):methods.reduce((n,m)=>n+m.parts.length,0))*(paired?2:1):0,
 libraries:methods.length?(mode==='combinatorial'?1:methods.length):0,
 variants:methods.length?(mode==='combinatorial'?perSite.reduce((n,d)=>n*d.variants,1n):perSite[0].variants):0n,
 stopProbability:methods.length?(mode==='combinatorial'?1-perSite.reduce((n,d)=>n*(1-d.stopProbability),1):perSite[0].stopProbability):0,perSite};
}
export function libraryStats(method,n,mode='separate',paired=false){return mixedLibraryStats(Array(n).fill(method),mode,paired);}
export function generate(seq,frameId,start,methodId,prefix='primer',paired=false,custom='NNS',mode='separate',overrides={}){
 const starts=[...new Set(Array.isArray(start)?start:[start])].sort((a,b)=>a-b);
 const frame=frames(seq).find(f=>f.id===frameId);
 const configs=Object.fromEntries(starts.map(s=>[s,overrides[s]??{method:methodId,custom}]));
 const recipes=Object.fromEntries(starts.map(s=>[s,resolveMethod(configs[s].method,configs[s].custom)]));
 if(!starts.length||!frame||Object.values(recipes).some(m=>!m)||starts.some(s=>!frame.codons.some(c=>c.start===s)))throw Error('Select complete codons and a supported method.');
 if(mixedLibraryStats(starts.map(s=>recipes[s]),mode,paired).oligos>MAX_OLIGOS)throw Error(`This design exceeds ${MAX_OLIGOS.toLocaleString()} oligos. Select fewer sites or use separate libraries.`);
 const safePrefix=prefix.trim().replace(/[^a-zA-Z0-9_.-]+/g,'_').slice(0,60)||'primer';
 const groups=mode==='combinatorial'?[starts]:starts.map(s=>[s]);
 return groups.flatMap(targets=>{
 const combinations=targets.reduce((all,pos)=>all.flatMap(parts=>recipes[pos].parts.map(p=>[...parts,p])),[[]]);
 return combinations.flatMap(parts=>{
 const codons=parts.map(p=>p[0]);
 const replacements=codons.map(c=>frame.strand===1?c:reverseComplement(c));
 let sequence=seq;targets.forEach((pos,i)=>{sequence=sequence.slice(0,pos)+replacements[i]+sequence.slice(pos+3);});
 const ratio=parts.reduce((v,p)=>v*p[1],1),fraction=parts.reduce((v,p,i)=>v*p[1]/recipes[targets[i]].total,1);
 const positions=targets.map(s=>s+1).join('-');
 const methodNames=targets.map(s=>configs[s].method==='other'?recipes[s].name:configs[s].method);
 const methodLabel=new Set(targets.map(s=>configs[s].method)).size===1?configs[targets[0]].method:'mixed';
 const name=`${safePrefix}_nt${positions}_${methodLabel}_${codons.join('-')}`;
 const library=`${mode==='combinatorial'?'combined':'site'}_nt${positions}`;
 const row={name,sequence,codon:codons.join('/'),replacement:replacements.join('/'),ratio,fraction,pool:`${library}_input`,strand:'input',library,mode,starts:targets,start:targets.map(s=>s+1).join(','),frame:frameId,method:methodNames.join('/')};
 return paired?[row,{...row,name:`${name}_rc`,sequence:reverseComplement(sequence),replacement:replacements.map(reverseComplement).join('/'),pool:`${library}_complement`,strand:'complement'}]:[row];
 });
 });
}
export function tsv(rows,detailed=false,header=true){
 const cols=detailed?['Name','Sequence (5–3)','Pool','Coding codons (input-position order)','Ordered triplets (input-position order)','Molar parts','Pool fraction','Frame','Input starts (1-based)','Method']:['Name','Sequence'];
 const data=rows.map(r=>detailed?[r.name,r.sequence,r.pool,r.codon,r.replacement,r.ratio,r.fraction.toPrecision(12),r.frame,r.start,r.method]:[r.name,r.sequence]);
 return (header?[cols,...data]:data).map(r=>r.join('\t')).join('\n')+'\n';
}
