import {normalize,frames,guessFrame,generate,tsv,resolveMethod,distribution,mixedLibraryStats} from './primer.js';
const $=id=>document.getElementById(id);
let sequence='', all=[], frameId='+1', selected=new Set(), method='22c', rows=[], manual=false, overrides={};
const escape=s=>String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function resetStatus(){$('status').textContent='';$('copy-fallback').hidden=true;}
function inputChanged(){
 const result=normalize($('sequence').value);
 $('error').textContent=result.error;$('error').hidden=!result.error;
 $('sequence').setAttribute('aria-invalid',String(!!result.error));
 if(sequence!==result.sequence){selected=new Set();overrides={};manual=false;}
 sequence=result.sequence;all=sequence?frames(sequence):[];
 if(!manual)frameId=guessFrame(all);
 $('stats').textContent=sequence?`${sequence.length} bases · ${Math.round([...sequence].filter(b=>'GC'.includes(b)).length/sequence.length*100)}% GC`:'0 bases';
 $('empty-sequence').hidden=!!sequence;$('selection').hidden=!sequence;
 resetStatus();renderSelection();renderSiteControls();renderOutput();
}
function renderSelection(){
 if(!sequence)return;
 const frame=all.find(f=>f.id===frameId);
 $('frames').innerHTML=all.map(f=>`<button class="frame ${f.id===frameId?'selected':''}" data-frame="${f.id}" aria-pressed="${f.id===frameId}" ${f.codons.length?'':'disabled'}><strong>${f.strand===1?'→ Forward':'← Reverse'} ${f.offset+1}</strong><small class="${f.stops?'stop':''}">${f.codons.length?`${f.stops?`${f.stops} stop${f.stops===1?'':'s'}`:'Stop-free'} · ${f.codons.length} codons`:'No complete codon'}</small></button>`).join('');
 const free=all.filter(f=>f.codons.length&&f.stops===0).length;
 $('frame-note').textContent=`${manual?'Manually selected.':'Auto-selected '+frameId+'.'} ${free===1?'One frame is stop-free; confirm against your coding sequence.':free>1?`${free} frames are stop-free. This is an ambiguous guess; confirm the frame against your coding sequence.`:'Every frame contains a stop. No stop-free guess is available; confirm the frame manually.'}`;
 const ordered=[...frame.codons].sort((a,b)=>a.start-b.start);
 let html='',cursor=0;
 for(const c of ordered){
 if(c.start>cursor)html+=`<span class="partial" title="Incomplete codon; not selectable">${sequence.slice(cursor,c.start)}</span>`;
 html+=`<button class="codon ${selected.has(c.start)?'selected':''} ${c.aa==='*'?'stop':''}" data-start="${c.start}" aria-pressed="${selected.has(c.start)}" aria-label="Input bases ${c.start+1} to ${c.end}, coding codon ${c.triplet}, amino acid ${c.aa==='*'?'stop':c.aa}" title="Coding ${c.triplet} → ${c.aa}; input bases ${c.start+1}–${c.end}"><span class="position">${c.start+1}</span><span class="bases">${sequence.slice(c.start,c.end)}</span><span class="aa">${c.aa==='*'?'STOP':c.aa}</span></button>`;
 cursor=c.end;
 }
 if(cursor<sequence.length)html+=`<span class="partial" title="Incomplete codon; not selectable">${sequence.slice(cursor)}</span>`;
 $('codons').innerHTML=html;
 $('direction').textContent=frame.strand===1?'Coding direction → · amino acids below each triplet':'← Coding direction · amino acids read on the reverse complement';
 const targets=frame.codons.filter(c=>selected.has(c.start));
 $('clear-selection').disabled=!targets.length;
 $('target').textContent=targets.length?`${targets.length} selected: ${targets.map(c=>`${c.triplet} (${c.aa}) at nt ${c.start+1}–${c.end}`).join(' · ')}. ${targets.some(c=>Math.min(c.start,sequence.length-c.end)<12)?'Short flank: fewer than 12 bases on at least one side; check your PCR design.':''}`:'Choose one or more triplets above to randomize them.';
}
function customCodon(){return $('other-codon').value==='custom'?$('custom-codon').value:$('other-codon').value;}
function currentMethod(){return resolveMethod(method,customCodon());}
function renderMethod(){
 document.querySelectorAll('[data-method]').forEach(b=>{b.classList.toggle('selected',b.dataset.method===method);b.setAttribute('aria-pressed',String(b.dataset.method===method));});
 $('other-controls').hidden=method!=='other';
 $('custom-control').hidden=$('other-codon').value!=='custom';
 const recipe=currentMethod();
 $('codon-error').hidden=!!recipe;
 $('custom-codon').setAttribute('aria-invalid',String(!recipe));
 if(!recipe){
 $('method-note').textContent='Complete your custom codon to calculate the library and generate primers.';
 $('aa-count').textContent=$('codon-count').textContent=$('stop-count').textContent='—';
 $('amino-acids').replaceChildren();return;
 }
 const d=distribution(recipe);
 $('method-note').textContent=method==='22c'?'22 codons → 20 amino acids. Leu and Val each occur twice; the other amino acids occur once. No encoded stops at the target.':method==='20c'?'20 codons → 20 amino acids. Each amino acid occurs exactly once. No encoded stops at the target.':`${d.codons} codons → ${d.aminoAcids} amino acids. ${d.stopCodons?`Stop probability: ${(d.probabilities['*']*100).toFixed(2)}%.`:'No encoded stops at the target.'}`;
 renderChart();
}
function siteConfig(pos){return overrides[pos]??{method,custom:customCodon()};}
function siteRecipe(pos){const c=siteConfig(pos);return resolveMethod(c.method,c.custom);}
function renderSiteControls(){
 const starts=[...selected].sort((a,b)=>a-b);
 const choices=[['default','Shared default'],['22c','22c-trick'],['20c','20c · Tang'],['NNK','NNK'],...Array.from($('other-codon').options).map(o=>[o.value,o.textContent])];
 $('site-methods').innerHTML=starts.length?`<div class="site-heading"><strong>Per-position methods</strong><button id="reset-methods" class="text-button">Reset to shared default</button></div>`+starts.map(pos=>{
 const c=overrides[pos],value=c?(c.method==='other'?(c.preset??'custom'):c.method):'default';
 return `<div class="site-method-row"><label for="site-${pos}">nt ${pos+1}–${pos+3}</label><select id="site-${pos}" data-site="${pos}">${choices.map(([v,label])=>`<option value="${v}" ${v===value?'selected':''}>${escape(label)}</option>`).join('')}</select>${value==='custom'?`<input data-custom-site="${pos}" aria-label="Custom codon at nt ${pos+1}" placeholder="e.g. NDT" maxlength="3" value="${escape(c.custom)}" spellcheck="false">`:''}<span id="site-error-${pos}" class="site-error" role="status"></span></div>`;
 }).join(''):'';
 const old=$('chart-site').value;
 $('chart-site').innerHTML='<option value="default">Shared default</option>'+starts.map(pos=>`<option value="${pos}">nt ${pos+1}–${pos+3}</option>`).join('');
 $('chart-site').value=old==='default'||starts.includes(Number(old))?old:(starts.length?String(starts[0]):'default');
}
function renderChart(){
 const value=$('chart-site').value,recipe=value==='default'?currentMethod():siteRecipe(Number(value));
 $('distribution-scope').textContent=value==='default'?'Shared default distribution':`Distribution at input nt ${Number(value)+1}–${Number(value)+3}`;
 if(!recipe){$('aa-count').textContent=$('codon-count').textContent=$('stop-count').textContent='—';$('amino-acids').replaceChildren();return;}
 const d=distribution(recipe);
 $('aa-count').textContent=d.aminoAcids;$('codon-count').textContent=d.codons;$('stop-count').textContent=d.stopCodons;
 const max=Math.max(...Object.values(d.probabilities));
 $('amino-acids').innerHTML=Object.entries(d.probabilities).map(([aa,p])=>`<div class="bar-row ${aa==='*'?'stop-bar':''}" aria-label="${aa==='*'?'Stop':aa}: ${(p*100).toFixed(2)}%"><span>${aa==='*'?'Stop':aa}</span><div class="bar-track"><div class="bar-fill" style="width:${p/max*100}%"></div></div><span>${(p*100).toFixed(2)}%</span></div>`).join('');
}

function renderOutput(){
 const mode=$('library-mode').value,paired=$('paired').checked,starts=[...selected].sort((a,b)=>a-b),recipes=starts.map(siteRecipe);
 const valid=recipes.every(Boolean),stats=valid?mixedLibraryStats(recipes,mode,paired):null;
 const formatNumber=n=>typeof n==='bigint'?n.toLocaleString():Number.isFinite(n)?n.toLocaleString():'Too many';
 const details=stats?(mode==='combinatorial'?`${formatNumber(stats.variants)} amino-acid combinations (excluding stops) · ${(stats.stopProbability*100).toFixed(2)}% probability of ≥1 stop at randomized sites.`:stats.perSite.map((d,i)=>`nt ${starts[i]+1}: ${formatNumber(d.variants)} amino acids, ${(d.stopProbability*100).toFixed(2)}% stops`).join(' · ')):'';
 $('library-summary').textContent=!valid?'Complete each position’s codon to calculate the library.':selected.size?`${selected.size} sites · ${formatNumber(stats.oligos)} oligos · ${stats.libraries} ${stats.libraries===1?'library':'separate libraries'} · ${details}`:'Select sites to preview library size and primer count.';
 starts.forEach((pos,i)=>{const el=$(`site-error-${pos}`);if(el)el.textContent=recipes[i]?'':'Enter three IUPAC symbols.';const input=document.querySelector(`[data-custom-site="${pos}"]`);if(input)input.setAttribute('aria-invalid',String(!recipes[i]));});
 renderChart();
 let error='';rows=[];
 if(sequence&&selected.size&&valid){try{rows=generate(sequence,frameId,starts,method,$('prefix').value,paired,customCodon(),mode,overrides);}catch(e){error=e.message;}}
 $('output-empty').textContent=error||(!valid?'Enter valid codons for all selected positions.':'Select one or more codons to generate your primer set.');
 $('output').hidden=!rows.length;$('output-empty').hidden=!!rows.length;
 $('copy').disabled=$('download').disabled=!rows.length;
 $('order-count').textContent=rows.length?`· ${rows.length} oligos`:'';
 if(!rows.length){$('rows').replaceChildren();return;}
 $('order-note').textContent=(mode==='separate'?'Each selected site is a separate library; all other positions retain the input sequence. Do not mix different site pools. ':'All selected sites are randomized together in one library. Component molar weights multiply across sites. ')+(paired?'Input and complement are separate strand pools. ':'')+'Codons and triplets are listed in ascending input-position order, including reverse frames. All sequences are 5′ → 3′.';
 $('rows').innerHTML=rows.map(r=>{
 const positions=r.starts.map(pos=>r.strand==='input'?pos:sequence.length-pos-3).sort((a,b)=>a-b);
 let highlighted='',cursor=0;for(const pos of positions){highlighted+=r.sequence.slice(cursor,pos)+`<mark>${r.sequence.slice(pos,pos+3)}</mark>`;cursor=pos+3;}highlighted+=r.sequence.slice(cursor);
 return `<tr><td>${escape(r.name)}</td><td>${escape(r.pool)}</td><td class="sequence">${highlighted}</td><td><code>${r.codon}</code></td><td>${r.ratio}</td><td>${(r.fraction*100).toPrecision(4)}%</td></tr>`;
 }).join('');
}
function exported(){return tsv(rows,$('format').value==='detailed',$('headers').checked);}
$('sequence').addEventListener('input',inputChanged);
$('example').addEventListener('click',()=>{$('sequence').value='CATGTTCATGGTGGCTGGACCGAATGGC';inputChanged();$('sequence').focus();});
$('frames').addEventListener('click',e=>{const b=e.target.closest('[data-frame]');if(!b)return;frameId=b.dataset.frame;manual=true;selected=new Set();overrides={};resetStatus();renderSelection();renderSiteControls();renderOutput();});
$('auto-frame').addEventListener('click',()=>{frameId=guessFrame(all);manual=false;selected=new Set();overrides={};resetStatus();renderSelection();renderSiteControls();renderOutput();});
$('codons').addEventListener('click',e=>{const b=e.target.closest('[data-start]');if(!b)return;const start=Number(b.dataset.start);if(selected.has(start)){selected.delete(start);delete overrides[start];}else selected.add(start);resetStatus();renderSelection();renderSiteControls();renderOutput();$('codons').querySelector(`[data-start="${start}"]`)?.focus({preventScroll:true});});
$('clear-selection').addEventListener('click',()=>{selected.clear();overrides={};resetStatus();renderSelection();renderSiteControls();renderOutput();});
$('library-mode').addEventListener('change',()=>{resetStatus();renderOutput();});
for(const b of document.querySelectorAll('[data-method]'))b.addEventListener('click',()=>{method=b.dataset.method;resetStatus();renderMethod();renderOutput();});
for(const id of ['other-codon','custom-codon'])$(id).addEventListener('input',()=>{resetStatus();renderMethod();renderOutput();});
for(const id of ['prefix','paired'])$(id).addEventListener('input',()=>{resetStatus();renderOutput();});
for(const id of ['format','headers'])$(id).addEventListener('change',resetStatus);
$('copy').addEventListener('click',async()=>{try{await navigator.clipboard.writeText(exported());$('status').textContent=`Copied ${rows.length} primers as TSV.`;}catch{const fallback=$('copy-fallback');fallback.hidden=false;fallback.value=exported();fallback.focus();fallback.select();$('status').textContent='Clipboard access is unavailable. The TSV is selected below; press Ctrl+C or ⌘C, or download the file.';}});
$('download').addEventListener('click',()=>{const blob=new Blob([exported()],{type:'text/tab-separated-values;charset=utf-8'});const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=`${rows[0].name}_primers.tsv`;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);$('status').textContent=`Downloaded ${rows.length} primers as TSV.`;});
renderMethod();inputChanged();

$('site-methods').addEventListener('change',e=>{
 if(!e.target.matches('[data-site]'))return;
 const pos=Number(e.target.dataset.site),value=e.target.value;
 if(value==='default')delete overrides[pos];
 else if(['22c','20c','NNK'].includes(value))overrides[pos]={method:value};
 else overrides[pos]={method:'other',custom:value==='custom'?'':value,preset:value};
 $('chart-site').value=String(pos);resetStatus();renderSiteControls();renderOutput();
});
$('site-methods').addEventListener('input',e=>{
 if(!e.target.matches('[data-custom-site]'))return;
 const pos=Number(e.target.dataset.customSite);overrides[pos].custom=e.target.value;
 $('chart-site').value=String(pos);resetStatus();renderOutput();
});
$('site-methods').addEventListener('click',e=>{if(e.target.id==='reset-methods'){overrides={};resetStatus();renderSiteControls();renderOutput();}});
$('chart-site').addEventListener('change',renderChart);
