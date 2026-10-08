/* Tick layout only; scientific curve values are never rounded. */
function niceAxis(minimum,maximum,includeZero=false){
 if(includeZero){minimum=Math.min(0,minimum);maximum=Math.max(0,maximum);}
 if(minimum===maximum){const pad=Math.abs(minimum)*.1||1;minimum-=pad;maximum+=pad;}
 const raw=(maximum-minimum)/5,base=Math.pow(10,Math.floor(Math.log10(raw))),fraction=raw/base;
 const step=(fraction<=1?1:fraction<=2?2:fraction<=5?5:10)*base;
 const low=Math.floor(minimum/step)*step,high=Math.ceil(maximum/step)*step,ticks=[];
 for(let i=0;i<=Math.round((high-low)/step);i++)ticks.push(low+i*step);
 return {low,high,step,ticks};
}
function tickLabel(value,step){
 const decimals=Math.max(0,-Math.floor(Math.log10(step)));
 return (Math.abs(value)<step*1e-8?0:value).toLocaleString('en-US',{maximumFractionDigits:Math.min(decimals,12),useGrouping:false});
}
/* Offline rendering only: all scientific samples come from production Python. */
function additionalPlots(){
 const host=byId('additionalPlots');host.innerHTML='';
 const titles={positions:'Normalized piston positions',heat:'Heat transfer: gas → wall',pressures:'Gas pressures',temperatures:'Cylinder gas temperatures',flows:'Signed passage mass flows',velocity:'Normalized piston velocity per cycle radian',acceleration:'Normalized piston acceleration per cycle radian squared'};
 const colors=['#17658b','#b74424','#087660','#773baa','#986700','#b32573','#247070','#666666'];
 for(const [kind,title] of Object.entries(titles)){
  const candidates=rows.filter(r=>r.plots?.[kind]);if(!candidates.length)continue;
  const card=document.createElement('section');card.className='card';host.appendChild(card);
  const heading=document.createElement('h2');heading.textContent=title;card.appendChild(heading);
  const available=[];
  for(const r of candidates){
   const p=r.plots[kind],note=document.createElement('p');note.className=p.unavailable?'warning':'muted';
   note.textContent=r.candidate_id.slice(0,12)+' · '+(p.unavailable?'Unavailable: '+p.unavailable:p.method)+(r.plots_runtime_compatible===false?' Stored runtime differs: these curves use the current code.':'');
   if(p.replay)note.textContent+=' Endpoint relative drift: '+p.replay.endpoint_relative_drift.toPrecision(3)+'. Backend: '+(p.replay.backend.actual_backend||'recorded in plot data')+'.';
   card.appendChild(note);if(!p.unavailable)for(const series of p.series){if(kind==='temperatures'&&!['small gas','large gas'].includes(series.label))continue;available.push({record:r,plot:p,series});}
  }
  if(!available.length)continue;
  let low=Infinity,high=-Infinity;for(const a of available)for(const v of a.series.values){low=Math.min(low,v);high=Math.max(high,v);}
  if(!Number.isFinite(low)||!Number.isFinite(high))continue;
  const axis=niceAxis(low,high,kind==='heat'||kind==='flows'||kind==='velocity'||kind==='acceleration');low=axis.low;high=axis.high;
  const w=1000,h=400,l=100,b=55,t=20,right=20,x=v=>l+v/360*(w-l-right),y=v=>h-b-(v-low)/(high-low)*(h-b-t);
  let svg='<svg role="img" aria-label="'+esc(title)+'" viewBox="0 0 '+w+' '+h+'"><path d="M'+l+' '+t+' V'+(h-b)+' H'+(w-right)+'" fill="none" stroke="#526575"/>';
  for(let i=0;i<=6;i++)svg+='<text x="'+x(i*60)+'" y="'+(h-b+22)+'" text-anchor="middle" font-size="13">'+i*60+'</text>';
  for(const v of axis.ticks)svg+='<text x="'+(l-8)+'" y="'+(y(v)+4)+'" text-anchor="end" font-size="13">'+tickLabel(v,axis.step)+'</text>';
  if(low<=0&&high>=0)svg+='<line class="zero-line" x1="'+l+'" x2="'+(w-right)+'" y1="'+y(0)+'" y2="'+y(0)+'" stroke="#526575" stroke-width="1.5" stroke-dasharray="4 3"/>';
  for(const [i,a] of available.entries()){
   const label=a.record.candidate_id.slice(0,12)+' · '+a.series.label;
   const points=a.plot.angle.map((angle,j)=>x(angle)+','+y(a.series.values[j])).join(' ');
   svg+='<polyline points="'+points+'" fill="none" stroke="'+colors[i%colors.length]+'" stroke-width="2"'+(candidates.indexOf(a.record)%2?' stroke-dasharray="7 4"':'')+'><title>'+esc(label)+'</title></polyline>';
  }
  svg+='<text x="550" y="395" text-anchor="middle">Solver cycle angle [deg]</text><text transform="translate(20 190) rotate(-90)" text-anchor="middle">'+esc(available[0].plot.unit)+'</text></svg>';
  const drawing=document.createElement('div');drawing.innerHTML=svg;card.appendChild(drawing);
  const legend=document.createElement('ul');for(const [i,a] of available.entries()){const li=document.createElement('li');li.style.color=colors[i%colors.length];li.textContent=a.record.candidate_id.slice(0,12)+' · '+a.series.label+(candidates.indexOf(a.record)%2?' (dashed)':' (solid)');legend.appendChild(li);}card.appendChild(legend);
 }
}
function mechanismPlots(){
 const host=byId('mechanismPlots');
 for(const row of rows){
  const motion=row.plots?.mechanisms;if(!motion)continue;
  const card=document.createElement('section');card.className='card';host.appendChild(card);
  const heading=document.createElement('h2');heading.textContent='Mechanisms · '+row.candidate_id.slice(0,12);card.appendChild(heading);
  const note=document.createElement('p');note.className=motion.unavailable?'warning':'muted';
  note.textContent=motion.unavailable?'Unavailable: '+motion.unavailable:motion.method+' Animation speed is illustrative. Layout: '+motion.layout+'.';
  if(row.plots_runtime_compatible===false)note.textContent+=' Stored runtime differs; this animation uses the current code.';
  card.appendChild(note);
  if(motion.data_uri || motion.url){
   const image=document.createElement('img');image.src=motion.data_uri || motion.url;
   image.alt='Whole machine with replayed gas temperatures and effective valve states';
   image.style.width='100%';image.style.height='auto';card.appendChild(image);
  }
 }
}
additionalPlots();mechanismPlots();
