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
 const titles={positions:'Normalized piston positions',heat:'Heat transfer: gas → wall',pressures:'Gas pressures',temperatures:'Gas and wall temperatures',flows:'Signed passage mass flows',velocity:'Normalized piston velocity per cycle radian',acceleration:'Normalized piston acceleration per cycle radian squared'};
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
   card.appendChild(note);if(!p.unavailable)for(const series of p.series)available.push({record:r,plot:p,series});
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
  const note=document.createElement('p');note.className='muted';note.textContent='Stored local frames in crank-radius units; no common-shaft layout inferred. Study angle before operation reversal. Animation speed is illustrative, not machine speed.';card.appendChild(note);
  const panels=[];
  const container=document.createElement('div');container.className='plots';card.appendChild(container);
  for(const side of ['small','large']){
   const data=motion.sides[side],panel=document.createElement('div');container.appendChild(panel);
   const label=document.createElement('h3');label.textContent=side+(data.crank_radius_m?' · crank radius '+data.crank_radius_m+' m':' · dimensionless geometry');panel.appendChild(label);
   if(!data.joints){const text=document.createElement('p');text.textContent='No linkage geometry for this motion family.';panel.appendChild(text);continue;}
   let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;
   for(const positions of Object.values(data.joints))for(const [x,y] of positions){xmin=Math.min(xmin,x);xmax=Math.max(xmax,x);ymin=Math.min(ymin,y);ymax=Math.max(ymax,y);}
   const scale=Math.min(520/(xmax-xmin||1),320/(ymax-ymin||1));
   const project=p=>[300+(p[0]-(xmin+xmax)/2)*scale,200-(p[1]-(ymin+ymax)/2)*scale];
   const drawing=document.createElement('div');panel.appendChild(drawing);panels.push({data,project,drawing});
  }
  if(!panels.length)continue;
  const button=document.createElement('button');button.textContent='Play';card.appendChild(button);
  const slider=document.createElement('input');slider.type='range';slider.min=0;slider.max=motion.angle.length-1;slider.value=0;slider.step=1;slider.setAttribute('aria-label','Mechanism study angle');card.appendChild(slider);
  const label=document.createElement('span');card.appendChild(label);
  let playing=false,last=null,frame=0;
  function draw(){const index=Number(slider.value);label.textContent=' '+(motion.angle[index]*180/Math.PI).toFixed(1)+'°';
   for(const {data,project,drawing} of panels){let svg='<svg role="img" aria-label="Linkage joint positions" viewBox="0 0 600 400">';
    for(const [a,b] of data.links){const p=project(data.joints[a][index]),q=project(data.joints[b][index]);svg+='<line x1="'+p[0]+'" y1="'+p[1]+'" x2="'+q[0]+'" y2="'+q[1]+'" stroke="#17658b" stroke-width="4"/>';}
    for(const [name,positions] of Object.entries(data.joints)){const p=project(positions[index]);svg+='<circle cx="'+p[0]+'" cy="'+p[1]+'" r="5" fill="#b74424"/><text x="'+(p[0]+8)+'" y="'+(p[1]-8)+'">'+esc(name)+'</text>';}
    drawing.innerHTML=svg+'</svg>';
   }
  }
  function tick(now){if(!playing)return;if(last!==null){frame=(frame+(now-last)*(motion.angle.length-1)/6000)%(motion.angle.length-1);slider.value=Math.floor(frame);draw();}last=now;requestAnimationFrame(tick);}
  button.onclick=()=>{playing=!playing;button.textContent=playing?'Pause':'Play';last=null;if(playing)requestAnimationFrame(tick);};
  slider.oninput=()=>{frame=Number(slider.value);draw();};draw();
 }
}
additionalPlots();mechanismPlots();
