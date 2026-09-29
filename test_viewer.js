// State/interaction regression checks with a small DOM/Plotly test double.
// This checks application logic, not real browser layout or Plotly rendering.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const html = fs.readFileSync('viewer_template.html', 'utf8');
const defaults = {speed:1,shift:0,count:6,length:1,duration:1,reference:true};
const supported = JSON.parse(fs.readFileSync('validated_combinations.json','utf8')).filter(r=>r.passed).map(r=>r.parameters);
const elements = new Map();
function element(id) {
  if (!elements.has(id)) elements.set(id, {value:'',checked:true,hidden:false,open:false,
    textContent:'',disabled:false,setAttribute(k,v){this[k]=v},removeAllListeners(){},on(){}});
  return elements.get(id);
}
let result = null, status = {status:'idle'}, requests = [], fail = false;
const context = vm.createContext({console, document:{getElementById:element,addEventListener(){}},
  window:{},setTimeout,requestAnimationFrame(){},
  Plotly:{async react(el,data,layout){el.data=data;el.layout=layout},async restyle(el,update){el.update=update},async relayout(el,update){el.layoutUpdate=update}},
  fetch:async(path,options)=>{
    requests.push(path);
    if(path==='/api/solve') {
      const p=JSON.parse(options.body);
      status=fail?{status:'failed',job_id:'job',parameters:p,error:'test failure'}:{status:'complete',job_id:'job',parameters:p,run_id:'run'};
      if(!fail) result={x:[0,1],t:[0,.5,1],u:[[-2,-1],[-1,-.5],[0,.1]],exact:p.reference?[[-2,-1],[-1,-.5],[0,.1]]:null,report:{parameters:p,run_id:'run'}};
      return {ok:true,json:async()=>({job_id:'job'})};
    }
    return {ok:path==='/api/status'||!!result,json:async()=>path==='/api/status'?status:result||{error:'Not found'}};
  }});
const script = html.split('<script>').at(-1).split('</script>')[0].replace('/*__DATA__*/',JSON.stringify({live:true,token:'test',defaults,supported}));
const flush=()=>new Promise(resolve=>setImmediate(resolve));
(async()=>{
  vm.runInContext(script,context);await flush();
  assert.equal(element('wave-speed').value,1);
  await element('solve').onclick();await flush();
  assert.match(element('status').textContent,/Solve passed/);
  assert.equal(element('curve').data.length,2);
  const before=requests.length;
  element('time').value='2';element('time').oninput();await flush();
  assert.equal(requests.length,before,'time scrubbing must not call backend');
  assert.match(element('time-value').textContent,/1.0000/);
  element('overlay').checked=false;await element('overlay').onchange();
  assert.equal(element('curve').data[1].visible,false);
  element('wave-speed').value='2';element('shift').value='.5';element('count').value='5';element('duration').value='1';element('count').onchange();
  assert.equal(element('solve').disabled,true);
  const displayed=element('displayed').textContent;
  element('original').onclick();assert.equal(element('solve').disabled,false);
  assert.equal(element('displayed').textContent,displayed,'reset must not relabel the old result');
  fail=true;await element('solve').onclick();await flush();
  assert.match(element('status').textContent,/failed/);
  assert.equal(element('displayed').textContent,displayed);
  fail=false;element('reference').checked=false;
  await element('solve').onclick();await flush();
  assert.equal(element('curve').data.length,1);
  assert.equal(element('error').disabled,true);
  assert.equal(element('overlay').disabled,true);
  console.log('Viewer state checks passed (no real-browser rendering).');
})().catch(e=>{console.error(e);process.exitCode=1});
