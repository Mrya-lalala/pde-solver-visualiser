// Browser adapter behavior: strict inputs, job isolation, failure/timeout retention.
const assert=require('node:assert/strict');
const fs=require('node:fs');
const vm=require('node:vm');
const workers=[],timers=new Map();let next=0;
const context=vm.createContext({window:{},document:{baseURI:'https://example.test/pde/',
  getElementById:()=>({textContent:''})},URL,
  crypto:{randomUUID:()=>String(++next)},
  setTimeout:fn=>{const id=++next;timers.set(id,fn);return id;},
  clearTimeout:id=>timers.delete(id),
  Worker:class {
    constructor(){workers.push(this);}
    postMessage(message){this.message=message;}
    terminate(){this.terminated=true;}
  }});
vm.runInContext(fs.readFileSync('browser_client.js','utf8'),context);
const rh={speed:1,shift:0,count:6,length:1,duration:1};
const included={report:{run_id:'saved'}};
const initial={saved_result:included,supported:[rh],fin_supported:[]};
const create=()=>context.window.createBrowserSolver(initial);
const submit=(api,p={...rh,reference:true})=>api('/api/solve',{method:'POST',body:JSON.stringify(p)});

(async()=>{
  const api=create(),other=create();
  await assert.rejects(submit(api,{...rh,reference:true,count:4}),/Unsupported/);
  await assert.rejects(submit(api,{...rh,reference:true,extra:1}),/Unsupported/);
  assert.equal(workers.length,0,'rejected inputs must not start a numerical worker');
  const job=await submit(api);
  await assert.rejects(submit(api),/already running/);
  assert.equal((await other('/api/status')).status,'idle','visitors must have separate jobs');
  const computed={report:{run_id:'computed'}};
  workers[0].onmessage({data:{job:job.job_id,type:'result',payload:computed}});
  assert.equal((await api('/api/status')).status,'complete');
  assert.equal(await api('/api/result'),computed);
  assert.equal(await other('/api/result'),included);
  const failed=await submit(api);
  workers[0].onmessage({data:{job:failed.job_id,type:'result',error:'Independent validation failed'}});
  assert.equal((await api('/api/status')).status,'failed');
  assert.equal(await api('/api/result'),computed,'failed solve must retain the last successful result');
  const retry=await submit(api);
  workers[0].onmessage({data:{job:job.job_id,type:'result',payload:{report:{run_id:'stale'}}}});
  assert.equal((await api('/api/status')).job_id,retry.job_id);
  assert.equal((await api('/api/status')).status,'running','late messages must not replace a new job');
  [...timers.values()][0]();
  assert.equal((await api('/api/status')).status,'failed');
  assert.equal(await api('/api/result'),computed,'timeout must retain the last successful result');
  assert.equal(workers[0].terminated,true,'timed-out work must stop');
  await submit(api);
  assert.equal(workers.length,2,'retry after timeout must use a fresh worker');
  console.log('Browser adapter checks passed.');
})().catch(error=>{console.error(error);process.exitCode=1;});
