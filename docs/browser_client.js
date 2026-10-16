// Adapt the existing local API contract to a tab-local numerical worker.
// Each visitor has independent jobs/results; failed solves retain the last result.
'use strict';
window.createBrowserSolver = function(initial) {
  let state={status:'idle'}, latest=initial.saved_result, worker=null, timeout=null;
  const fail=message=>{
    if(state.status!=='running') return;
    clearTimeout(timeout);
    state={...state,status:'failed',error:message};
    worker?.terminate(); worker=null;
  };
  function start(parameters) {
    if(state.status==='running') throw Error('A solve is already running.');
    const dovetail=parameters.problem==='dovetail_fin';
    const combinations=dovetail?initial.dovetail_supported:parameters.problem==='porous_fin'?initial.fin_supported:initial.supported;
    const keys=dovetail
      ?['problem','basis','nc','nr','m_squared','taper','q','epsilon_g','theta_a','wet_parameter','power','angle_degrees','phi1','phi2','phi3','count','cells','duration','reference']
      :parameters.problem==='porous_fin'
      ?['problem','nc','nr','n1_squared','wet_term','power','duration','cells','count','reference']
      :['speed','shift','count','length','duration','reference'];
    if(Object.keys(parameters).length!==keys.length || keys.some(k=>!(k in parameters)) ||
       typeof parameters.reference!=='boolean' || (dovetail && parameters.basis!=='legendre') ||
       !combinations.some(p=>Object.keys(p).every(k=>p[k]===parameters[k]))) {
      throw Error('Unsupported parameter combination. Choose a tested preset.');
    }
    const job=crypto.randomUUID();
    state={status:'running',job_id:job,parameters};
    if(!worker) {
      try {const url=new URL('browser_worker.js',document.baseURI);if(initial.site_version)url.searchParams.set('v',initial.site_version);worker=new Worker(url);}
      catch(error) {fail(error.message); throw error;}
      worker.onerror=event=>{event.preventDefault();fail('Browser worker failed. Please retry or use a current desktop browser.');};
      worker.onmessage=({data})=>{
        if(data.job!==state.job_id || state.status!=='running') return;
        if(data.type==='progress') {
          document.getElementById('status').textContent=data.message+' Previous result remains available.';
          return;
        }
        clearTimeout(timeout);
        if(data.error) state={...state,status:'failed',error:data.error,report:data.report};
        else {
          latest=data.payload;
          state={...state,status:'complete',run_id:latest.report.run_id,report:latest.report};
        }
      };
    }
    // Includes first-use library downloads; a timeout only discards the new job.
    timeout=setTimeout(()=>fail('The solve timed out after five minutes. Previous result retained; retry on a faster connection or device.'),300000);
    worker.postMessage({job,parameters});
    return {job_id:job};
  }
  return async function(path,options) {
    if(path==='/api/status') return state;
    if(path==='/api/result' && latest) return latest;
    if(path==='/api/solve' && options?.method==='POST') return start(JSON.parse(options.body));
    throw Error('Result not available.');
  };
};
