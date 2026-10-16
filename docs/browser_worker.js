// The deployed site executes the original Python modules in an isolated worker.
// Pinned runtime: no server, approximation shortcut, or exact-formula substitute.
'use strict';
const runtimeURL = 'https://cdn.jsdelivr.net/pyodide/v0.28.3/full/';
let runtime;

async function initialize(job) {
  if (!runtime) runtime = (async () => {
    self.postMessage({type:'progress', job, message:'Loading the Python runtime. First use downloads the numerical libraries; please keep this tab open.'});
    importScripts(runtimeURL + 'pyodide.js');
    const py = await loadPyodide({indexURL:runtimeURL});
    await py.loadPackage(['numpy', 'scipy', 'sympy'], {
      messageCallback:message=>self.postMessage({type:'progress',job,message}),
    });
    py.FS.mkdirTree('/app');
    await Promise.all(['rosenau_hyman.py', 'porous_fin.py', 'dovetail_fin.py', 'solver_service.py',
      'validated_combinations.json', 'validated_fin_combinations.json', 'validated_dovetail_combinations.json'].map(async name=>{
      const url=new URL('python/'+name,self.location.href);
      const version=new URL(self.location.href).searchParams.get('v');
      if(version)url.searchParams.set('v',version);
      const response = await fetch(url);
      if (!response.ok) throw Error('Could not load solver file: ' + name);
      py.FS.writeFile('/app/' + name, await response.text());
    }));
    await py.runPythonAsync(`
import sys, json, time
sys.path.insert(0, '/app')
from solver_service import solve_request

def web_solve(parameters_json):
    parameters = json.loads(parameters_json)
    started = time.perf_counter()
    arrays, report = solve_request(parameters)
    report['execution_environment'] = 'Python / Pyodide 0.28.3 in browser Web Worker'
    report['elapsed_seconds'] = time.perf_counter() - started
    if not report['validation_passed']:
        return json.dumps({'report': report, 'error': 'Independent numerical validation failed. Previous result retained.'}, allow_nan=False)
    payload = dict(x=arrays['x'].tolist(), t=arrays['t'].tolist(), u=arrays['u'].T.tolist(),
                   exact=arrays['exact'].T.tolist() if 'exact' in arrays else None,
                   comparison=arrays['comparison'].T.tolist() if 'comparison' in arrays else None,
                   report=report)
    return json.dumps({'payload': payload}, allow_nan=False)
`);
    return py;
  })().catch(error=>{runtime=null; throw error;});
  return runtime;
}

self.onmessage = async ({data}) => {
  const {job, parameters} = data;
  try {
    const py = await initialize(job);
    self.postMessage({type:'progress',job,message:'Computing the '+(parameters.basis==='legendre'?'Legendre-wavelet':'Hermite')+' numerical solution and independent residual checks...'});
    py.globals.set('parameters_json', JSON.stringify(parameters));
    const result = JSON.parse(await py.runPythonAsync('web_solve(parameters_json)'));
    self.postMessage({type:'result',job,...result});
  } catch (error) {
    self.postMessage({type:'result',job,error:'Python computation or download failed: ' + error.message});
  }
};
