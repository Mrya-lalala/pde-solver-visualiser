"""Independent accuracy, boundary and API checks for the PDF's fin PDE."""
import json
from pathlib import Path
import re
import tempfile
import threading
import time
import unittest
from urllib.request import Request,urlopen
from urllib.error import HTTPError

import numpy as np
import porous_fin as fin
import local_server as server_module


class FinTests(unittest.TestCase):
    def test_manufactured_solution(self):
        # A compatible smooth problem with a known solution and derived forcing
        # exercises signs, spatial derivatives, constraints and BDF time accuracy.
        p=dict(fin.DEFAULTS,cells=8,duration=.4)
        model=fin.HermiteFin(p['cells'],p['count'])
        def source(x,t):
            h=x*(2-x);hp=2-2*x
            return h+2*t+2*x/(1-x*x)*t*hp+fin.reaction(1+t*h,p)/np.sqrt(1-x*x)
        result=model.integrate(p,initial=np.zeros_like(model.initial),source=source)
        xs=np.linspace(0,1,101);ts=np.linspace(0,.4,21)
        numerical=1+model.operators(xs)[0] @ result.sol(ts)
        expected=1+xs[:,None]*(2-xs[:,None])*ts
        self.assertLess(float(np.max(np.abs(numerical-expected))),2e-7)

    def test_supported_matrix(self):
        for row in server_module.FIN_COMBINATIONS:
            p=dict(row['parameters'],reference=True)
            with self.subTest(p=p):
                arrays,report=fin.compute(p)
                self.assertEqual(report['validation_passed'],row['passed'])
                if row['passed']:
                    self.assertLess(report['max_abs_comparison_difference_after_start'],2e-5)
                    self.assertLess(report['base_max_abs_error'],1e-10)
                    self.assertLess(report['tip_max_abs_derivative'],1e-10)
                self.assertNotIn('exact',arrays)

    def test_refinement_reference_and_startup(self):
        a,r=fin.compute()
        b,s=fin.compute(dict(fin.DEFAULTS,reference=False))
        for key in b:np.testing.assert_array_equal(a[key],b[key])
        self.assertNotIn('comparison',b)
        self.assertTrue(s['validation_passed'])
        np.testing.assert_array_equal(a['u'][1:,0],0)
        self.assertEqual(a['u'][0,0],1)
        residuals=[]
        for cells in (8,12,16):
            _,rr=fin.compute(dict(fin.DEFAULTS,cells=cells,reference=False))
            residuals.append(rr['independent_weighted_strong_pde_residual'])
        self.assertTrue(all(a>b for a,b in zip(residuals,residuals[1:])))
        centres,fine=fin.finite_volume(fin.DEFAULTS,grid=800)
        x,t=a['x'],a['t'];values=fine.sol(t)
        comparison=np.column_stack([np.interp(x,np.r_[0,centres,1],np.r_[1,u,u[-1]]) for u in values.T])
        self.assertLess(float(np.max(abs(a['u'][:,t>=.05]-comparison[:,t>=.05]))),1e-5)
        self.assertTrue(np.isfinite(a['efficiency']).all())

    def test_invalid_and_unsupported(self):
        for change in ({'cells':12.0},{'power':1.5},{'duration':0}, {'nr':float('nan')},
                       {'reference':'yes'},{'nc':True},{'wet_term':-1}):
            with self.subTest(change=change),self.assertRaises(ValueError):
                fin.compute(dict(fin.DEFAULTS,**change))
        with self.assertRaises(ValueError):
            server_module.validate_request(dict(fin.DEFAULTS,problem='porous_fin',count=4))

    def test_api_solve_and_offline_export(self):
        with tempfile.TemporaryDirectory() as folder:
            server=server_module.make_server(0,folder)
            thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
            base=f'http://127.0.0.1:{server.server_port}'
            def get(path):
                with urlopen(base+path,timeout=10) as response:return json.load(response)
            try:
                p=dict(fin.DEFAULTS,problem='porous_fin')
                req=Request(base+'/api/solve',data=json.dumps(p).encode(),headers={
                    'Content-Type':'application/json','Origin':base,'X-App-Token':server.app.token})
                with urlopen(req,timeout=10) as response:self.assertEqual(response.status,202)
                for _ in range(300):
                    state=get('/api/status')
                    if state['status']!='running':break
                    time.sleep(.02)
                self.assertEqual(state['status'],'complete')
                payload=get('/api/result')
                self.assertIsNone(payload['exact'])
                self.assertIn('not exact',payload['report']['comparison_label'])
                saved=Path(folder)/state['run_id']
                with np.load(saved/'solution.npz') as data:
                    np.testing.assert_array_equal(np.asarray(payload['u']).T,data['u'])
                    np.testing.assert_array_equal(np.asarray(payload['comparison']).T,data['comparison'])
                html=(saved/'interactive.html').read_text(encoding='utf-8')
                embedded=json.loads(re.search(r'const INITIAL=(.*?);\s*const \$=',html,re.S).group(1))
                self.assertEqual(embedded,payload)
            finally:
                server.shutdown();server.server_close();server.app.worker.shutdown();thread.join()


if __name__=='__main__':unittest.main()
