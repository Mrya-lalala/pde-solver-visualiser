"""Independent mathematical checks of the new PDE and Legendre-wavelet method."""
import unittest
import json
from pathlib import Path
import tempfile
import threading
import time
from urllib.request import Request,urlopen
import numpy as np
import dovetail_fin as fin


class DovetailTests(unittest.TestCase):
    def test_supported_service_combinations_and_basis_rejection(self):
        import solver_service as service
        for row in service.DOVETAIL_COMBINATIONS:
            if not row['passed']:continue
            p=dict(row['parameters'],reference=True,problem='dovetail_fin',basis='legendre')
            with self.subTest(parameters=p):
                arrays,report=service.solve_request(p)
                self.assertTrue(report['validation_passed'])
                self.assertEqual(report['basis'],'legendre')
                self.assertNotIn('exact',arrays)
        default=dict(fin.DEFAULTS,problem='dovetail_fin',basis='legendre')
        for change in ({'basis':'hermite'},{'count':4},{'taper':-.2},{'extra':1}):
            with self.subTest(change=change),self.assertRaises(ValueError):service.validate_request(dict(default,**change))

    def test_local_api_runs_legendre_and_exports_arrays(self):
        import local_server as service
        with tempfile.TemporaryDirectory() as directory:
            server=service.make_server(0,directory)
            thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
            base=f'http://127.0.0.1:{server.server_port}'
            def get(path):
                with urlopen(base+path,timeout=10) as response:return json.load(response)
            try:
                p=dict(fin.DEFAULTS,problem='dovetail_fin',basis='legendre')
                request=Request(base+'/api/solve',data=json.dumps(p).encode(),headers={
                    'Content-Type':'application/json','Origin':base,'X-App-Token':server.app.token})
                with urlopen(request,timeout=10) as response:self.assertEqual(response.status,202)
                deadline=time.monotonic()+15
                while time.monotonic()<deadline:
                    status=get('/api/status')
                    if status['status']!='running':break
                    time.sleep(.02)
                self.assertEqual(status['status'],'complete')
                result=get('/api/result')
                folder=Path(directory)/status['run_id']
                with np.load(folder/'solution.npz') as data:np.testing.assert_array_equal(data['u'],np.asarray(result['u']).T)
                self.assertIn('Legendre',result['report']['method'])
                self.assertTrue((folder/'interactive.html').exists())
            finally:
                server.shutdown();server.server_close();server.app.worker.shutdown(wait=True);thread.join()

    def test_wavelets_are_local_orthonormal_and_constraints_hold(self):
        model=fin.LegendreWavelets(8,6,-.1)
        x,w=model.quadrature(12);raw=model.raw(x)[0]
        np.testing.assert_allclose(raw.T @ (w[:,None]*raw),np.eye(48),atol=4e-14)
        self.assertTrue(np.all(np.count_nonzero(raw,axis=1)<=6))
        v,d,_=model.operators(np.array([0.,1.]))
        np.testing.assert_allclose(d[0],0,atol=2e-12)
        np.testing.assert_allclose(v[1],0,atol=2e-14)
        for point in model.edges[1:-1]:
            left=model.operators(np.array([point-1e-12]))
            right=model.operators(np.array([point+1e-12]))
            np.testing.assert_allclose(left[0],right[0],atol=2e-9)
            np.testing.assert_allclose(left[1],right[1],atol=2e-6)

    def test_manufactured_solution_exercises_full_equation(self):
        p=dict(fin.DEFAULTS,duration=.4)
        model=fin.LegendreWavelets(p['cells'],p['count'],p['taper']);b=fin.coefficients(p)
        # theta=1+t*(1-X^2) satisfies theta_X(0,t)=0 and theta(1,t)=1.
        def source(x,t):
            a=1-p['taper']*x;theta=1+t*(1-x*x)
            diffusion=b['diffusion']*(-2*t*a+2*p['taper']*t*x)
            generation=b['generation']*a*(1+p['epsilon_g']*(theta-p['theta_a']))
            return (a*(1-x*x)-diffusion-generation+fin.losses(theta,p))/a
        result=model.integrate(p,initial=np.zeros_like(model.initial),source=source)
        x=np.linspace(0,1,101);t=np.linspace(0,.4,21)
        computed=1+model.operators(x)[0] @ result.sol(t)
        expected=1+(1-x[:,None]**2)*t
        self.assertLess(float(np.max(abs(computed-expected))),5e-8)

    def test_independent_analytic_limit_and_published_table(self):
        p=dict(fin.DEFAULTS,**fin.PRESETS['linear_benchmark']);p['duration']=5
        arrays,report=fin.compute(p)
        exact=fin.analytic_benchmark(arrays['x'],arrays['t'])
        check=arrays['t']>=.05
        self.assertTrue(report['validation_passed'])
        self.assertLess(float(np.max(abs(arrays['u'][:,check]-exact[:,check]))),2e-6)
        published=np.array([.648054,.661058,.700593,.768245,.866730,1])
        samples=np.interp(np.linspace(0,1,6),arrays['x'],arrays['u'][:,-1])
        # The paper's approximate tabulated values differ from the analytic
        # limit by up to 8e-7; they are not exact six-decimal targets.
        self.assertLess(float(np.max(abs(samples-published))),1e-6)

    def test_refinement_reference_independence_and_initial_data(self):
        a,r=fin.compute();b,s=fin.compute(dict(fin.DEFAULTS,reference=False))
        self.assertTrue(r['validation_passed']);self.assertTrue(s['validation_passed'])
        for key in b:np.testing.assert_array_equal(a[key],b[key])
        self.assertNotIn('comparison',b);self.assertNotIn('exact',a)
        np.testing.assert_array_equal(a['u'][:-1,0],0);self.assertEqual(a['u'][-1,0],1)
        raw=fin.LegendreWavelets(8,6,fin.DEFAULTS['taper']).raw(a['x'])[0]
        np.testing.assert_allclose(1+raw @ a['wavelet_coefficients'][:,1:],a['u'][:,1:],atol=2e-14)
        residuals=[]
        for cells in (4,8,16):
            _,rr=fin.compute(dict(fin.DEFAULTS,cells=cells,reference=False))
            residuals.append(rr['independent_weighted_strong_pde_residual'])
        self.assertTrue(all(a>b for a,b in zip(residuals,residuals[1:])))
        centres,reference=fin.finite_volume(fin.DEFAULTS,grid=800)
        values=reference.sol(a['t'])
        comparison=np.column_stack([np.interp(a['x'],np.r_[0,centres,1],np.r_[row[0],row,1]) for row in values.T])
        self.assertLess(float(np.max(abs(a['u'][:,a['t']>=.05]-comparison[:,a['t']>=.05]))),2e-6)

    def test_invalid_numeric_inputs_and_unity_material_limit(self):
        for change in ({'taper':1},{'count':6.0},{'power':.5},{'cells':True},
                       {'theta_a':1},{'nc':float('nan')},{'phi1':-.1},{'reference':'yes'}):
            with self.subTest(change=change),self.assertRaises(ValueError):
                fin.compute(dict(fin.DEFAULTS,**change))
        ratios=fin.material_ratios(dict(fin.DEFAULTS,phi1=0,phi2=0,phi3=0))
        for value in ratios.values():self.assertAlmostEqual(value,1)


if __name__=='__main__':unittest.main()
