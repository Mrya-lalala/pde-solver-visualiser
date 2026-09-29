"""Numerical regressions and local API failure/validation checks (no dependencies)."""
import json
from pathlib import Path
import re
import tempfile
import threading
import time
import unittest
from unittest.mock import patch
from urllib.request import Request, urlopen
from urllib.error import HTTPError

import numpy as np
import local_server as app
import rosenau_hyman as rh


class NumericalTests(unittest.TestCase):
    def test_supported_matrix_and_reference_independence(self):
        for p in app.SUPPORTED:
            with self.subTest(p=p):
                arrays, report = app.solve_request(dict(p, reference=True))
                without, no_report = app.solve_request(dict(p, reference=False))
                self.assertTrue(report['validation_passed'])
                self.assertTrue(no_report['validation_passed'])
                self.assertLess(report['max_abs_solution_error'], 1e-6)
                self.assertNotIn('exact', without)
                for key in without:
                    np.testing.assert_array_equal(arrays[key], without[key])

    def test_original_and_speed_two_regression(self):
        arrays, report = app.solve_request(app.DEFAULTS)
        self.assertLess(report['max_abs_solution_error'], 3e-10)
        self.assertAlmostEqual(arrays['u'][0, 0], -8/3)
        restored, _, _ = rh.compute()
        np.testing.assert_array_equal(arrays['u'], restored['u'])
        changed, _ = app.solve_request(dict(app.DEFAULTS, speed=2))
        with np.load('results/solution.npz', allow_pickle=False) as old:
            np.testing.assert_array_equal(changed['u'], old['u'])

    def test_invalid_inputs_and_failed_combination(self):
        for change in ({'speed':float('nan')}, {'shift':float('inf')}, {'speed':True},
                       {'count':6.0}, {'length':0}, {'duration':2}, {'reference':'yes'},
                       {'extra':1}, {'speed':2,'shift':.5,'count':5}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                app.validate_request(dict(app.DEFAULTS, **change))
        f,g,e=rh.travelling_wave(2,.5)
        _,report,_=rh.compute(5,1,1,f,g,e)
        self.assertFalse(report['validation_passed'])


class ServiceTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.server=app.make_server(0,self.temp.name)
        self.thread=threading.Thread(target=self.server.serve_forever,daemon=True)
        self.thread.start()
        self.url=f'http://127.0.0.1:{self.server.server_port}'

    def tearDown(self):
        self.server.shutdown();self.server.server_close()
        self.server.app.worker.shutdown(wait=True)
        self.thread.join();self.temp.cleanup()

    def request(self,path,params=None,headers=None):
        options={} if params is None else dict(data=json.dumps(params).encode(),headers={
            'Content-Type':'application/json','Origin':self.url,
            'X-App-Token':self.server.app.token,**(headers or {})})
        with urlopen(Request(self.url+path,**options),timeout=10) as response:
            return json.load(response)

    def wait(self):
        deadline=time.monotonic()+15
        while time.monotonic()<deadline:
            state=self.request('/api/status')
            if state['status']!='running':return state
            time.sleep(.02)
        self.fail('Solve timed out')

    def test_real_solve_export_and_failure_retention(self):
        self.request('/api/solve',app.DEFAULTS)
        state=self.wait();self.assertEqual(state['status'],'complete')
        old=self.request('/api/result')
        folder=Path(self.temp.name)/state['run_id']
        with np.load(folder/'solution.npz') as arrays:
            np.testing.assert_array_equal(np.asarray(old['u']).T,arrays['u'])
        html=(folder/'interactive.html').read_text(encoding='utf-8')
        embedded=json.loads(re.search(r'const INITIAL=(.*?);\s*const \$=',html,re.S).group(1))
        self.assertEqual(embedded,old)
        with patch.object(app,'solve_request',side_effect=RuntimeError('Injected solver failure')):
            self.request('/api/solve',dict(app.DEFAULTS,speed=2))
            self.assertEqual(self.wait()['status'],'failed')
        self.assertEqual(self.request('/api/result'),old)
        arrays,report=app.solve_request(app.DEFAULTS);report['validation_passed']=False
        with patch.object(app,'solve_request',return_value=(arrays,report)):
            self.request('/api/solve',app.DEFAULTS)
            self.assertEqual(self.wait()['status'],'failed')
        self.assertEqual(self.request('/api/result'),old)

    def test_rejections_and_busy(self):
        for headers in ({'Origin':'http://example.com'},{'X-App-Token':'wrong'}):
            with self.assertRaises(HTTPError) as caught:self.request('/api/solve',app.DEFAULTS,headers)
            self.assertEqual(caught.exception.code,403)
        with self.assertRaises(HTTPError) as caught:self.request('/api/solve',dict(app.DEFAULTS,count=99))
        self.assertEqual(caught.exception.code,400)
        entered,release=threading.Event(),threading.Event()
        def blocked(p):
            entered.set();release.wait(5);raise RuntimeError('Test completion')
        with patch.object(app,'solve_request',side_effect=blocked):
            self.request('/api/solve',app.DEFAULTS);self.assertTrue(entered.wait(2))
            try:
                with self.assertRaises(HTTPError) as caught:self.request('/api/solve',app.DEFAULTS)
                self.assertEqual(caught.exception.code,409)
            finally:release.set()
            self.wait()


if __name__=='__main__':
    unittest.main()
