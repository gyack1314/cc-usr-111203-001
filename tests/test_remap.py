"""球面守恒重映射核对:跨日界线、整圈经度、极点、沿岸掩码变化、共享源格与守恒。"""
import json
import math
import subprocess
import sys
import unittest
from climate_data.domain import cell_area
from climate_data.io import load
from climate_data.remap import remap

H10 = 2 * math.sin(math.radians(10))  # -10..10纬度带的sin差


def cell(cid, west, width, south=-10.0, north=10.0):
    return {'id': cid, 'west': west, 'width': width, 'south': south, 'north': north}


def frame(fid, values, valid, variance):
    return {'id': fid, 'values': values, 'valid': valid, 'variance': variance}


def dataset(source, target, frames, quantity='intensive', min_coverage=0.8, pairs=()):
    return {'radius': 1.0, 'source': source, 'target': target, 'frames': frames,
            'quantity': quantity, 'min_coverage': min_coverage,
            'covariance_pairs': [list(p) for p in pairs], 'geometry_tolerance': 1e-12}


class GeometryTests(unittest.TestCase):
    def test_dateline_crossing(self):
        v = dataset([cell('s0', 350, 20)],
                    [cell('t0', -10, 20), cell('t1', 10, 10)],
                    [frame('f', {'s0': 5}, {'s0': True}, {'s0': 2})], min_coverage=1.0)
        r = remap(v)
        overlaps = r['geometry']['overlaps']
        self.assertEqual(len(overlaps), 1)  # t1仅接边界,零面积不进稀疏结构
        self.assertEqual((overlaps[0]['source'], overlaps[0]['target']), ('s0', 't0'))
        self.assertAlmostEqual(overlaps[0]['area'], cell_area(cell('s0', 350, 20), 1.0), places=15)
        t = r['frames'][0]['targets']
        self.assertEqual(t['t0']['mean'], 5)  # west=-10与west=350是同一格
        self.assertIsNone(t['t1']['mean'])
        self.assertEqual(t['t1']['variance_status'], 'no_valid_coverage')

    def test_full_circle_polar_ring(self):
        r = remap(load('data/full-polar-ring.json')[0])
        expected = 2 * math.pi * (1 - math.sin(math.radians(80)))
        self.assertAlmostEqual(r['geometry']['overlaps'][0]['area'], expected, places=15)
        t = r['frames'][0]['targets']['polar']
        self.assertEqual(t['mean'], 3)
        self.assertEqual(t['variance'], 2)
        self.assertAlmostEqual(t['coverage_ratio'], 1.0, places=15)

    def test_polar_area_exact_not_planar(self):
        area = cell_area(cell('p', 0, 360, 80, 90), 1.0)
        self.assertAlmostEqual(area, 2 * math.pi * (1 - math.sin(math.radians(80))), places=15)
        self.assertGreater(area, 0)  # 极点处零长度经边不等于零面积格

    def test_edge_touch_zero_area(self):
        v = dataset([cell('s0', 170, 20)], [cell('t0', 190, 20)],
                    [frame('f', {'s0': 1}, {'s0': True}, {'s0': 1})])
        r = remap(v)
        self.assertEqual(r['geometry']['overlaps'], [])
        t = r['frames'][0]['targets']['t0']
        self.assertIsNone(t['mean'])
        self.assertEqual(t['coverage_integral'], 0)
        self.assertAlmostEqual(r['frames'][0]['conservation']['residual'], 0.0, places=15)

    def test_self_overlap_rejected(self):
        with self.assertRaises(ValueError):  # 源格互相重叠
            remap(dataset([cell('s0', 170, 20), cell('s1', 180, 20)], [cell('t0', 170, 40)],
                          [frame('f', {'s0': 1, 's1': 2}, {'s0': True, 's1': True}, {'s0': 1, 's1': 1})]))
        with self.assertRaises(ValueError):  # 目标格互相重叠
            remap(dataset([cell('s0', 170, 40)], [cell('t0', 170, 20), cell('t1', 180, 20)],
                          [frame('f', {'s0': 1}, {'s0': True}, {'s0': 1})]))
        with self.assertRaises(ValueError):  # 整圈经度格与任何其他源格重叠
            remap(dataset([cell('s0', 0, 360), cell('s1', 10, 10)], [cell('t0', 0, 360)],
                          [frame('f', {'s0': 1, 's1': 2}, {'s0': True, 's1': True}, {'s0': 1, 's1': 1})]))

    def test_sparse_geometry(self):
        r = remap(load('data/normal.json')[0])
        keys = {(e['source'], e['target']) for e in r['geometry']['overlaps']}
        self.assertEqual(keys, {('s0', 't0'), ('s0', 't1'), ('s1', 't1')})  # 无(s1,t0),不是完整矩阵


class IntensiveTests(unittest.TestCase):
    def test_normal_frames_and_mask_change(self):
        r = remap(load('data/normal.json')[0])
        a10 = math.radians(10) * H10
        day1, day2 = r['frames']
        t0 = day1['targets']['t0']
        self.assertAlmostEqual(t0['mean'], 10.0, places=12)
        self.assertAlmostEqual(t0['variance'], 4.0, places=12)
        t1 = day1['targets']['t1']
        self.assertAlmostEqual(t1['mean'], 50 / 3, places=12)
        self.assertAlmostEqual(t1['variance'], 40 / 9, places=12)
        self.assertAlmostEqual(t1['valid_coverage_area'], 3 * a10, places=15)
        weights = {c['source']: c['weight'] for c in t1['contributions']}
        self.assertAlmostEqual(weights['s0'], 1 / 3, places=12)
        self.assertAlmostEqual(weights['s1'], 2 / 3, places=12)
        cov = day1['covariances'][0]
        self.assertEqual(cov['status'], 'ok')
        self.assertAlmostEqual(cov['value'], 4 / 3, places=12)  # 共享s0,协方差非对角
        # 沿岸掩码变化:day2的s1无效,分母必须按当前帧重算
        t1d2 = day2['targets']['t1']
        self.assertAlmostEqual(t1d2['valid_coverage_area'], a10, places=15)
        self.assertIsNone(t1d2['mean'])
        self.assertEqual(t1d2['variance_status'], 'below_min_coverage')
        self.assertAlmostEqual(t1d2['coverage_integral'], 10 * a10, places=12)
        self.assertIsNone(day2['covariances'][0]['value'])
        self.assertEqual(day2['covariances'][0]['status'], 'below_min_coverage')
        for f in r['frames']:
            self.assertAlmostEqual(f['conservation']['residual'], 0.0, places=12)
            self.assertEqual(f['conservation']['domain'], 'overlap')

    def test_shared_source_covariance(self):
        # 源格s0横跨两个相邻目标,目标互不重叠;两目标各有一个专属源格
        v = dataset([cell('s0', 5, 10, 0, 10), cell('s1', 0, 5, 0, 10), cell('s2', 15, 5, 0, 10)],
                    [cell('t0', 0, 10, 0, 10), cell('t1', 10, 10, 0, 10)],
                    [frame('f', {'s0': 2, 's1': 1, 's2': 4},
                           {'s0': True, 's1': True, 's2': True},
                           {'s0': 2, 's1': 1, 's2': 3})],
                    min_coverage=1.0, pairs=[('t0', 't1')])
        f = remap(v)['frames'][0]
        t = f['targets']
        self.assertAlmostEqual(t['t0']['mean'], 1.5, places=12)
        self.assertAlmostEqual(t['t1']['mean'], 3.0, places=12)
        self.assertAlmostEqual(t['t0']['variance'], 0.75, places=12)
        self.assertAlmostEqual(t['t1']['variance'], 1.25, places=12)
        # 共享s0:Cov = (1/2)(1/2)*2 = 0.5,不是只对角
        self.assertAlmostEqual(f['covariances'][0]['value'], 0.5, places=12)

    def test_unknown_variance_not_estimable(self):
        r = remap(load('data/unknown-variance.json')[0])
        day1, day2 = r['frames']
        t0 = day1['targets']['t0']
        self.assertAlmostEqual(t0['mean'], 10.0, places=12)  # 均值照出
        self.assertIsNone(t0['variance'])  # 不确定性不可估,不补零
        self.assertEqual(t0['variance_status'], 'unknown_source_variance')
        self.assertIsNone(day1['covariances'][0]['value'])
        self.assertEqual(day1['covariances'][0]['status'], 'unknown_source_variance')
        self.assertAlmostEqual(day2['targets']['t0']['variance'], 4.0, places=12)
        self.assertEqual(day2['targets']['t0']['variance_status'], 'ok')

    def test_dateline_file(self):
        r = remap(load('data/dateline.json')[0])
        day1, day2 = r['frames']
        self.assertAlmostEqual(day1['targets']['t0']['coverage_ratio'], 0.75, places=12)
        self.assertAlmostEqual(day1['targets']['t0']['mean'], 6.0, places=12)
        self.assertAlmostEqual(day1['targets']['t1']['mean'], 10.5, places=12)
        self.assertAlmostEqual(day1['targets']['t1']['variance'], 1.1875, places=12)
        self.assertAlmostEqual(day1['covariances'][0]['value'], 0.25, places=12)
        self.assertIsNone(day2['targets']['t0']['mean'])
        self.assertEqual(day2['targets']['t0']['variance_status'], 'no_valid_coverage')
        self.assertAlmostEqual(day2['targets']['t1']['mean'], 12.0, places=12)
        self.assertAlmostEqual(day2['targets']['t1']['variance'], 2.0, places=12)
        self.assertIsNone(day2['covariances'][0]['value'])
        for f in r['frames']:
            self.assertAlmostEqual(f['conservation']['residual'], 0.0, places=12)


class ExtensiveTests(unittest.TestCase):
    def test_source_totals_file(self):
        r = remap(load('data/source-totals.json')[0])
        day1, day2 = r['frames']
        self.assertAlmostEqual(day1['targets']['t0']['allocated_total'], 5.0, places=12)
        self.assertAlmostEqual(day1['targets']['t0']['variance'], 1.0, places=12)
        self.assertAlmostEqual(day1['targets']['t1']['allocated_total'], 25.0, places=12)
        self.assertAlmostEqual(day1['targets']['t1']['variance'], 10.0, places=12)
        self.assertAlmostEqual(day1['covariances'][0]['value'], 1.0, places=12)  # 共享s0
        self.assertEqual(set(day1['unallocated_by_source']), {'s0', 's1'})
        for x in day1['unallocated_by_source'].values():
            self.assertAlmostEqual(x, 0.0, places=12)
        c = day1['conservation']
        self.assertAlmostEqual(c['valid_source_total'], 30.0, places=12)
        self.assertAlmostEqual(c['residual'], 0.0, places=12)
        # day2:s1无效,其总量不参与分配,也不挤进已覆盖区域
        self.assertAlmostEqual(day2['targets']['t1']['allocated_total'], 5.0, places=12)
        self.assertAlmostEqual(day2['conservation']['valid_source_total'], 10.0, places=12)
        self.assertNotIn('s1', day2['unallocated_by_source'])
        self.assertAlmostEqual(day2['conservation']['residual'], 0.0, places=12)

    def test_unallocated_share_not_squeezed(self):
        v = dataset([cell('s0', 0, 20, 0, 10)], [cell('t0', 0, 10, 0, 10)],
                    [frame('f', {'s0': 10}, {'s0': True}, {'s0': 4})], quantity='extensive')
        f = remap(v)['frames'][0]
        self.assertAlmostEqual(f['targets']['t0']['allocated_total'], 5.0, places=12)  # 不再归一化到10
        self.assertAlmostEqual(f['targets']['t0']['variance'], 1.0, places=12)
        self.assertAlmostEqual(f['unallocated_by_source']['s0'], 5.0, places=12)
        c = f['conservation']
        self.assertAlmostEqual(c['allocated_total'] + c['unallocated_total'],
                               c['valid_source_total'], places=12)
        self.assertAlmostEqual(c['residual'], 0.0, places=12)

    def test_extensive_unknown_variance(self):
        v = dataset([cell('s0', 0, 20, 0, 10)], [cell('t0', 0, 10, 0, 10)],
                    [frame('f', {'s0': 10}, {'s0': True}, {'s0': None})], quantity='extensive')
        t = remap(v)['frames'][0]['targets']['t0']
        self.assertAlmostEqual(t['allocated_total'], 5.0, places=12)
        self.assertIsNone(t['variance'])
        self.assertEqual(t['variance_status'], 'unknown_source_variance')


class ConfigTests(unittest.TestCase):
    def test_invalid_configs(self):
        fr = [frame('f', {'s0': 1}, {'s0': True}, {'s0': 1})]
        base = dataset([cell('s0', 0, 10)], [cell('t0', 0, 10)], fr)
        with self.assertRaises(ValueError):
            remap(dict(base, quantity='average'))
        with self.assertRaises(ValueError):
            remap(dict(base, min_coverage=1.5))
        with self.assertRaises(ValueError):
            remap(dict(base, covariance_pairs=[['t0', 't9']]))
        with self.assertRaises(ValueError):
            remap(dict(base, covariance_pairs=[['t0', 't0']]))

    def test_invalid_domain_file(self):
        with self.assertRaises(ValueError):
            remap(load('data/invalid-domain.json')[0])

    def test_duplicate_pair_deduped(self):
        v = dataset([cell('s0', 0, 10)], [cell('t0', 0, 5), cell('t1', 5, 5)],
                    [frame('f', {'s0': 1}, {'s0': True}, {'s0': 1})],
                    pairs=[('t0', 't1'), ('t1', 't0')])
        self.assertEqual(len(remap(v)['frames'][0]['covariances']), 1)


class CliTests(unittest.TestCase):
    def test_remap_command(self):
        p = subprocess.run([sys.executable, '-m', 'climate_data.cli', 'remap', 'data/normal.json'],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        out = json.loads(p.stdout)
        self.assertEqual(len(out['source_sha256']), 64)
        self.assertEqual(len(out['remap']['frames']), 2)
        self.assertEqual(out['remap']['quantity'], 'intensive')


if __name__ == '__main__':
    unittest.main()
