"""球面守恒重映射:稀疏交叠权重、格内平均与总量分配、不确定性传播。

几何关系只保存正面积交叠(稀疏结构),不分配完整源乘目标矩阵;几何权重跨帧复用,
但有效掩码逐帧变化,分母按当前帧掩码重算,不缓存上一时刻的结果。
"""
import math
from .domain import cell_area, check_frames
from .io import number, positive, unique

QUANTITIES = ('intensive', 'extensive')
DEFAULT_TOLERANCE = 1e-12


def _lon_intervals(west, width):
    """经度规范到[0,360)的区间列表;width=360为整圈,允许跨日界线。"""
    if width >= 360.0:
        return [(0.0, 360.0)]
    w = math.fmod(west, 360.0)
    if w < 0.0:
        w += 360.0
    end = w + width
    if end <= 360.0:
        return [(w, end)]
    return [(w, 360.0), (0.0, end - 360.0)]


def _cells_overlap(a, b, radius, tol):
    """两格交叠面积;仅接边界为零面积,纬度带与经度段分别求交。"""
    south = max(a['south'], b['south'])
    north = min(a['north'], b['north'])
    if north - south <= tol:
        return 0.0
    height = math.sin(math.radians(north)) - math.sin(math.radians(south))
    width = 0.0
    for a0, a1 in a['lon']:
        for b0, b1 in b['lon']:
            lo = max(a0, b0)
            hi = min(a1, b1)
            if hi - lo > tol:
                width += hi - lo
    return radius * radius * math.radians(width) * height


def _read_cells(rows, radius):
    cells = {}
    for row in rows:
        cells[row['id']] = {
            'south': number(row['south']),
            'north': number(row['north']),
            'area': cell_area(row, radius),
            'lon': _lon_intervals(number(row['west']), positive(row['width'])),
        }
    return cells


def _reject_self_overlap(cells, radius, tol, label):
    """同一套网格内部重叠必须拒绝,不靠归一化抹掉双计。"""
    ids = list(cells)
    limit = tol * radius * radius
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if _cells_overlap(cells[a], cells[b], radius, tol) > limit:
                raise ValueError(f'{label}网格自身重叠:{a}与{b}')


def build_geometry(v):
    """构建稀疏交叠结构:只保存正面积的(源,目标)对及其面积。"""
    radius = positive(v['radius'])
    tol = positive(v.get('geometry_tolerance', DEFAULT_TOLERANCE))
    unique(v['source'], 'id')
    unique(v['target'], 'id')
    sources = _read_cells(v['source'], radius)
    targets = _read_cells(v['target'], radius)
    _reject_self_overlap(sources, radius, tol, '源')
    _reject_self_overlap(targets, radius, tol, '目标')
    limit = tol * radius * radius
    overlaps = []
    by_source = {sid: [] for sid in sources}
    by_target = {tid: [] for tid in targets}
    for sid, scell in sources.items():
        for tid, tcell in targets.items():
            area = _cells_overlap(scell, tcell, radius, tol)
            if area > limit:
                entry = {'source': sid, 'target': tid, 'area': area}
                overlaps.append(entry)
                by_source[sid].append(entry)
                by_target[tid].append(entry)
    return {
        'radius': radius,
        'tolerance': tol,
        'sources': sources,
        'targets': targets,
        'overlaps': overlaps,
        'by_source': by_source,
        'by_target': by_target,
        'covered_source': {sid: sum(e['area'] for e in by_source[sid]) for sid in sources},
        'covered_target': {tid: sum(e['area'] for e in by_target[tid]) for tid in targets},
    }


def _covariance_pairs(v, target_ids):
    pairs = []
    seen = set()
    for item in v.get('covariance_pairs', []):
        if not isinstance(item, list) or len(item) != 2:
            raise ValueError('协方差对须为两个目标id')
        a, b = item
        if a not in target_ids or b not in target_ids:
            raise ValueError('协方差对引用未知目标')
        if a == b:
            raise ValueError('协方差对须为不同目标')
        key = frozenset((a, b))
        if key in seen:
            continue
        seen.add(key)
        pairs.append((a, b))
    return pairs


def _frame_context(f, geom):
    valid = {sid: f['valid'][sid] for sid in geom['sources']}
    values = {sid: f['values'][sid] for sid in geom['sources']}
    variance = {sid: f['variance'][sid] for sid in geom['sources']}
    return valid, values, variance


def _weighted_variance(terms):
    """独立源方差按权重平方和传播;任一贡献源方差未知则不可估,不补零。"""
    total = 0.0
    for weight, var in terms:
        if var is None:
            return None, 'unknown_source_variance'
        total += weight * weight * var
    return total, 'ok'


def _intensive_covariance(a, b, out, geom, valid, variance):
    if out[a]['mean'] is None or out[b]['mean'] is None:
        no_cov = out[a]['valid_coverage_area'] == 0.0 or out[b]['valid_coverage_area'] == 0.0
        return {'pair': [a, b], 'value': None,
                'status': 'no_valid_coverage' if no_cov else 'below_min_coverage'}
    other = {e['source']: e['area'] for e in geom['by_target'][b] if valid[e['source']]}
    va = out[a]['valid_coverage_area']
    vb = out[b]['valid_coverage_area']
    value = 0.0
    for e in geom['by_target'][a]:
        sid = e['source']
        if not valid[sid] or sid not in other:
            continue
        var = variance[sid]
        if var is None:
            return {'pair': [a, b], 'value': None, 'status': 'unknown_source_variance'}
        value += (e['area'] / va) * (other[sid] / vb) * var
    return {'pair': [a, b], 'value': value, 'status': 'ok'}


def _intensive_frame(f, geom, min_coverage, pairs):
    valid, values, variance = _frame_context(f, geom)
    out = {}
    for tid, tcell in geom['targets'].items():
        contribs = [e for e in geom['by_target'][tid] if valid[e['source']]]
        valid_cov = sum(e['area'] for e in contribs)
        integral = sum(e['area'] * values[e['source']] for e in contribs)
        ratio = valid_cov / tcell['area']
        entry = {
            'full_area': tcell['area'],
            'geometric_coverage_area': geom['covered_target'][tid],
            'valid_coverage_area': valid_cov,
            'coverage_ratio': ratio,
            'coverage_integral': integral,
            'contributions': [
                {'source': e['source'], 'overlap_area': e['area'], 'weight': e['area'] / valid_cov}
                for e in contribs
            ],
        }
        if valid_cov > 0.0 and ratio + geom['tolerance'] >= min_coverage:
            entry['mean'] = integral / valid_cov
            entry['variance'], entry['variance_status'] = _weighted_variance(
                [(e['area'] / valid_cov, variance[e['source']]) for e in contribs])
        else:
            entry['mean'] = None
            entry['variance'] = None
            entry['variance_status'] = 'no_valid_coverage' if valid_cov == 0.0 else 'below_min_coverage'
        out[tid] = entry
    covariances = [_intensive_covariance(a, b, out, geom, valid, variance) for a, b in pairs]
    deposited = sum(t['coverage_integral'] for t in out.values())
    expected = sum(values[sid] * geom['covered_source'][sid] for sid in geom['sources'] if valid[sid])
    return {
        'id': f['id'],
        'targets': out,
        'covariances': covariances,
        'conservation': {
            'domain': 'overlap',
            'deposited_integral': deposited,
            'source_integral_on_overlap': expected,
            'residual': deposited - expected,
        },
    }


def _extensive_covariance(a, b, geom, valid, variance):
    sources = geom['sources']
    other = {e['source']: e['area'] for e in geom['by_target'][b] if valid[e['source']]}
    value = 0.0
    for e in geom['by_target'][a]:
        sid = e['source']
        if not valid[sid] or sid not in other:
            continue
        var = variance[sid]
        if var is None:
            return {'pair': [a, b], 'value': None, 'status': 'unknown_source_variance'}
        value += (e['area'] / sources[sid]['area']) * (other[sid] / sources[sid]['area']) * var
    return {'pair': [a, b], 'value': value, 'status': 'ok'}


def _extensive_frame(f, geom, pairs):
    sources = geom['sources']
    valid, values, variance = _frame_context(f, geom)
    out = {}
    for tid, tcell in geom['targets'].items():
        contribs = [e for e in geom['by_target'][tid] if valid[e['source']]]
        weights = {e['source']: e['area'] / sources[e['source']]['area'] for e in contribs}
        var, status = _weighted_variance([(weights[e['source']], variance[e['source']]) for e in contribs])
        out[tid] = {
            'full_area': tcell['area'],
            'geometric_coverage_area': geom['covered_target'][tid],
            'valid_coverage_area': sum(e['area'] for e in contribs),
            'allocated_total': sum(weights[e['source']] * values[e['source']] for e in contribs),
            'variance': var,
            'variance_status': status,
            'contributions': [
                {'source': e['source'], 'overlap_area': e['area'], 'weight': weights[e['source']]}
                for e in contribs
            ],
        }
    unallocated = {}
    for sid, scell in sources.items():
        if not valid[sid]:
            continue
        remainder = max(0.0, scell['area'] - geom['covered_source'][sid])
        unallocated[sid] = values[sid] * remainder / scell['area']
    covariances = [_extensive_covariance(a, b, geom, valid, variance) for a, b in pairs]
    allocated = sum(t['allocated_total'] for t in out.values())
    valid_total = sum(values[sid] for sid in sources if valid[sid])
    unallocated_total = sum(unallocated.values())
    return {
        'id': f['id'],
        'targets': out,
        'unallocated_by_source': unallocated,
        'covariances': covariances,
        'conservation': {
            'domain': 'overlap',
            'valid_source_total': valid_total,
            'allocated_total': allocated,
            'unallocated_total': unallocated_total,
            'residual': allocated + unallocated_total - valid_total,
        },
    }


def remap(v):
    """执行球面守恒重映射,返回可JSON序列化的结果。"""
    quantity = v.get('quantity')
    if quantity not in QUANTITIES:
        raise ValueError('quantity须为intensive或extensive')
    min_coverage = number(v['min_coverage'])
    if not 0.0 <= min_coverage <= 1.0:
        raise ValueError('min_coverage须在[0,1]内')
    geom = build_geometry(v)
    check_frames(v['frames'], set(geom['sources']))
    pairs = _covariance_pairs(v, set(geom['targets']))
    frames = []
    for f in v['frames']:
        if quantity == 'intensive':
            frames.append(_intensive_frame(f, geom, min_coverage, pairs))
        else:
            frames.append(_extensive_frame(f, geom, pairs))
    return {
        'quantity': quantity,
        'min_coverage': min_coverage,
        'geometry': {
            'radius': geom['radius'],
            'source_areas': {sid: c['area'] for sid, c in geom['sources'].items()},
            'target_areas': {tid: c['area'] for tid, c in geom['targets'].items()},
            'overlaps': [{'source': e['source'], 'target': e['target'], 'area': e['area']}
                         for e in geom['overlaps']],
            'source_covered_area': dict(geom['covered_source']),
            'target_covered_area': dict(geom['covered_target']),
        },
        'frames': frames,
    }
