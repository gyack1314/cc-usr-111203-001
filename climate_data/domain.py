"""核对单元字段并计算每格的原始球面面积。"""
import math
from .io import unique,number,positive

def cell_area(cell,radius):
    west=number(cell['west']);width=positive(cell['width']);south=number(cell['south']);north=number(cell['north'])
    if width>360 or not -90<=south<north<=90:raise ValueError('格网边界超出范围')
    return radius*radius*math.radians(width)*(math.sin(math.radians(north))-math.sin(math.radians(south)))

def summarize(v):
    radius=positive(v['radius']);source=unique(v['source'],'id');unique(v['target'],'id')
    areas={r['id']:cell_area(r,radius) for r in v['source']};target_areas={r['id']:cell_area(r,radius) for r in v['target']}
    unique(v['frames'],'id');frames=[]
    for f in v['frames']:
        if set(f['values'])!=source or set(f['valid'])!=source or set(f['variance'])!=source:raise ValueError('每帧格引用需齐全')
        valid=0
        for sid in source:
            if type(f['valid'][sid]) is not bool:raise ValueError('有效性必须布尔')
            if f['values'][sid] is not None:number(f['values'][sid])
            if f['valid'][sid] and f['values'][sid] is None:raise ValueError('有效格缺少数值')
            if f['variance'][sid] is not None and number(f['variance'][sid])<0:raise ValueError('方差不能负')
            valid+=f['valid'][sid]
        frames.append({'id':f['id'],'valid_source_cells':valid})
    return {'source_areas':dict(sorted(areas.items())),'target_areas':dict(sorted(target_areas.items())),'frames':frames}
