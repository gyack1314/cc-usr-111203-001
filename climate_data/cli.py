"""显示已有原始摘要或执行球面守恒重映射,来源文件保持不动。"""
import argparse
import json
from .io import load
from .domain import summarize
from .remap import remap

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['inspect','remap']);parser.add_argument('file');a=parser.parse_args()
    v,digest=load(a.file)
    if a.command=='remap':
        print(json.dumps({'source_sha256':digest,'remap':remap(v)},ensure_ascii=False,indent=2,allow_nan=False))
    else:
        print(json.dumps({'source_sha256':digest,'summary':summarize(v)},ensure_ascii=False,indent=2,allow_nan=False))
if __name__=='__main__':main()
