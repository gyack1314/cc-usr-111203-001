# 海洋资料换网格的球面守恒重映射

这是独立的气候资料Python项目。`data/`保存合成资料，`climate_data.io`负责严格读取、公历日期和来源字节摘要，`climate_data.domain`提供现有描述统计，`climate_data.remap`执行球面守恒重映射（稀疏交叠权重、格内平均与总量分配、方差与目标对协方差传播、交叠域守恒核对），`climate_data.cli`将摘要与重映射结果输出为JSON。所有观测值、站点、人群、情景与金额都是虚构。字段和研究约定见[资料约定](docs/data-contract.md)。

当前Python 3.11及以上可直接使用，没有第三方Python依赖。项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q climate_data tests
python -m climate_data.cli inspect data/normal.json
python -m climate_data.cli remap data/normal.json
```

`inspect`显示当前原始摘要，不覆盖资料。`remap`输出稀疏交叠权重、逐帧目标估计（平均或总量分配）、逐项来源贡献、目标对协方差与守恒核对；几何跨帧复用，有效掩码逐帧重算分母，来源文件保持不动。`invalid-`前缀文件为不合法资料，其他变体使用相同字段保留业务边界。来源字节SHA-256随摘要输出，后续结果应同时标明输入身份、量纲和适用范围。

容器使用固定基础镜像，可运行`docker build -t climate-material .`，随后`docker run --rm climate-material`检查资料。`docker run --rm climate-material python -m unittest discover -s tests -v`执行已有检查。运行不依赖远端账号或实时服务。
