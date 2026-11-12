# 球面网格口径

west可为任意有限经度，width为向东角宽，不由east-west推断；width=360表示完整纬带。矩形纬边是等纬线，经边是经线，边界接触不算面积覆盖。radius的平方决定面积单位，演示radius=1取单位球面，实际使用须与物理量单位一致。面积本身是正数，不把极点处零长度边当成零面积格。

frames按id独立，values/valid/variance齐全引用源格，null方差表示未知；独立方差假设只跨源格，不跨目标格。quantity=intensive表示格内平均量，extensive表示源格总量。目标平均的权重为交叠面积除当前有效覆盖面积，达到min_coverage才有正式值；未达到时覆盖积分仍可展示但应明确不可作为完整目标估计。未知方差只污染使用该源格的目标及目标对。

现有工具给每格原始面积，不求交叠、不推送数值或协方差。

# 重映射输出口径

`remap`在`inspect`之外给出球面守恒重映射，来源文件保持不动。`geometry`只列正面积交叠（稀疏结构，不分配完整源乘目标矩阵），另给源/目标全面积与几何覆盖面积；同一套网格内部重叠直接拒绝，不做归一化。每帧独立计算：几何交叠跨帧复用，有效掩码与分母按当前帧重算，不沿用上一时刻。

intensive目标每项列`full_area`、`geometric_coverage_area`、`valid_coverage_area`、`coverage_ratio`、`coverage_integral`与`contributions`（逐源交叠面积及权重，权重=交叠/当前有效覆盖）。`coverage_ratio`达到`min_coverage`才返回`mean`，否则`mean`为null；覆盖积分始终展示但不代表完整目标估计。extensive目标列`allocated_total`=Σ有效源总量×交叠/源全面积，不再归一化；`unallocated_by_source`列出没有去向的份额，无效源总量不参与分配。

方差按权重平方和传播；来源方差未知（null）时相应不确定性不可估，`variance`为null且`variance_status=unknown_source_variance`，不补零。`covariances`按`covariance_pairs`输出：共享同一源格的目标不独立，协方差含非对角项；共享源方差未知则该对不可估；intensive的协方差要求两侧`mean`均达到阈值。`conservation`只在交叠覆盖域核对：intensive比较`deposited_integral`与`source_integral_on_overlap`，extensive核对`allocated_total`+`unallocated_total`=`valid_source_total`；没有覆盖的目标不参与守恒要求。`variance_status`取值：`ok`、`unknown_source_variance`、`below_min_coverage`、`no_valid_coverage`。
