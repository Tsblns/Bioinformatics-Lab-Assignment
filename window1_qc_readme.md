## window1_qc.py 目前实际做了什么

它输出的东西可以这样对应：

### 1）数据清点

对应文件：

- dataset_qc_summary.xlsx 里的 dataset_summary

- missingness/ 目录下各组学缺失统计

这部分回答的是：

- 每个数据集有多大

- 样本在行还是列

- 特征是什么

- 缺失大概怎样

------

### 2）样本对齐与过滤

对应文件：

- dataset_qc_summary.xlsx 里的 sample_filtering

这部分回答的是：

- 原始样本有多少

- 去掉非 -01 后还剩多少

- 去掉了多少正常/非原发样本

- 最终有多少肿瘤样本、多少病人

------

### 3）共同样本表

对应文件：

- master_samples_6omics.csv

- master_samples_4omics.csv

- sample_presence_6omics.csv

- clinical_survival_overlap.csv

- dataset_qc_summary.xlsx 里的 overlap_summary

这部分回答的是：

- 六组学共同样本有多少

- 四组学主分析共同样本有多少

- clinical / survival 跟这些共同样本能对上多少

------

### 4）后续预处理接口说明

对应文件：

- preprocess_rules.md

- window_input_spec.md

这部分其实已经在给窗口2以后铺路了：

- 哪些组学怎么预处理

- 默认读 master_samples_4omics.csv 还是 master_samples_6omics.csv

- 输出矩阵格式统一成什么样

------

## 所以你现在所在的位置

你现在已经完成了：

- 数据清点

- 样本过滤

- 共同样本对齐

- 主分析/扩展分析分流

- 后续窗口的输入规则定义

也就是说，窗口1已经完成。

------

## 下一步是什么

对，下一步就是窗口2之后的各组学预处理。

按你们现在这套方案，最合理的顺序是：

### 先并行做三条线

- 窗口2：RNA + miRNA + RPPA 预处理

- 窗口3：CNV + mutation 预处理

- 窗口4：methylation 预处理

### 然后再进入

- 窗口5：多组学整合与聚类

------

## 但这里要注意一个实际点

虽然 master_samples_6omics.csv 已经有了，主分析最好先走四组学：

- RNA

- miRNA

- CNV

- RPPA

因为你前面已经确认过：

- 甲基化当前是部分数据

- 突变当前也是部分数据

所以更稳的是：

- 窗口2先开工

- 窗口3做时区分 baseline 和扩展

- 窗口4作为扩展模块保留

- 窗口5优先先跑四组学主分析

------

## 一句话总结

是的，window1_qc.py 基本就是把

“数据清点 → 样本对齐 → 共同样本表 → 后续输入规范”

这一整段都做完了。

下一步就是进入窗口2/3/4的各组学预处理。