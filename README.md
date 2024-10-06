# Seu_eda_pmm

东南大学自动化学院eda精英挑战赛第七题，矢量拟合算法实现

See teamwork/test_demo.m for getting start.

<!-- PROJECT SHIELDS -->


## 目录

- [项目实现](#项目实现)
  - [优化部分](#优化部分)
  - [待实现](#待实现)
  - [未实现](#未实现)
- [作者](#作者)
- [版权说明](#版权说明)



### 项目实现

基础实现：带权重的松弛vf算法拟合S参数，输出一个mat文件

###### 优化部分

1. 字符索引加速读取文件
2. 谐振波判断初始阶数（factor因子需要测试）
3. 并行处理VF的QR分解
4. 极点残差融合
5. VF转PA预测最佳阶数
6. half_size S参数的被动性矩阵P
7. P端口网络的列空间压缩

###### 待实现
none

###### 未实现
1. 无源性和dc的同时实现
2. 大规模电路dc_preserve的优化（现有的二次规划构建方程过大，限制能使用的最大阶数）


### 作者

蔡雨洋：[email-redacted]

知乎:我在 ;     

 *属于刘知秋团队，同组成员除作者还有李奥*

### 版权说明

代码基于https://github.com/yezuochang/pmm 和Matrix Fitting 工具箱实现