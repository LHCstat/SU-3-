# 三阶幺正矩阵探索器 | U(3) / SU(3)

从九个反厄米生成元出发，利用矩阵指数的级数重组构造任意三阶幺正矩阵。

## 交互式可视化

下载仓库（Code → Download ZIP），解压后用浏览器打开 `index.html`。无需服务器或第三方 JavaScript 依赖。

- 九个参数滑块，实时生成复数矩阵。
- U(3) / SU(3) 切换，后者自动约束 α₃ = −α₁ − α₂。
- 矩阵元素模方热图、九个元素在复平面上的向量图。
- 显示 U†U、行列式及幺正性误差。
- 随机参数、单位矩阵重置和 JSON 复制。

## Python

```bash
pip install numpy
python unitary.py --test
python unitary.py --su3
python unitary.py --degrees --angles 30 20 45 0 60 0 0 0 0
```

```python
from unitary import unitary3
U = unitary3(0.3, 0.2, 0.7, phi13=1.0, special=True)
```

## 数学原理

对于 i < j，设 Eᵢⱼ 是只有 (i,j) 元素为 1 的矩阵：

- Aᵢⱼ = Eᵢⱼ − Eⱼᵢ
- Bᵢⱼ = i(Eᵢⱼ + Eⱼᵢ)
- Cⱼ = iEⱼⱼ

令 Gᵢⱼ = Aᵢⱼ cos φ − Bᵢⱼ sin φ，则 G³ = −G，因而

`U_ij = I + G_ij sin θ + G_ij² (1 − cos θ)`

最终 U = D U₂₃ U₁₃ U₁₂，其中 D = diag(exp(iα₁), exp(iα₂), exp(iα₃))。

该分解通过复数 Givens 消元覆盖 U(3)。SU(3) 要求三个对角相位之和为 0 (mod 2π)。当所有相位为零时退化为 SO(3)。

随机生成的是均匀选取参数的矩阵，**并非 Haar 均匀分布**。一般参数化也并非全局一一对应。
