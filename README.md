# 曝光融合（C++17 实现）

实现 **Mertens、Kautz、Van Reeth，《Exposure Fusion》，Pacific Graphics 2007** 的曝光融合算法。
纯 C++17，仅用 [stb](https://github.com/nothings/stb) 读写 PNG，不依赖 OpenCV。

## 编译

```bash
cd project
make                 # 生成 build/exposure_fusion
```

## 运行

```bash
./build/exposure_fusion \
  --inputs data/inputs/gan/o1.JPG,data/inputs/gan/u1.JPG \
  --outdir data/output/kitchen \
  --levels 10
```

参数：

| 参数 | 说明 |
| --- | --- |
| `--inputs a,b[,c...]` | 同一场景、不同曝光的照片（2 张及以上） |
| `--outdir DIR` | 输出目录 |
| `--levels N` | 金字塔层数（会按图像尺寸自动截断） |
| `--sigma` | 曝光适中度的高斯 σ，默认 0.2 |
| `--weight-blur N` | 权重图平滑次数，默认 6（消除中间的光晕/阴影） |
| `--no-pyramids` | 不保存金字塔层级图 |

输出：

```
fused.png            多分辨率融合结果
naive.png            直接加权平均（用于对比接缝/光晕）
weights/weight_i.png 每张图的归一化权重
measures/            对比度 / 饱和度 / 曝光适中度
pyramid/             权重的高斯金字塔 + 结果图的拉普拉斯金字塔
stats.json           分辨率、层数、耗时
```

2000×1328 的两张图、10 层金字塔，单核约 1 秒。

## 算法

对第 `i` 张图：

1. **对比度** `C = |Laplacian(gray(I))|`
2. **饱和度** `S = std(R, G, B)`
3. **曝光适中度** `E = ∏_c exp(−(I_c − 0.5)² / 2σ²)`

组合并在所有图片间归一化：

```
W_i = C_i^wc · S_i^ws · E_i^we        W_i ← W_i / Σ_j W_j
```

直接加权平均 `F = Σ W_i · I_i` 会出现接缝和光晕。

多分辨率融合：

```
权重的高斯金字塔：      G(W_i)^l
图像的拉普拉斯金字塔：  L(I_i)^l
逐层融合：              L_F^l = Σ_i G(W_i)^l · L(I_i)^l
重建：                  F = reconstruct(L_F)
```

工程里为了让 8-bit 照片稳定，做了三处小处理：

- 对比度、饱和度加一个很小的下限（0.02），平坦区域交给曝光适中度决定，
  避免 JPEG 噪声主导权重；
- 权重图在归一化前做 6 次平滑（`--weight-blur`），避免像素级闪烁；
- 金字塔层数按图像尺寸自动截断，最粗一层不小于 1 像素。

## 输入素材

`data/inputs/gan/o1.JPG`（过曝）与 `u1.JPG`（欠曝）是 **MEF-GAN** 测试集里的真实照片，
用于本项目的讲解视频；`data/inputs/belgium` 保留了 **MEFB** 的 Belgium House 对，作为第二个例子。

- MEFB：X. Zhang，《Benchmarking and Comparing Multi-exposure Image Fusion Algorithms》，
  Information Fusion 2021 — https://github.com/xingchenzhang/MEFB
- MEF-GAN：H. Xu 等 — https://github.com/hanna-xu/MEF-GAN
- 原论文：T. Mertens, J. Kautz, F. Van Reeth，《Exposure Fusion》，PG 2007

## 视频素材

`scripts/make_video_assets.py` 会把 C++ 输出转成 base64 JPEG 纹理
（`src/exposure/assets.json`），供 Remotion + Three.js 讲解视频使用。需要 Pillow：

```bash
PYTHONPATH=.pylibs /home/lchin/miniconda3/bin/python scripts/make_video_assets.py
```
