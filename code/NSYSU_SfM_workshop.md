NSYSU_SfM_workshop
================
Yan
2025-05-03

# 簡介 / Introduction

本教學示範如何使用 `habtools` 套件計算三維網格 (3D mesh) 的複雜度指標。
This tutorial demonstrates how to compute complexity metrics for 3D
meshes using the `habtools` package.

## 前置作業與資料載入 / Setup and Data Loading

**中文說明**：載入必要的套件並從 `data/Goose_S1_cleaned.ply`
範例檔案載入網格資料。 **English Explanation**: Load the required
packages and import the example mesh from `data/Goose_S1_cleaned.ply`.

``` r
# 載入套件 / Load packages
library(habtools)
library(rgl)
library(Rvcg)

# 載入範例網格資料 / Load example mesh
# 請確保檔案在專案的 data 資料夾中
mesh <- vcgPlyRead("data/Goose_S1_cleaned.ply")
```

    ## Removed 2 duplicate 3168 unreferenced vertices and 6 duplicate faces

## 檢查網格 / Checking the Mesh

在計算任何指標前，先視覺化網格並檢查 z 軸方向是否正確。 Before
calculating any metrics, visualize the mesh and ensure the z orientation
is correct.

``` r
# 顯示三維網格視覺化 / Visualize the 3D mesh
plot3d(mesh)
```

## 解析度分佈 / Resolution Distribution

**為何需要此步驟？**
在處理三維網格資料時，各頂點之間的間距（解析度）往往不一致。如果不先確認和調整解析度分佈，後續的指標計算（例如分形維度和
Rugosity）可能會因解析度差異而產生偏差。此外，解析度分佈也有助於選擇適當的
voxelSize 參數以進行均勻重網格。 **Why this step?** Meshes often have
variable vertex spacing. Inspecting the distribution of distances
between adjacent vertices ensures that subsequent metric calculations
(e.g., fractal dimension and rugosity) are not biased by uneven
resolution. It also informs an appropriate choice of `voxelSize` for
uniform remeshing.

**中文說明**：計算網格中相鄰頂點間距並繪製直方圖以檢查解析度變異。
**English Explanation**: Compute distances between adjacent vertices in
the mesh and plot a histogram to inspect variation in resolution.

**中文說明**：計算網格中相鄰頂點間距並繪製直方圖以檢查解析度變異。
**English Explanation**: Compute distances between adjacent vertices in
the mesh and plot a histogram to inspect variation in resolution.

``` r
# 計算解析度向量 / Calculate vector of resolutions
resvec <- vcgMeshres(mesh)[[2]]

# 繪製直方圖 / Plot histogram of resolutions
hist(resvec, main = "Resolution Distribution", xlab = "Distance between vertices")
```

![](NSYSU_SfM_workshop_files/figure-gfm/resolution-histogram-1.png)<!-- -->

``` r
# 顯示解析度摘要統計 / Display summary statistics of resolutions
summary(resvec)
```

    ##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    ## 1.857e-06 1.258e-03 1.505e-03 1.703e-03 1.858e-03 4.430e-02

\$1 **何時需要均勻重網格？**
當解析度分佈非常寬（例如直方圖有長尾或多峰），表示網格上有極大或極小的頂點間距，這會影響後續的複雜度指標計算。此時建議進行均勻重網格，以獲得一致的解析度。相反地，如果直方圖大部分值集中在單一窄範圍，分佈均勻，則無需重新網格。
**When to remesh?** If the resolution histogram is very broad (e.g. long
tail or multimodal), indicating some vertex distances are much larger or
smaller than others, remeshing is recommended to achieve uniform
resolution. Conversely, if most values cluster tightly in a single
narrow range, the mesh is already uniform and no remeshing is needed.

``` r
# 均勻重網格 / Uniformly remesh mesh
mcap_uniform <- Rvcg::vcgUniformRemesh(mesh, silent = TRUE, multiSample = TRUE, voxelSize = min(resvec), mergeClost = TRUE)
```

    ## Error: vector::_M_default_append

``` r
Rvcg::vcgMeshres(mcap_uniform)[[1]]
```

    ## Error: object 'mcap_uniform' not found

``` r
# 驗證重網格後解析度 / Check resolution of the remeshed object
new_res <- vcgMeshres(mcap_uniform)[[2]]
```

    ## Error: object 'mcap_uniform' not found

``` r
hist(new_res, main = "Resolution Distribution", xlab = "Distance between vertices")
```

    ## Error: object 'new_res' not found

``` r
summary(new_res)
```

    ## Error: object 'new_res' not found
