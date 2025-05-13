# SfM 珊瑚生態調查工作坊 / SfM Coral Survey Workshop

## 專案介紹 / Project Overview

本專案為一場專注於運動恢復結構 (SfM) 3D 調查技術的教學與實作工作坊，涵蓋拍攝方法、Metashape 軟體操作、3D 模型生成與後續分析步驟。
This project is a hands-on workshop focused on Structure-from-Motion (SfM) 3D survey techniques, covering filming methods, Metashape software operations, 3D model generation, and subsequent analysis steps.

---

## 行程表 / Schedule

| 時間 / Time    | 行程 / Activity                                                                                          |
| ------------ | ------------------------------------------------------------------------------------------------------ |
| **Day 1**    |                                                                                                        |
| 2:00–2:20 PM | 開場：為何需要 3D 調查、SfM 原理簡介 / Introduction: Why 3D Survey & SfM Basics                 |
| 2:20–2:50 PM | 群體尺度拍攝方式說明與調查練習 / Colony-Scale Survey Tutorial                       |
| 3:00–3:30 PM | 棲地尺度拍攝方式說明與調查練習 / Habitat-Scale Survey Tutorial                         |
| 3:30–4:30 PM | [Metashape 軟體操作](docs/NSYSU-MetashapeSOP-2505.pdf) / Metashape Software Demo                            |
| 4:30–4:45 PM | 簡易 3D 模型動畫與輸出展示 / Simple 3D Model Export & Visualization                     |
| 4:45–4:55 PM | 小結與方法設備差異討論、相關資源分享 / Summary & Discussion                                                              |
| 4:55–5:00 PM | Q\&A                                                                                                   |
| **Day 2**    |                                                                                                        |
| 2:00–2:20 PM | 如何分析 3D 調查結果、實際應用案例分享 / Analyzing 3D Survey Data & Case Studies       |
| 2:20–2:50 PM | 地理參照與 DEM／Orthomosaic 生成 / Georeferencing & DEM/Orthomosaic Generation |
| 3:00–3:10 PM | Metashape 完善群體尺度模型 / Improving Colony-Scale Model                   |
| 3:10–3:50 PM | [R 語言群體尺度複雜度分析](code/NSYSU_SfM_workshop.md) / R-Based Colony-Scale Analysis                     |
| 4:00–4:40 PM | [R + QGIS 棲地尺度複雜度分析](code/NSYSU_dem_analysis.md) / R + QGIS Habitat-Scale Analysis                 |
| 4:40–4:55 PM | 進階分析方法分享 / Advanced Analysis Methods                                                                   |
| 4:55–5:00 PM | Q\&A                                                                                                   |

---

## 要帶的事項 / To Bring

* 個人電腦（附充電線）/ Laptop (with charger)
* 相機或能拍照的手機 / Camera or smartphone capable of taking photos
* 讀卡機或傳輸線 / Card reader or USB cable to transfer photos to computer
* 想要掃描成 3D 模型的最愛小物 / Favorite object to scan into a 3D model

---

## 需要先下載的軟體 / Required Software

* [Metashape Professional 2.1.4](https://www.agisoft.com/downloads/installer/)
* [R & RStudio](https://posit.co/download/rstudio-desktop/)

  * Windows 使用者請同時 [下載並安裝 Rtools](https://cran.r-project.org/bin/windows/Rtools/)
* 安裝以下 R 套件 / Install the following R packages:

  ```r
  install.packages(c("habtools", "rgl", "Rvcg", "raster", "sf", "dplyr", "ggplot2"))
  ```
* [QGIS](https://qgis.org/en/site/forusers/download.html)
* [Sketchfab 帳號註冊 / Sign up for Sketchfab](https://sketchfab.com/signup)

### 非必要但值得下載嘗試 / Optional but Recommended Software

* [Meshlab](https://www.meshlab.net/) — 一款開源 3D 網格處理工具，可用於檢視、編輯與優化 3D 網格模型 / An open-source 3D mesh processing tool for viewing, editing, and optimizing 3D mesh models.

* [CloudCompare](https://www.danielgm.net/cc/) — 一個點雲與網格比較與處理軟體，可用於點雲配準、濾波與分析 / A point cloud and mesh processing software for registration, filtering, and analysis.

---

## 課程簡報 / Course slides

* 第一天 / [Day 1](docs/DAY1-3Dreef_survey_intro.pdf) - Intro to 3D Coral Reef Survey

* 第二天 / Day 2 - 3D Analysis and Case Studies

---

## 簡易 R 操作指南 / Quick R Usage Guide

1. 開啟 RStudio / Open RStudio
2. 在 Console 輸入 `install.packages(...)` 安裝套件
3. 使用 `library(habtools)` 等命令載入套件 / Load packages with `library()`

### 相關學習資源 / Learning Resources

* [R Cheat Sheet](https://iqss.github.io/dss-workshops/R/Rintro/base-r-cheat-sheet.pdf) — 一頁式 R 基本語法速查表 / One-page R syntax cheat sheet
* [R Tutorial for Beginners](https://www.w3schools.com/r/default.asp) — 入門線上教學資源 / Online beginner tutorial

---

## 簡易 QGIS 操作指南 / Quick QGIS Usage Guide

[簡報slide](docs/NSYSU-QGIS-guide.pdf)

1. 載入 Orthomosaic 與 DEM：

   * 選擇「圖層」→「新增圖層」→「新增光柵圖層」，載入正射影像與 DEM 檔案
     Add Orthomosaic & DEM:
   * Layer → Add Layer → Add Raster Layer, then select your orthomosaic and DEM files.

2. 創建 Polygon 進行珊瑚群體圈選：

   * 使用「編輯圖層」工具啟用編輯，再選取「新增要素」，繪製多邊形圈選珊瑚群體
     Create Polygons for Coral Selection:
   * Toggle editing with the Toggle Editing tool, then use the Add Feature tool to draw polygons around coral colonies.

3. 計算面積與邊長：

   * 在圖層面板右鍵點選多邊形圖層 → 開啟屬性表 → 點擊「欄位計算器」
   * 新增字段並使用運算符 `"$area"` 與 `"$perimeter"` 計算面積與周長
     Calculate Area & Perimeter:
   * Right-click your polygon layer in Layers panel → Open Attribute Table → Toggle Field Calculator
   * Create new fields using expressions `"$area"` and `"$perimeter"` to compute metrics.

---

## 給初次github使用者的指引 / Download & Setup guide for github beginner

* 在 GitHub 專案頁面點擊 **Code** → **Download ZIP**
  On the GitHub repo page, click **Code** → **Download ZIP**

* 將下載的 ZIP 檔解壓縮到您想存放的資料夾
  Unzip the downloaded archive to your desired folder

* 開啟 RStudio，選擇 **File** → **Open Project** 或 **Open Directory**，導覽至剛剛解壓縮的專案資料夾並打開
  In RStudio, use **File** → **Open Project** or **Open Directory**, then navigate to and open the unzipped project folder

* 設定工作目錄：在 RStudio Console 中輸入 `setwd("/path/to/your/project")`，或使用 RStudio 上方選單 **Session** → **Set Working Directory** → **To Project Directory**
  Set the working directory by running `setwd("/path/to/your/project")` in the Console, or via **Session** → **Set Working Directory** → **To Project Directory**

* Day 2 的 R 分析前，請先在 `code/` 資料夾中開啟並測試以下檔案是否可以正常執行：
  Before Day 2 R analysis, please open and test the following files in the `code/` directory:

  * `NSYSU_SfM_workshop.Rmd`
  * `NSYSU_dem_analysis.Rmd`

* 如有軟體安裝或環境設定問題，我將於工作坊開始前 30 分鐘抵達現場協助解決
  If you encounter any software installation or environment setup issues, I will arrive 30 minutes before the workshop to assist.



---


*README.md generated for 2025 NSYSU\_SfM\_workshop. By Guan-Yan Chen 陳冠言.*
