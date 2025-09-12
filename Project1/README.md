### 四個結果圖總結

| 區塊位置 | 模型                     | 圖片內容    | 設定條件                            |
| ---- | ---------------------- | ------- | ------------------------------- |
| 左上   | Gaussian Rayleigh      | PDF 直方圖 | T=1、N=8、v=60km/h                |
| 右上   | Rayleigh fading（含 LOS） | 時域變化圖   | T=1、N=8、v=60km/h、k=2 (K-factor) |
| 左下   | Jakes Model            | PDF 直方圖 | 同上                              |
| 右下   | Jakes Model            | 時域變化圖   | 同上                              |

---

### ✅ 差異分析

#### 1. 左上 Gaussian Rayleigh PDF vs 左下 Jakes PDF

* **左上（Gaussian Rayleigh）** 符合理論 Rayleigh 分布，滑順單峰形狀。
* **左下（Jakes 模擬）** 會看到抖動、多峰、不連續，這是因為：

  * 散射角數 N=8 太小，不夠逼近中心極限定理。
  * 相位未完全隨機，可能造成 constructive/destructive interference。

✅ **解法**：提高 N 至 30 或 100，會更像 Rayleigh 理論曲線。

#### 2. 右上 vs 右下 的時域波形比較

* **右上** 明顯是 Rician Fading 的樣子（因為有 K=2）

  * 包絡波動「不會到 0」→ 有一條穩定直視徑（LOS）在。
  * 因此 dips 較淺，顯示的是 **Rician fading**。

* **右下** 則是 Jakes Rayleigh 模型：

  * 波動範圍更大（有 deep fading），代表沒有直視徑，只有多徑。

---

## 📡 Rayleigh Fading Channel 模擬與分析筆記

### 🎯 模擬目的

* 比較 **Jakes Model** 與 **Gaussian Rayleigh** 的 fading 特性
* 分析在改變不同參數（T、v、N）下，fading 行為與分布的變化

---

## 🔍 模擬模型說明

### ✅ Gaussian Rayleigh Model

* 實作：實部與虛部為獨立高斯分佈 N(0,σ²/2)
* 結果：|c| 服從 Rayleigh 分布
* 屬於「**統計建模**」，非物理建模

### ✅ Jakes Model

* 實作：模擬 N 條路徑，每條有不同入射角與 Doppler shift
* 相位角 \$,\phi\_k \sim \mathcal{U}(0,2\pi)\$，頻率偏移 \$\omega\_k = 2\pi f\_D \cos(\theta\_k)\$
* 屬於「**物理建模**」

---

## 🧪 實驗內容與觀察

### 1️⃣ 改變 **T = 1 → 1.5s**（延長觀察時間）

| 條件 | N=8，v=60km/hr |
| -- | ------------- |

#### 結果：

* 分布圖更平滑 → 更接近理論 Rayleigh 分布
* fading time waveform 變得更連續（取樣更密）
* Gaussian Rayleigh 衰減變平坦
* Jakes model 衰減更密集但平滑

👉 **延長 T 可減少統計誤差，更接近真實分布**

---

### 2️⃣ 改變速度 **v = 60 → 45 km/hr**

| 條件 | N=8，T=1s |
| -- | -------- |

#### 結果：

* 較慢速度 → Doppler shift 減少
* fading 的波動頻率變小
* deep fading 間隔變大

👉 **v 小 → fading 頻率變小 → 衰減間距變大**

---

### 3️⃣ 改變多徑數量 **N = 8 → 20 → 30**

| 條件 | T=1s，v=60km/hr |
| -- | -------------- |

#### 結果：

* Jakes model 在 N=20/30 時分布圖明顯更接近 Rayleigh
* time waveform 更穩定，波動更平滑
* Deep fading 減少

👉 **N ↑ → 相位總合越趨近高斯 → 中心極限定理效果 → 越像 Rayleigh**

---

### 4️⃣ 附加 fading waveform（k=2）

* 顯示不同組合下的 fading magnitude 隨時間變化
* 比較 Gaussian vs Jakes model

#### 結果：

* Gaussian: 看起來較隨機
* Jakes: 頻率固定感強，波峰密集

---

## 🧠 總結表格

| 參數變化    | 結果觀察                          | 結論            |
| ------- | ----------------------------- | ------------- |
| **T ↑** | 分布收斂穩定，波動平滑                   | 統計性質穩定，接近理論   |
| **v ↑** | Doppler 增，加快衰減頻率              | fade 間隔變短，較擁擠 |
| **N ↑** | Jakes 趨近 Rayleigh，waveform 平滑 | 相位混合多，類高斯     |

---

### ✅ 模擬驗證了：

* Jakes Model 可合理逼近 Rayleigh 分布（在 N 夠大時）
* Gaussian Rayleigh 是統計建模但能快速模擬
* 通道參數對 fading 行為有顯著影響

---

> 若需加入公式推導或圖表說明，也可額外補充進來。
