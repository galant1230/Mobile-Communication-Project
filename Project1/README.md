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
### 📊 模擬參數變化分析與觀察結果

---

#### 1. 改變 \$T = 1.5\$ （觀察總取樣時間延長的影響）

* 當 \$T\$ 增加，總時間延長，都卜勒偏移 \$\omega\_d = 2\pi f\_D\$ 並未改變，但因為時間軸變長，**每個樣本間的時間間距 \$\Delta t\$ 相對更密集**，會使得整體觀測時間的頻率分佈更完整。
* 對於 **Jakes 模型** 而言，這有助於彌補「散射角不夠多」造成的 PDF 抖動問題。

**觀察結果：**

* Gaussian Rayleigh 模型的 PDF 更平滑
* Jakes 模型 的分佈也變得更接近 Rayleigh 分佈
* 時域衰減圖中，fading 趨勢變得更密集、頻率更高

---

#### 2. 改變 \$v = 45\$ km/hr（觀察速度變化影響 Doppler Spread）

* 都卜勒偏移頻率 \$f\_D = \frac{v}{c}f\_c\$ 隨著 \$v\$ 降低而變小
* 表示 **channel fading 的頻率降低**，造成：

  * 訊號衰減之間的時間間隔變長（fading 間距變大）
  * 深度衰減（deep fading）出現頻率變少

**在 PDF 圖形上，接近理想 Rayleigh 分佈，原因如下：**

* 較小的 Doppler Spread 造成頻率干擾較少
* 訊號 envelope 穩定性較高

---

#### 3. 改變路徑數 \$N = 20\$、\$30\$（觀察中心極限定理效果）

* Jakes Model 中，每條路徑對應一個相位：

  $\alpha_m = \frac{2\pi m}{N},\quad m = 0, 1, ..., N-1$

* 當 \$N\$ 增加時：

  * 相位組合的數量更多、覆蓋更均勻
  * 有助於滿足 **中心極限定理（CLT）** → Re/Im 接近高斯分佈

**觀察結果：**

* Jakes 模型的分佈曲線越來越像 Gaussian Rayleigh
* 通道衰減也更穩定，出現極端 fading 的機率降低

---

### 🔍 總結整理表格

| 參數變化     | 對 Rayleigh PDF 的影響        | 對 Fading Waveform 的影響 |
| -------- | ------------------------- | --------------------- |
| 增加 \$T\$ | 資料點數增多 → 更平滑              | 更高時間解析度，變化頻繁          |
| 降低 \$v\$ | Doppler spread 變小 → 更穩定   | Fading 間距變大、變緩        |
| 增加 \$N\$ | 滿足 CLT → 趨近理想 Rayleigh 分佈 | 波形更穩定、深衰減減少           |

---

### 📡 Doppler Shift 與模型關聯

在 Jakes 模型中，每條路徑的相位變化如下：

$$
h(t) = \sum_{m=1}^{N} a_m e^{j(\omega_D \cos(\theta_m) t + \phi_m)}
$$

其中：

* \$\omega\_D = \frac{2\pi v}{\lambda}\$ 是最大都卜勒角頻率
* \$\cos(\theta\_m)\$ 為路徑入射角對移動方向的投影
* \$\phi\_m\$ 為隨機相位

對應到 MATLAB 程式碼：

```matlab
h(i,:) = cn * exp(j * (cos(rfa_m(i)) * w_m * n' + phi_m));
```

| 程式碼變數           | 對應公式符號              | 意義說明         |
| --------------- | ------------------- | ------------ |
| `w_m`           | \$\omega\_D\$       | 最大都卜勒頻率      |
| `cos(rfa_m(i))` | \$\cos(\theta\_m)\$ | 入射角與移動方向夾角餘弦 |
| `n'`            | \$t\$               | 時間軸取樣點       |
| `phi_m`         | \$\phi\_m\$         | 多徑隨機相位       |

### 🧠 實體意涵：

* **移動越快（v ↑）→ 都卜勒偏移越大 → 相位變化速度越快 → fading 更快更密**
* **入射角方向越對準移動方向（cos θ 趨近 1）→ 造成最大 Doppler 效果**

---

