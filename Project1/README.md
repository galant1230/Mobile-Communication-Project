# 📡 1. Rayleigh Fading Channel 模擬與分析筆記

## 🎯 模擬目的

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
## 📊 模擬參數變化分析與觀察結果

### 1. 改變 \$T = 1.5\$ （觀察總取樣時間延長的影響）

* 當 \$T\$ 增加，總時間延長，都卜勒偏移 \$\omega\_d = 2\pi f\_D\$ 並未改變，但因為時間軸變長，**每個樣本間的時間間距 \$\Delta t\$ 相對更密集**，會使得整體觀測時間的頻率分佈更完整。
* 對於 **Jakes 模型** 而言，這有助於彌補「散射角不夠多」造成的 PDF 抖動問題。

**觀察結果：**

* Gaussian Rayleigh 模型的 PDF 更平滑
* Jakes 模型 的分佈也變得更接近 Rayleigh 分佈
* 時域衰減圖中，fading 趨勢變得更密集、頻率更高

---

### 2. 改變 \$v = 45\$ km/hr（觀察速度變化影響 Doppler Spread）

* 都卜勒偏移頻率 \$f\_D = \frac{v}{c}f\_c\$ 隨著 \$v\$ 降低而變小
* 表示 **channel fading 的頻率降低**，造成：

  * 訊號衰減之間的時間間隔變長（fading 間距變大）
  * 深度衰減（deep fading）出現頻率變少

**在 PDF 圖形上，接近理想 Rayleigh 分佈，原因如下：**

* 較小的 Doppler Spread 造成頻率干擾較少
* 訊號 envelope 穩定性較高

---

### 3. 改變路徑數 \$N = 20\$、\$30\$（觀察中心極限定理效果）

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
# 📈 2. Rician Channel 與參數變化分析

Rician fading impulse response:

$$
h(t) = Ae^{j\omega_n t} + E_0 C_m \sum_{m=1}^{N} e^{j(\omega_n t + \psi_m)}
$$

Rician factor：

$$
K = \frac{A^2}{\sigma^2} \Rightarrow A = \sqrt{K\sigma^2}
$$

---

### 各參數分析：

1. **振幅 \$|h(t)|\$**

   * K 越大 → LOS 成分 A 越大 → 振幅越大
   * 分佈圖向右偏移，模擬曲線更像 Gaussian（近似常數包絡）

2. **頻率響應**

   * K = 5 時頻率響應最平穩，接近 Gaussian
   * K = 10 時過度偏移，造成過大 response
   * K = 2 時雖然存在 Rayleigh-like 多變性，但偏移不足

3. **T（觀測時間長度）**

   * 只有在 K = 10 時影響明顯（過強 LOS 使響應過大）
   * 其他 K 值下，改變 T 無明顯變化，可能因 T 改變幅度小

4. **V（速度）**

   * 速度下降 → 都卜勒減少 → 對應曲線更接近理想分佈
   * 在 K = 2, 10 下，速度降低反而使得 response 偏移加劇

5. **K（Rician 因子）**

   * K 越大 → LOS 成分越強 → Rician 分佈中 Gaussian 成分主導
   * 當 \$A^2 \gg \text{Var}(X), \text{Var}(Y)\$，則包絡接近以 A 為均值的高斯分佈

6. **N（多徑路徑數）**

   * \$C\_m = \sqrt{1/N}\$
   * N 增加 → 多徑成分下降 → LOS A 影響更明顯
   * 最終使得分佈更接近以 A 為中心的 Gaussian 分佈

---

### 📘 小結：Jakes vs Rician 模型參數影響

| 模型     | 關鍵參數  | 效果說明                        |
| ------ | ----- | --------------------------- |
| Jakes  | \$N\$ | 多徑平均 → CLT 生效 → 趨近 Rayleigh |
| Jakes  | \$v\$ | 都卜勒擴展變大，小變化更頻繁              |
| Jakes  | \$T\$ | 時間解析度提升，PDF 更平滑             |
| Rician | \$K\$ | 越大 → LOS 強 → 分佈趨近 Gaussian  |
| Rician | \$N\$ | 多徑減弱 → 強調直射 LOS             |
| Rician | \$v\$ | 影響 Doppler spread，頻率響應寬度變化  |

---
# 📘 3. BER vs. Eb/No 模擬比較筆記

## 🔬 實驗目的 (b)

比較兩種通道模型對 BPSK 傳輸系統的 BER（位元錯誤率）影響：

1. **Jakes model**：\$c\[k] = h(t)|\_{t=kT}\$
2. **Gaussian Rayleigh model**：\$c\[k] = n\_R\[k] + jn\_I\[k]\$

並畫出 BER vs. Eb/No 曲線，觀察性能差異。

### 🧪 實驗觀察與討論 (b-i)

1. Jakes Model 對應到 SNR = 15 時 BER 約為 \$10^{-2}\$
2. Gaussian Rayleigh 對應 SNR = 14 即可達相同 BER
   → 代表 Jakes model **在高 SNR 效能略優**

📌 可能原因：

* Jakes Model 中的 channel gain 是時間相關的，對同一筆資料的不同路徑有相干性
* Gaussian Rayleigh 是獨立生成的雜訊，對接收端來說有更高的變異性 → 錯誤率稍高

### 🌟 K-Factor = 5 實驗 (b-ii)

Rician 通道情況：K=5，表示有強直視徑（LOS）分量

**N=50 時觀察：**

* Jakes Rayleigh (K=5) 與 Gaussian Rayleigh (K=5) 非常接近

**N=100 時觀察：**

* Gaussian Rayleigh (K=5) 略優

📌 解釋：K 增加 → Rician 分佈趨近高斯 → Rayleigh 模型也能模擬出類似效果

### 🔁 (iii) Diversity Order 實驗 (L=2, L=3)

圖中展示多路徑接收時 (L=2,3)，使用 MRC（Maximal Ratio Combining）合併接收：

* L=2 時：兩模型表現差不多
* L=3 時：Jakes model 優於 Gaussian Rayleigh

📌 解釋：

* **Diversity order** 越高 → 越能減少 deep fading 機率
* Jakes Model 因為有結構性相干性 → 在 L 增加時更能展現出優勢

---

## ✍️ MRC 筆記

### 📌 定義

MRC（最大比率合併）是一種合併多天線或多路徑訊號的方法，將各通道的接收訊號乘上其共軛增益再相加：

$$
y = \sum_{i=1}^{L} h_i^* x_i
$$

其中 \$h\_i\$ 是通道增益，\$x\_i\$ 是第 \$i\$ 路的接收訊號。

### ✅ 優點

* 可最佳化 SNR
* 增強系統可靠度與穩定性
* 在 fading channel 中尤其有效（特別是 Rayleigh、Rician）

### 🔁 實驗中的應用

```matlab
% MRC 合併決策
% 對於 Jakes 模型與 Gaussian Rayleigh 模型都使用

% 決策前處理：使用 conj(h) * x
decision_k = real(conj(c_k1).*x_k1 + conj(c_k2).*x_k2);
decision_k_Gauss = real(conj(c_k1).*x_k1_Gauss + conj(c_k2).*x_k2_Gauss);
```

---

## ✅ 總結表格：比較各項實驗結果

| 實驗項目         | 通道模型比較                     | 優勢模型          | 說明                      |
| ------------ | -------------------------- | ------------- | ----------------------- |
| b(i)         | Jakes vs Gaussian Rayleigh | Jakes (高 SNR) | 有時間相干性，fading 平穩        |
| b(ii) K=5    | Jakes vs Gaussian (Rician) | Gaussian 略優   | Rician 趨近 Gaussian，差異縮小 |
| b(iii) L=2,3 | 多路徑 MRC                    | L=3 時 Jakes 優 | 多分支下 Jakes 穩定性更高        |

