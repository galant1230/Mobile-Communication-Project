# (II) 實驗二：BER 曲線繪製 (理想通道)

## 先備知識 - 匹配濾波器

### 🔹 等效基頻模型與匹配濾波 (Matched Filter)

#### 1. RF (實際訊號) 模型

在實際無線通道中，傳送端送出的訊號形式為：

$s(t) = \sqrt{\tfrac{2E_b}{T_b}} \cos(2\pi f_c t + \phi)$

* \$f\_c\$：載波頻率
* \$T\_b\$：一個 bit 的持續時間
* \$E\_b\$：每個 bit 的能量
* \$\phi\$：相位 (BPSK 時是 0 或 \$\pi\$)

這是實際物理層送出的「連續時間 + 載波」訊號。

---

#### 2. 接收端處理：解調 + 匹配濾波

接收端拿到：
$r(t) = s(t) + n(t)$
其中 \$n(t)\$ 是 AWGN 雜訊。

我們要判斷這個符號是 bit 0 還是 bit 1。

如果只是「隨便取一個瞬間值」來判斷，會很不準，因為雜訊在每個時間點都不同。

👉 所以我們要想辦法「整合整個 bit 期間的能量」，把訊號累積起來，這樣雜訊的影響相對變小。

處理步驟：

1. **解調 (Demodulation)**：把訊號搬到基頻，消除載波 cos 成分。
2. **匹配濾波 (Matched Filter)**：
   * 匹配濾波器是一個線性濾波器，脈衝響應設定為 傳送符號波形的時間反轉版：
      * 濾波器脈衝響應 \$h(t) = s(T\_b - t)\$。
      * 目標：最大化「取樣時刻的訊號對雜訊比 (SNR)」。
      * 等效於把整個 bit 的能量積分起來，雜訊部分因隨機性而部分抵消。
---
### 🔹 為什麼匹配濾波器能產生最大相關

#### 2.1. 濾波器輸出定義

接收端濾波器輸出：
$y(t) = \int r(\tau) \ h(t-\tau) \ d\tau$

若選擇濾波器脈衝響應：
$h(t) = s(T_b - t)$

#####  2.1.1 匹配濾波數學推導

- 代入 \$h(t)\$

匹配濾波器的脈衝響應：
$h(t) = s(T_b - t)$

代入卷積公式：
$y(t) = \int r(\tau) \ s(T_b - (t - \tau)) \ d\tau$

- 取樣在 \$t = T\_b\$

關注取樣時刻 \$t = T\_b\$：
$y(T_b) = \int r(\tau) \ s(T_b - (T_b - \tau)) \ d\tau$

裡面可以化簡：
$s(T_b - (T_b - \tau)) = s(\tau)$

因此：
$y(T_b) = \int r(\tau) \ s(\tau) \ d\tau$

#### ✅ 結果意義

* 在 \$t=T\_b\$ 時，濾波器輸出 = 接收訊號 \$r(t)\$ 與原始訊號 \$s(t)\$ 的 **內積 (相關值)**。
* 若 \$r(t) = s(t)\$（無雜訊），輸出 = 訊號能量 \$E\_b\$。
* 有雜訊時，輸出 = 訊號相關值 + 雜訊項。

則在 \$t=T\_b\$ 時：
$y(T_b) = \int r(\tau) \ s(\tau) \ d\tau$

---

#### 2. 2. 與互相關的關係

* 互相關定義：
  $(r \star s)(t) = \int r(\tau) \ s(\tau+t) \ d\tau$
* 當 \$t=0\$：
  $(r \star s)(0) = \int r(\tau) \ s(\tau) \ d\tau$

👉 可以看出，當 \$t=T\_b\$，匹配濾波器輸出其實就是 \$r\$ 與 \$s\$ 的 **相關值**。

---

#### 2. 3. 為什麼是最大化？

若輸入 \$r(t) = s(t)\$（無雜訊）：
$y(T_b) = \int_0^{T_b} s(\tau) s(\tau) \ d\tau = \int_0^{T_b} |s(\tau)|^2 \ d\tau$

這就是 **訊號能量**：
$y(T_b) = E_b$

* 當訊號與模板完全對齊 → 輸出 = 最大能量。
* 若不是同一個符號 (或有雜訊/相位差) → 相關值會小於 \$E\_b\$。

---

#### 2. 4. 直覺理解

* 濾波器輸出 = 測量「輸入與模板的相似度」。
* 設定 \$h(t) = s(T\_b - t)\$ → 在取樣時刻得到互相關值。
* 當輸入真的是 \$s(t)\$ → 相關值 = 能量最大。
* 所以匹配濾波器能在取樣點 **最大化 SNR**：訊號被累積最大，雜訊因隨機性部分抵消。

---

#### 3. 等效基頻模型

經過下變頻 + 匹配濾波後，訊號簡化為：

$y = \pm \sqrt{E_b} + n, \quad n \sim \mathcal{N}(0, N_0/2)$

* \$+\sqrt{E\_b}\$：對應 bit = 0
* \$-\sqrt{E\_b}\$：對應 bit = 1
* \$n\$：經過濾波後的等效高斯雜訊

👉 這就是 **等效基頻模型**：不再需要 cos 波形，只保留符號能量與雜訊。

---

#### 4. 為什麼符號是 ±√Eb？

* 定義：每 bit 能量 \$E\_b = A^2 T\_b\$。
* 假設 \$T\_b=1\$，則 \$A = \sqrt{E\_b}\$。
* 因此符號點必須是：
  $s \in \{+\sqrt{E_b}, -\sqrt{E_b}\}$
* 若設定 \$E\_b=1\$，就簡化為 ±1。

---

#### 5. 匹配濾波的目的

1. **最大化訊號能量**：把一整個 bit 的能量集中到一個取樣點。
2. **最小化雜訊干擾**：在 AWGN 通道下可證明是最佳化的。
3. **轉換形式**：把連續時間訊號變成單一數值，方便判決 bit。

---

### ✅ 總結

* **RF 實際模型**：用 cos 表示，因為硬體真的在送載波。
* **等效基頻模型**：化簡成 ±√Eb，方便數學分析 BER。
* **匹配濾波**：最佳化接收器設計，確保在取樣點 SNR 最高。

---
```matlab
%---------------------先宣告題目所需的條件----------------------------%
clear all; close all; clc;
M = 64;                 % Number of FFT points
sc = 52;                % Number of data sub-carriers
ofdmbit = 52;           % Number of bits per OFDM symbol
                        % 由於是使用BPSK來調變，因此BPSK的sub-carrier也是52
Nsymbol = 100;          % Number of symbols

%---------------------使用BPSK symbol進行調變--------------------------%
Data = randi([0,1], ofdmbit*Nsymbol,1);   % 產生隨機bits
dk   = 2*Data - 1;                        % BPSK調變: 0→-1, 1→+1
dk_sym = reshape(dk, sc, Nsymbol);        % Serial → Parallel

% 把 subcarrier 放入題目所給的 data subcarrier
% 注意在 data subcarrier 上 32 與 34 之間沒有第 33 號 sub-carrier
Dk = zeros(M, Nsymbol);
Dk(7:32, :)  = dk_sym(1:26, :);
Dk(34:59, :) = dk_sym(27:52, :);

%--------------------- M-point IDFT (轉換到時域) -----------------------%
IDFT = [];
for i = 1:Nsymbol
    IDFT = [IDFT ifft(Dk(:, i), M)];
end

%--------------------- adding 16 samples (Cyclic prefix) ----------------%
D_cp = zeros(M+16, Nsymbol);
D_cp(1:16, :)   = IDFT(49:64, :);
D_cp(17:80, :)  = IDFT(1:64, :);

%--------------------- receiving signal z[n] = x[n]*h[n] + n[n] ---------%
h = 1;

%--------------------- parallel to serial conversion --------------------%
x_n = reshape(D_cp, 1, (M+16)*Nsymbol);
xh  = conv(x_n, h);

%--------------------- awgn addition ------------------------------------%
Eb = mean(abs(x_n).^2);
for Eb_N0 = 0:2:8
    % theoretical BER
    BPSK_Pb(Eb_N0/2+1) = 0.5*erfc(sqrt(10^(Eb_N0/10)));

    Np = Eb * 10^(-(Eb_N0/10));
    ns = sqrt((64/52)*Np/2) * (randn(1,length(D_cp(:,:))*Nsymbol) ...
       + 1i*randn(1,length(D_cp(:,:))*Nsymbol));
    z_nxh = xh(1:end) + ns;

    %----------------- cyclic prefix removal ----------------------------%
    receiver = reshape(z_nxh, 80, Nsymbol);
    rn = receiver(17:end, :);

    %----------------- M-point DFT -------------------------------------%
    DFT = [];
    for i = 1:Nsymbol
        DFT = [DFT fft(rn(:, i), M)];
    end

    %----------------- pass through equalizer --------------------------%
    Yk64 = [DFT(7:32,:) ; DFT(34:59,:)];
    Yk52 = [Yk64(:,1).' Yk64(:,2).'].' ;

    %----------------- scatterplot (Before Equalizer) ------------------%
    scatterplot(Yk52);

    H_n = fft(h, M);
    Dk = [];
    for i = 1:Nsymbol
        Dk = [Dk DFT(:, i)./H_n.'];
    end
    Dkreceiver = [Dk(7:32,:) ; Dk(34:59,:)];

    %----------------- serial to parallel conversion -------------------%
    Dkr  = reshape(Dkreceiver, 1, ofdmbit*Nsymbol);
    Dkre = sign(real(Dkr));

    ber(Eb_N0/2+1) = length(find(Dkre ~= dk)) / (52*Nsymbol);
end

%----------------- 畫 BER 曲線 ------------------------------------------%
semilogy((0:2:8), ber, 'r-o');
hold on;
semilogy((0:2:8), BPSK_Pb, 'b-*');
hold on;
grid on;
title('Bit error probability curve for BPSK using OFDM');
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('Simulation BER','Theoretical BER');

%----------------- scatterplot (After Equalizer) -----------------------%
scatterplot(Dkr(1,1:104));
```

---
## 1. AWGN
```matlab
Eb = mean(abs(x_n).^2);
for Eb_N0 = 0:2:8
    % theoretical BER
    BPSK_Pb(Eb_N0/2+1) = 0.5*erfc(sqrt(10^(Eb_N0/10)));

    Np = Eb * 10^(-(Eb_N0/10));
    ns = sqrt((64/52)*Np/2) * (randn(1,length(D_cp(:,:))*Nsymbol) ...
       + 1i*randn(1,length(D_cp(:,:))*Nsymbol));
    z_nxh = xh(1:end) + ns;
```

---

### 1.1 Eb

程式逐行解釋

```matlab
Eb = mean(abs(x_n).^2);
```

* `x_n` 是發送端 OFDM 訊號（IFFT + CP 後展平的一維訊號）。
* `abs(x_n).^2`：每個取樣點的能量。
* `mean(...)`：取平均能量。
* 在 BPSK 下，這等效於 **每 bit 能量 \$E\_b\$**。

---

```matlab
for Eb_N0 = 0:2:8
```

* 設定模擬的 \$E\_b/N\_0\$ 測試範圍：0, 2, 4, 6, 8 dB。
* 用來觀察 BER 隨 SNR 的變化。

---
### 1.2 BPSK Pre-Knowledge 
### 🔹 BPSK 符號能量與一個 bit 的關係
#### 1.2.1. BPSK 符號 = 一個 bit

在 BPSK (Binary Phase Shift Keying) 中：

* bit = 0 → \$+1\$（或 \$+\sqrt{E\_b}\$ ）
* bit = 1 → \$-1\$（或 \$-\sqrt{E\_b}\$ ）

👉 每個符號只代表 **一個 bit**，不是兩個 bit。
（對比 QPSK：一個符號代表 2 bits，因為星座圖有 4 個點。）

---

#### 1.2.2. 為什麼寫成 ±√Eb？

我們希望「每個 bit 的能量 = \$E\_b\$ 」。

能量定義：
$E_b = \int_{0}^{T_b} |s(t)|^2 dt$

假設符號在整個 bit 期間振幅固定 = \$A\$：
$E_b = A^2 \cdot T_b$

若取 \$T\_b = 1\$，就有：
$A^2 = E_b \quad \Rightarrow \quad A = \sqrt{E_b}$

因此 BPSK 的符號點應該是：
$s \in \{+\sqrt{E_b}, -\sqrt{E_b}\}$

教科書常簡化設定 \$E\_b = 1\$，這樣就寫成 ±1。
👉 此時 ±1 的能量 = 1，也就是「每 bit 的能量 = 1」。

---

#### 1.2.3. 為什麼不是兩個 bit？

* BPSK 只有 **兩個狀態 (±)**，所以只能表示 2 個可能值 → 也就是 **1 個 bit**。
* 若要表示 2 bits 需要 4 個點（例如 QPSK，星座圖在四象限）。

### 1.2. ✅ 總結

* BPSK：一個符號 = 1 個 bit，符號點 = ±√Eb。
* 若 \$E\_b = 1\$，則符號點 = ±1。
* ±1 並不是兩個 bit，而是 **同一個 bit 的兩種可能值**。

---

### 1.2.4 BPSK 在 AWGN 下的理論 BER

在 **AWGN 通道**下，BPSK 的錯誤率公式為：

(i). **系統模型**：

   * 傳送 bit = 0 → \$+\sqrt{E\_b}\$
   * 傳送 bit = 1 → \$-\sqrt{E\_b}\$
   * 接收端：
     
     $y = \pm\sqrt{E_b} + n, \quad n \sim \mathcal{N}(0, N_0/2)$

(ii). **錯誤機率**：

   使用零門檻判決，錯誤率為：
   
   $P_b = Q\!\left(\sqrt{\tfrac{2E_b}{N_0}}\right)$

(iii). **Q-function 與 erfc 的關係**：
   
   $Q(x) = \tfrac{1}{2}\,\text{erfc}\!\left(\tfrac{x}{\sqrt{2}}\right)$

   因此：
   
   $P_b = \tfrac{1}{2}\,\text{erfc}\!\left(\sqrt{\tfrac{E_b}{N_0}}\right)$

---

### 🔹 MATLAB 程式對應

```matlab
% theoretical BER
BPSK_Pb(Eb_N0/2+1) = 0.5*erfc(sqrt(10^(Eb_N0/10)));
```

* `10^(Eb_N0/10)`：把 dB 轉換成線性比例。
* `sqrt(...)`：計算 \$\sqrt{E\_b/N\_0}\$ 。
* `erfc(...)`：互補誤差函數，對應公式。
* `0.5*erfc(...)`：對應 \$P\_b = \tfrac{1}{2},\text{erfc}(\sqrt{E\_b/N\_0})\$ 。

---

### 🔹 為什麼寫 `BPSK_Pb(Eb_N0/2+1)`？

迴圈設定：

```matlab
for Eb_N0 = 0:2:8
```

* `Eb_N0` 代表真實的 SNR(dB)：0, 2, 4, 6, 8。
* MATLAB 陣列索引 **不能從 0 開始**，所以不能直接用 `BPSK_Pb(Eb_N0)`。

解法：

* 除以 2：把 {0,2,4,6,8} → {0,1,2,3,4}
* 再 +1：變成 {1,2,3,4,5}

這樣剛好對應陣列索引，結果存在 `BPSK_Pb(1:5)` 裡。


---

```matlab
Np = Eb * 10^(-(Eb_N0/10));
```

* 根據公式 \$ N\_0 = E\_b / (E\_b/N\_0)\$ ，計算 **雜訊功率**。
* 在這裡 `Np` 代表 AWGN 的平均功率。

---

```matlab
ns = sqrt((64/52)*Np/2) * (randn(1,length(D_cp(:,:))*Nsymbol) ...
       + 1i*randn(1,length(D_cp(:,:))*Nsymbol));
```

* 產生複數 AWGN：`randn + 1i*randn`。
* 除以 `sqrt(2)` 讓實部、虛部各自有變異數 \$Np/2\$，總功率為 \$Np\$。
* 前面的 `(64/52)`：修正因子，因為雖然有 64 個子載波，但只有 52 個用來傳資料，需做功率正規化。
* 結果 `ns` 為符合功率分配的 AWGN。

---

```matlab
z_nxh = xh(1:end) + ns;
```

* `xh`：通過通道的訊號 (這裡 h=1，所以相當於原始訊號)。
* `ns`：生成的 AWGN。
* 兩者相加 → 接收端訊號 `z_nxh`。



