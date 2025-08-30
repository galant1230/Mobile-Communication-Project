# 📡 Mobile Communication Course Project #2

## 系統架構介紹

以下為本專案的 **OFDM 收發機 (Transceiver)** 系統架構圖：  

![OFDM System Block Diagram](1.jpg)

---

### 🔼 發射端 (Transmitter)

1. **Modulator**  
   - 將輸入的 bits 調變為 BPSK 符號：
 
$$
D_k \in \{-1, +1\}
$$  

2. **M-Point IDFT (IFFT)**  
   - 將頻域符號轉換至時域，得到 $D[n]$。  
   - 目的：讓多個子載波保持正交，避免干擾。  

3. **Add cyclic extension (CP)**  
   - 在符號前加上循環字首 (Cyclic Prefix)，形成 $x[n]$。  
   - 功能：將通道中的「線性卷積」轉為「循環卷積」，抵抗符號間干擾 (ISI)。  

4. **I-Q Mod. DAC**  
   - 轉換為 I/Q 類比訊號（模擬中可省略）。  

5. **Channel $h[n]$**  
   - 模擬無線通道效應，包括多路徑衰落與加性白雜訊 (AWGN)。  

---

### 🔽 接收端 (Receiver)

1. **I-Q Mod. ADC**  
   - 將類比訊號轉回數位訊號，得到受通道影響的 $z[n]$。  

2. **Remove cyclic extension**  
   - 去掉循環字首，還原有效的 OFDM 符號，得到 $r[n]$。  

3. **M-Point DFT (FFT)**  
   - 將時域訊號轉回頻域，得到：
    
$$
Y_k = H_k D_k + N_k
$$  

4. **Frequency domain equalizer**  
   - 使用 Zero-Forcing 等化器：
   
$$
\hat{D}_k = \frac{Y_k}{H_k}
$$  

   - 目的：補償通道的增益與相位失真。  

5. **Demodulator**  
   - 將等化後的符號 $\hat{D}_k$ 判決回原始 bits。  

---
## Code

```matlab
%% Mobile Communication Project #2 - OFDM Simulation
clear; clc; close all;

%% 系統參數
M = 64;                 % FFT 點數
sc = 52;                % 資料子載波數
ofdm_bit = 52;          % 每個 OFDM symbol 的 bits 數 (BPSK = 52)
Nsymbol = 100;          % 模擬的 OFDM 符號數
cp_len = 16;            % 循環字首長度

%% 產生隨機 bits 並做 BPSK 調變
Data = randi([0 1], ofdm_bit * Nsymbol, 1);
dk = 2*Data - 1;                            % BPSK: 0→-1, 1→+1
dk_sym = reshape(dk, sc, Nsymbol);          % Serial → Parallel

%% 插入子載波 (52 個資料子載波)
Dk = zeros(M, Nsymbol);
Dk(7:32,:)  = dk_sym(1:26,:);               % 左 26 個子載波
Dk(34:59,:) = dk_sym(27:52,:);              % 右 26 個子載波

%% IFFT (轉回時域)
IDFT = ifft(Dk, M);

%% 加入循環字首 (CP)
D_cp = [IDFT((M-cp_len+1):M,:); IDFT];      % 前 16 點複製到最前面
x_n = reshape(D_cp, 1, (M+cp_len)*Nsymbol);

%% 通道模型 (多路徑衰落)
h = [0.5-0.5j, 0, 0.15+0.12j, 0, 0, -0.1+0.05j];
h = h ./ norm(h);                           % 正規化
xh = conv(x_n, h);

%% 加入 AWGN (Eb/N0 = 20 dB)
Eb_N0 = 20;
Eb = mean(abs(x_n).^2);
Np = Eb*10^(-Eb_N0/10);
z = awgn(xh(1:end-5), Eb_N0, 'measured');

%% 接收端：移除循環字首 (CP)
receiver = reshape(z, M+cp_len, Nsymbol);
rn = receiver(cp_len+1:end,:);              % 去掉 CP

%% FFT (轉回頻域)
DFT = fft(rn, M);

%% 取出 52 個有效子載波 (等化前)
Yk = [DFT(7:32,:); DFT(34:59,:)];
Yk_serial = reshape(Yk, 1, sc*Nsymbol);

%% 繪製等化前星座圖
figure;
scatterplot(Yk_serial(1:104));              % 取前 2 個 OFDM symbol
title('Constellation Before Equalization');

%% 頻域等化 (Zero-Forcing)
Hn = fft(h, M);
Dk_eq = DFT ./ Hn.';                        % Zero-Forcing: Yk / Hk

%% 取出 52 個子載波 (等化後)
Dk_receiver = [Dk_eq(7:32,:); Dk_eq(34:59,:)];
Dkr_serial = reshape(Dk_receiver, 1, sc*Nsymbol);

%% 繪製等化後星座圖
figure;
scatterplot(Dkr_serial(1:104));
title('Constellation After Equalization');
```
---
## 1. OFDM 系統參數與公式說明

```matlab
M = 64;                 % FFT 點數
sc = 52;                % 資料子載波數
ofdm_bit = 52;          % 每個 OFDM symbol 的 bits 數 (BPSK = 52)
Nsymbol = 100;          % 模擬的 OFDM 符號數
cp_len = 16;            % 循環字首長度
```

---

### 1.1 FFT 點數 \$M = 64\$

在 OFDM 中，一個符號是透過 \$M\$ 點 IFFT 轉換成時域訊號。公式為：

$$
x[n] = \frac{1}{M} \sum_{k=0}^{M-1} D_k \, e^{j2\pi kn/M}, \quad n=0,1,\dots,M-1
$$

* **\$M=64\$** → 表示每個 OFDM 符號包含 **64 個子載波 (subcarriers)**。
* **\$D\_k\$** → 第 \$k\$ 個子載波上的符號。

#### 🔹 一個 OFDM Symbol 是什麼？

* **Symbol**：在數位通訊裡，symbol 指的是 **調變後的一個訊號單位**。

  * 例如：BPSK → 1 bit 對應 1 symbol。
  * QPSK → 2 bits 對應 1 symbol。
  * 16-QAM → 4 bits 對應 1 symbol。

* **在 OFDM 中**：
  一個 OFDM symbol 並不是單純的一個調變符號，而是 **一整個 IFFT 輸出的時域訊號 (長度 \$M\$)**。

  * 它由 \$M\$ 個子載波組合而成。
  * 每個子載波上承載一個調變符號 (這裡是 BPSK)。

👉 換句話說：**一個 OFDM Symbol = 所有子載波的調變符號同時送出的一個時域區塊**。

#### 🔹 為什麼一個 Symbol 是 64 個子載波？

因為我們使用 \$M=64\$ 點 IFFT，頻域有 \$64\$ 個「bin」位置。

* 每個 bin 代表一個子載波 (subcarrier)。
* 把資料 \$D\_k\$ 放進去後，做 IFFT → 同時轉成一段時域波形。

所以：

* 一個 OFDM 符號就是 64 個子載波「同時傳送」的組合波。
* 這 64 個子載波中，只有 **52 個載波承載資料** (\$sc=52\$)。
* **其餘 12 個子載波用途**：

  * **DC 子載波 (1 個)**：位於頻譜中心 (\$k=0\$)，通常設為 0，不承載資料，避免直流偏移 (DC offset) 造成干擾。
  * **Guard Bands (11 個)**：分布在頻譜兩端，用來隔離相鄰頻道，避免邊界干擾 (Adjacent Channel Interference)。

#### 🔹 為什麼 BPSK 調變後是在頻域？

* 在傳統單載波系統中，BPSK 調變會直接變成時域波形 (例如 ±1 的方波)。
* **但在 OFDM 中，我們先把 BPSK 調變後的符號放到頻域的子載波位置 (\$D\_k\$)。**
* 接著利用 IFFT 把整個頻域訊號轉換成時域訊號 → 這才是實際要發送的 OFDM 符號。

所以流程是：

1. **比特 → 符號 (BPSK 調變)**：把 bits 映射成 ±1 (這時還只是「符號」，還沒分配頻率)。
2. **符號 → 子載波 (\$D\_k\$)**：把這些 ±1 塞進頻域的子載波 bin。
3. **頻域 → 時域 (IFFT)**：把整個頻域訊號轉成時域訊號，這就是要送出去的 OFDM symbol。

### ✅ 總結：

* 一個 **OFDM Symbol** 是一整段時域訊號，長度 \$M=64\$。
* 它由 **64 個子載波**組成，其中只有 52 個承載資料。
* 另外 **1 個 DC 子載波不承載資料**，**兩端的 Guard bands 隔離鄰頻**。
* BPSK 調變後的 ±1 先放在頻域子載波位置，再經過 IFFT 轉成時域，才會送到通道。

---

### 1.2. 資料子載波數 $sc = 52$

在 $M=64$ 的子載波裡，只有 **52 個承載資料**，其餘 **12 個作為保護子載波 (Guard bands) 與 DC**。  

---

### 1.3 每個 OFDM Symbol 的 Bits 數（BPSK）

因為採用 BPSK 調變，每個子載波承載 1 bit，對應到頻域符號集合：D\_k ∈ {-1, +1}。因此在本模擬（52 個資料子載波）下：

* 每個 OFDM symbol 的 bits = sc × 1 = 52

一般化：若改成 M\_mod-QAM，則每個 OFDM symbol 的 bits = sc × log2(M\_mod)（例如 QPSK=2、16-QAM=4）。

和程式對應：

* sc = 52，ofdm\_bit = 52（BPSK → 每個子載波 1 bit）。
* dk = 2\*Data - 1 進行 BPSK 映射（0→-1、1→+1）。
* `Dk(7:32,:)` 與 `Dk(34:59,:)` 將 52 個資料載波填入頻域（跳過 DC=索引 33，並在兩端留 Guard bands）。

---

### 1.4 符號數 Nsymbol = 100

模擬傳送 100 個 OFDM 符號：

* 總 bits 數 N\_bits = ofdm\_bit × Nsymbol = 52 × 100 = 5200。

同時計算時間長度與開銷：

* 每個 OFDM 符號（含 CP）取樣數：M + cp\_len = 64 + 16 = 80。
* 總取樣數（含 CP）：80 × 100 = 8000。
* 總取樣數（去 CP）：64 × 100 = 6400。
* CP 開銷：cp\_len / M = 16 / 64 = 25%（亦即有效頻寬/吞吐下降 25%）。

注意：若有導頻（pilots）、通道編碼、或子載波關閉，實際有效資料速率會再降低。本模擬未使用導頻，直接以 52 個資料子載波傳送。

---

### 1.5 循環字首長度 $len_{cp} = 16$

在時域 OFDM 符號前，會複製最後 $len_{cp} = 16$ 個樣本加到最前面：

$$
x_{cp}[n] =
\begin{cases}
x[n + M - len_{cp}], & 0 \leq n < len_{cp} \\
x[n - len_{cp}], & len_{cp} \leq n < M + len_{cp}
\end{cases}
$$

#### 🔹 CP 的目的與條件

1. **避免符號間干擾 (ISI)**

   * 若通道最長延遲擴散 (最大有效路徑延遲) 為 \$L-1\$，則只要
     $len_{cp} \geq L-1$
     即可確保前一個符號的尾巴不會干擾到目前符號。

2. **把線性卷積變成循環卷積**

   * 去掉 CP 後的 \$M\$ 點區段滿足循環卷積關係。
   * 在頻域可以寫成：
     $Y[k] = H[k] \cdot D[k] + W[k]$
   * 其中：

     * \$H\[k]\$ = 通道頻率響應 (對 \$h\[\ell]\$ 做 \$M\$ 點 DFT)。
     * \$D\[k]\$ = 發送符號。
     * \$W\[k]\$ = 雜訊。
   * 這使得 **每個子載波可獨立等化**。

---
## 2. 隨機 bits 產生與 BPSK 調變
```matlab
%% 產生隨機 bits 並做 BPSK 調變
Data = randi([0 1], ofdm_bit * Nsymbol, 1);    
% randi([imin imax], m, n)：產生一個 m×n 的隨機整數矩陣
% 這裡 randi([0 1], ofdm_bit * Nsymbol, 1) → 產生 (5200×1) 的隨機 0/1 序列
% (5200 = 52 × 100)

dk = 2*Data - 1;                            
% BPSK 映射：0 → -1, 1 → +1
% 映射結果是一串 {-1, +1} 的序列（實數）

dk_sym = reshape(dk, sc, Nsymbol);          
% reshape(A, m, n)：將矩陣 A 重新排成 m×n 大小（不改變元素總數）
% 這裡 reshape(dk, 52, 100) → 把長度 5200 的向量，重組成 52×100 矩陣
% → 每一欄 (52×1) 對應一個 OFDM symbol（52 個子載波）
```
---
## 3. 插入子載波 (52 個資料子載波)
```matlab
%% 插入子載波 (52 個資料子載波)
Dk = zeros(M, Nsymbol);      
% 先建立一個 M×Nsymbol 的零矩陣 (64×100)
% → 相當於頻域有 64 個子載波位置 (每欄是一個 OFDM symbol)

Dk(7:32,:)  = dk_sym(1:26,:);    
% 將前 26 個資料符號填入索引 7–32 的子載波
% → 這是「左半邊」的 26 個子載波 (避開最左邊 6 個 Guard band)

Dk(34:59,:) = dk_sym(27:52,:);   
% 將後 26 個資料符號填入索引 34–59 的子載波
% → 這是「右半邊」的 26 個子載波 (避開索引 33 的 DC 與最右 5 個 Guard band)
```
---
## 4. IFFT (轉回時域)

```matlab
%% IFFT (轉回時域)
IDFT = ifft(Dk, M);
```

### 🔹 程式解釋

* `ifft(X, M)`：對矩陣 `X` 的每一欄做 \$M\$ 點 **逆快速傅立葉轉換 (IFFT)**。
* 這裡 `Dk` 是 \$64 \times 100\$ 的頻域資料矩陣：

  * 每一欄代表一個 OFDM symbol 在頻域的 64 個子載波。
* `ifft(Dk, M)` 會輸出一個 \$64 \times 100\$ 的矩陣 `IDFT`：

  * 每一欄是經過 \$M=64\$ 點 IFFT 後的**時域 OFDM symbol**。

### 🔹 為什麼 IFFT = IDFT？

在數學定義上，**\$M\$ 點逆離散傅立葉轉換 (IDFT)** 定義為：

$$
x[n] = \frac{1}{M}\sum_{k=0}^{M-1} D[k] e^{j2\pi kn/M}, \quad n=0,1,\dots,M-1
$$

而 MATLAB 的 `ifft` 函數計算公式正好就是：

$$
\texttt{ifft}(D) = \frac{1}{M} \cdot IDFT(D)
$$

因此 **`ifft` 實現的就是數學上的 IDFT**，並且已經自動包含了 \$1/M\$ 的歸一化因子。

### 🔹 小範例

```matlab
X = [1 1 1 1];        % 頻域輸入
x = ifft(X, 4)        % 做 4 點 IFFT
```

結果：

```
x = 0.25    0.25    0.25    0.25
```

這對應數學公式：

* 因為所有頻域係數都 =1 → 時域輸出為常數序列。

✅ 總結：
* **數學上的 IDFT** 與 **MATLAB 的 `ifft`** 完全等價。
* 在 OFDM 中，`ifft` 的作用就是把頻域的 \$D\_k\$（子載波符號）轉換成一段時域 OFDM 符號 \$x\[n]\$。

---
## 1.5 加入循環字首 (\$len\_{cp} = 16\$)

```matlab
%% 加入循環字首 (CP)
D_cp = [IDFT((M-cp_len+1):M,:); IDFT];   % 前 16 點複製到最前面
x_n = reshape(D_cp, 1, (M+cp_len)*Nsymbol);
```

### 🔹 程式解釋

1. **`IDFT((M-cp_len+1):M,:)`**

   * 從每個 OFDM 符號 (每一欄) 的尾端，取出最後 \$len\_{cp}\$ 個樣本。
   * 這裡 \$len\_{cp} = 16\$，所以取出第 49–64 點。
   * 注意：因為 \$M=64\$，所以 index `49:64` 剛好在範圍內，不會有「沒值」的情況。

2. **`[IDFT((M-cp_len+1):M,:); IDFT]`**

   * 把取出的 16 點樣本放到原本的符號前面，形成「循環字首 (Cyclic Prefix)」。
   * 新的每欄長度 = \$64+16=80\$ 點。

3. **`reshape(D_cp, 1, (M+cp_len)*Nsymbol)`**

   * 將矩陣展平成一列向量，方便後續做通道傳輸。
   * 長度 = \$(M+len\_{cp}) \times N\_{symbol} = 80 \times 100 = 8000\$。

### 🔹 為什麼要加 CP？

* **避免符號間干擾 (ISI)**：確保多路徑延遲不會影響到下一個符號。
* **把線性卷積轉換為循環卷積**：

  * 通道卷積在時域表現為循環結構。
  * 去掉 CP 後，FFT 可以把接收信號還原成：
    $Y[k] = H[k]\,D[k] + W[k]$
  * 每個子載波獨立等化，簡化接收端設計。

✅ 總結：

* \$len\_{cp}\$：循環字首長度，這裡為 16。
* `D_cp`：在每個 OFDM symbol 前面加上循環字首。
* `x_n`：把所有帶 CP 的符號展平成傳輸序列。
* 這一步確保 OFDM 能在多路徑通道下維持子載波正交性。

---
## 1.6 通道模型 (多路徑衰落)

```matlab
%% 通道模型 (多路徑衰落)
h = [0.5-0.5j, 0, 0.15+0.12j, 0, 0, -0.1+0.05j];
h = h ./ norm(h);                           % 正規化
xh = conv(x_n, h);                          % 通道卷積
```

### 🔹 通道脈衝響應 \$h\$

* `h` 表示通道的 **脈衝響應 (Impulse Response)**，即多路徑效應：

  * \$h\[0] = 0.5 - 0.5j\$ → **主路徑 (直射路徑)**，能量最強並帶有相位旋轉。
  * \$h\[1] = 0\$ → 此延遲沒有路徑。
  * \$h\[2] = 0.15 + 0.12j\$ → **次要反射路徑**，能量較弱並帶有不同相位。
  * \$h\[3] = 0\$, \$h\[4] = 0\$ → 沒有能量的延遲位置。
  * \$h\[5] = -0.1 + 0.05j\$ → **更長延遲的路徑**，能量最弱。

👉 這樣的設定模擬了一個 **6 taps 離散多路徑通道**，只有第 0、2、5 tap 有效。

### 🔹 為什麼要正規化？

```matlab
h = h ./ norm(h);
```

* `norm(h)` = \$\sqrt{\sum |h\[n]|^2}\$ = 通道能量。
* 除以 `norm(h)` → 讓通道的平均能量 = 1。
* 功能：
  * 避免通道能量過大或過小影響 SNR 控制。
  * 確保接收端的功率尺度與發射端一致。

### 🔹 通道卷積

```matlab
xh = conv(x_n, h);
```

* 在無線通訊中，通道作用相當於 **線性卷積**：
  $y[n] = (x * h)[n]$
* MATLAB `conv` → 對傳輸信號 `x_n` 與通道 `h` 做卷積。
* 輸出 `xh` = 信號經過多路徑通道後的時域波形。
  
✅ 總結：

* `h`：模擬 3 條有效路徑的多路徑通道。
* `h./norm(h)`：能量正規化，方便控制 SNR。
* `conv(x_n,h)`：模擬信號通過多路徑通道的效果。

---
## 1.7 加入 AWGN (雜訊)

```matlab
Eb_N0 = 20;                        % 設定 Eb/N0 = 20 dB
Eb = mean(abs(x_n).^2);            % 平均 bit 能量 (BPSK 下符號能量 = bit 能量)
Np = Eb*10^(-Eb_N0/10);            % 由 Eb/N0 換算得到雜訊功率

z = awgn(xh(1:end-5), Eb_N0, 'measured');
% 在通道輸出訊號上加入 AWGN
% 'measured' 表示依據訊號的實際功率來決定雜訊能量
```

### 🔹 MATLAB `awgn` 函數用法

* 語法：

  ```matlab
  y = awgn(x, SNRdB, 'measured')
  ```

  * `x`：輸入訊號。
  * `SNRdB`：目標訊號雜訊比 (dB)。
  * `'measured'`：自動量測輸入訊號功率，並根據目標 SNR 計算雜訊功率。

* 功能：在輸入訊號 `x` 上加入符合指定 SNR 的 **白雜訊 (Gaussian Noise)**。

* 如果不加 `'measured'`，`awgn` 預設訊號功率為 0 dBW，可能導致錯誤的 SNR。

---

### 🔹 為什麼 `Np = Eb*10^(-Eb_N0/10)`？（公式推導與單位說明）

在通訊系統中， \$E\_b/N\_0\$  常以 dB 表示：

$$
\frac{E_b}{N_0} = 10^{\tfrac{(E_b/N_0)_{dB}}{10}}
$$


將其轉換為線性比值：

$$
\frac{E_b}{N_0} = 10^{\tfrac{(E_b/N_0)_{dB}}{10}}
$$

因此：

$$
N_0 = \frac{E_b}{E_b/N_0} = \frac{E_b}{10^{\tfrac{(E_b/N_0)_{dB}}{10}}}
$$

再整理：

$$
N_p = E_b \cdot 10^{- (E_b/N_0)_{dB} / 10}
$$

### 🔹 單位說明

* \$E\_b\$：單 bit 能量，單位為焦耳/bit (J/bit)。在程式中用 `mean(abs(x_n).^2)` 近似，因為 BPSK 下符號能量 = bit 能量。
* \$N\_0\$：雜訊功率譜密度，單位 W/Hz (瓦特每赫茲)。
* \$N\_p\$：對應的雜訊功率，單位 W (瓦特)。

👉 結論：`Np = Eb*10^(-Eb_N0/10)` 就是依據 **dB 定義公式**，將 \$E\_b/N\_0\$ (dB) 換算成雜訊功率，並搭配 MATLAB `awgn` 函數來產生對應的高斯白雜訊。

---

### 🔹 \$E\_b/N\_0\$ 與 SNR 的比較

* **\$E\_b/N\_0\$ (bit energy to noise spectral density ratio)**

  * 單一 bit 的能量與單位頻寬雜訊功率密度的比值。
  * 偏重於理論分析，常用於比較不同調變方式的效能。

* **SNR (Signal-to-Noise Ratio)**

  * 接收訊號平均功率與雜訊平均功率的比值。
  * 偏重於實際量測，與系統頻寬及符號率有關。


### 🔹 兩者的關係

一般公式：

$$
SNR = \frac{E_b}{N_0} \cdot \log_2(M)
$$

其中 \$M\$ 為調變星座點數。

* **BPSK (\$M=2\$)**：\$\log\_2(2) = 1\$ → \$SNR = E\_b/N\_0\$。
* **QPSK (\$M=4\$)**：\$SNR = 2 \cdot (E\_b/N\_0)\$。
* **16-QAM (\$M=16\$)**：\$SNR = 4 \cdot (E\_b/N\_0)\$。

✅ 總結：

* 在 **BPSK** 下，\$E\_b/N\_0\$ 與 SNR 幾乎相同。
* 在 **多進制調變 (QPSK, QAM)** 下，SNR 比 \$E\_b/N\_0\$ 高，差距為 \$\log\_2(M)\$ 倍。

---
