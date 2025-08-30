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

### 1.3 每個 OFDM Symbol 的 Bits 數

由於採用 **BPSK 調變**，每個子載波可以承載 1 bit：

* BPSK 映射：

  * bit = 0 → \$D\_k = -1\$
  * bit = 1 → \$D\_k = +1\$

因此：

$$
\text{每個 OFDM symbol 的 bits 數} = sc \times 1 = 52
$$

換句話說，在 \$M=64\$ 個子載波中，**只有 52 個子載波承載資料**，因此一個 OFDM symbol 最多傳送 52 個 bits。

---

### 1.4 符號數 \$N\_{symbol} = 100\$

假設模擬傳送 **100 個 OFDM 符號**，則總傳送 bit 數為：

$$
N_{bits} = ofdm\_bit \times N_{symbol} = 52 \times 100 = 5200
$$

這表示：在模擬中，系統會隨機產生並傳送 **5200 個 bits**，以進行性能分析 (例如 BER)。

---

### 1.5 循環字首長度 \$cp\_len = 16\$

在 OFDM 中，為了避免通道造成的延遲擾動 (ISI)，會在每個時域 OFDM 符號的前端加上 **循環字首 (Cyclic Prefix, CP)**。方法是：

把最後 \$cp\_len=16\$ 個樣本複製到最前面，形成新的序列：

$$
x_{cp}[n] =
\begin{cases}
  x[n + M - cp\_len], & 0 \leq n < cp\_len \\
  x[n - cp\_len], & cp\_len \leq n < M+cp\_len
\end{cases}
$$

### 🔹 循環字首的功能

1. **避免符號間干擾 (ISI)**：多路徑通道會造成延遲疊加，CP 的存在確保延遲不會影響到下一個符號。
2. **將線性卷積轉為循環卷積**：

   * 在無 CP 時，通道效應為 **線性卷積**：
     $y[n] = (x * h)[n] + w[n]$
   * 加上 CP 後，卷積可以近似為 **循環卷積**：
     $y[n] = (x \circledast h)[n] + w[n]$

   這樣在頻域就能得到簡單關係：

   $Y_k = H_k \cdot D_k + N_k$

   * \$H\_k\$：通道在第 \$k\$ 個子載波的頻率響應。
   * \$D\_k\$：發送符號 (BPSK)。
   * \$N\_k\$：雜訊。

### 🔹 優點

有了這個頻域模型：

* 通道等化變得非常簡單 (只需做 **除法 \$D\_k = Y\_k / H\_k\$**)。
* 系統能在多路徑環境下維持正確的頻域解調。

---

✅ 總結：

* **1.3**：BPSK → 每個子載波 1 bit，一個 OFDM symbol = 52 bits。
* **1.4**：模擬 100 個符號 → 總共 5200 bits。
* **1.5**：CP (16 個樣本) 讓 **頻域關係 \$Y\_k = H\_k D\_k + N\_k\$ 成立**，等化更簡單。

---

### 1.3 每個 OFDM Symbol 的 Bits 數（BPSK）

因為採用 BPSK 調變，每個子載波承載 1 bit，對應到頻域符號集合：D\_k ∈ {-1, +1}。因此在本模擬（52 個資料子載波）下：

* 每個 OFDM symbol 的 bits = sc × 1 = 52

一般化：若改成 M\_mod-QAM，則每個 OFDM symbol 的 bits = sc × log2(M\_mod)（例如 QPSK=2、16-QAM=4）。

和你的程式對應：

* sc = 52，ofdm\_bit = 52（BPSK → 每個子載波 1 bit）。
* dk = 2\*Data - 1 進行 BPSK 映射（0→-1、1→+1）。
* Dk(7:32,:) 與 Dk(34:59,:) 將 52 個資料載波填入頻域（跳過 DC=索引 33，並在兩端留 Guard bands）。

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

### 1.5 循環字首長度 cp\_len = 16

在時域 OFDM 符號前，複製最後 cp\_len 個樣本加到最前面：

x\_cp\[n] =

* x\[n + M - cp\_len], 當 0 ≤ n < cp\_len；
* x\[n - cp\_len], 當 cp\_len ≤ n < M + cp\_len。

目的與條件：

1. 避免符號間干擾（ISI）：若通道最長延遲擴散（最大有效路徑延遲）為 L-1，則只要 cp\_len ≥ L-1，即可確保前一個符號的尾巴不會干擾到目前符號。
2. 把線性卷積變成循環卷積：去掉 CP 後的 M 點區段滿足循環卷積，取 M 點 DFT 後得到子載波分離關係：Y\[k] = H\[k] \* D\[k] + W\[k]。其中 H\[k] 為通道的頻率響應（對 h\[l] 做 M 點 DFT）。這使得每個子載波可獨立等化。

和你的程式對應：

* 通道 h = \[0.5-0.5j, 0, 0.15+0.12j, 0, 0, -0.1+0.05j]，有效 taps 在 l = 0, 2, 5，故最大延遲 = 5。→ cp\_len = 16 ≥ 5，因此滿足條件，ISI 可被 CP 吸收，頻域模型成立。
* 頻域等化：Hn = fft(h, M)，Dk\_eq = DFT ./ Hn.' 為 Zero-Forcing (ZF) 等化，即 D\_hat\[k] = Y\[k] / H\[k]。小提醒：當 |H\[k]| 很小時，ZF 會放大雜訊；若要抑制雜訊，可改用 MMSE 等化：D\_hat\[k] = (H\*\[k] / (|H\[k]|^2 + sigma^2)) \* Y\[k]。

索引示意（64 點 FFT）：

* 左側 Guard：索引 1–6
* 資料：7–32（26 條）
* DC：33（不載資料）
* 資料：34–59（26 條）
* 右側 Guard：60–64
  → 共 52 條資料子載波 + 1 條 DC + 11 條 Guard bands。

---

