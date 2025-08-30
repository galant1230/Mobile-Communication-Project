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

* ** \$M=64\$ ** → 表示每個 OFDM 符號包含 **64 個子載波 (subcarriers)**。
* ** \$D\_k\$ ** → 第 \$k\$ 個子載波上的符號。

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

### 🔹 為什麼一個 Symbol 是 64 個子載波？

因為我們使用 \$M=64\$ 點 IFFT，頻域有 \$64\$ 個「bin」位置。

* 每個 bin 代表一個子載波 (subcarrier)。
* 把資料 \$D\_k\$ 放進去後，做 IFFT → 同時轉成一段時域波形。

所以：

* 一個 OFDM 符號就是 64 個子載波「同時傳送」的組合波。
* 這 64 個子載波中，只有 **52 個載波承載資料** (sc=52)，其餘 12 個是 DC + Guard band。

---

### 🔹 為什麼 BPSK 調變後是在頻域？

* 在傳統單載波系統中，BPSK 調變會直接變成時域波形 (例如 ±1 的方波)。
* **但在 OFDM 中，我們先把 BPSK 調變後的符號放到頻域的子載波位置 (\$D\_k\$)。**
* 接著利用 IFFT 把整個頻域訊號轉換成時域訊號 → 這才是實際要發送的 OFDM 符號。

所以流程是：

1. **比特 → 符號 (BPSK 調變)**：把 bits 映射成 ±1 (這時還只是「符號」，還沒分配頻率)。
2. **符號 → 子載波 (\$D\_k\$)**：把這些 ±1 塞進頻域的子載波 bin。
3. **頻域 → 時域 (IFFT)**：把整個頻域訊號轉成時域訊號，這就是要送出去的 OFDM symbol。

---

✅ 總結：

* 一個 **OFDM Symbol** 是一整段時域訊號，長度 \$M=64\$。
* 它由 **64 個子載波**組成，其中只有 52 個承載資料。
* BPSK 調變後的 ±1 先放在頻域子載波位置，再經過 IFFT 轉成時域，才會送到通道。


---

### 1.2. 資料子載波數 $sc = 52$

在 $M=64$ 的子載波裡，只有 **52 個承載資料**，其餘 **12 個作為保護子載波 (Guard bands) 與 DC**。  

### 1.3. 每個 OFDM Symbol 的 Bits 數 $ofdm\_bit = 52$

因為採用 **BPSK 調變**，每個子載波承載 **1 bit**：

$$
D_k \in \{-1, +1\}
$$

所以：  

$$
\text{每個 OFDM symbol bits 數} = sc \times 1 = 52
$$

### 1.4. 符號數 $Nsymbol = 100$

模擬傳送 $100$ 個 OFDM 符號：  

$$
N_{bits} = ofdm\_bit \times Nsymbol = 52 \times 100 = 5200
$$

### 1.5. 循環字首長度 $cp\_len = 16$

在時域 OFDM 符號前，複製最後 $16$ 個樣本加到最前面：  

$$
x_{cp}[n] =
\begin{cases}
x[n+M-cp\_len], & 0 \leq n < cp\_len \\
x[n-cp\_len], & cp\_len \leq n < M+cp\_len
\end{cases}
$$

功能：將通道的「**線性卷積**」轉為「**循環卷積**」，確保 FFT 還原正確。  

接收端移除 CP 後，頻域模型成立：  

$$
Y_k = H_k D_k + N_k
$$

---

