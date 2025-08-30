# 📡 行動通訊課程專案 #2

## 題目說明
在這個專案中，我們要透過模擬實現一個 **OFDM 收發機**。  
為了簡化模擬，DAC 與 ADC 被省略，使用的是 **離散時間基頻訊號模型**。  
假設已經達成完美的 **時間與頻率同步**，並且 **通道脈衝響應已知**。  

依據上圖流程，資訊位元會先被調變，再輸入至 $M$ 點 IDFT (IFFT)。  
加入循環字首 (Cyclic Prefix, CP) 之後，訊號 $x[n]$ 被送出。  

## 系統參數
- FFT 點數： $M = 64$  
- 資料子載波數： $52$  
- 循環字首長度： $16$ samples  
- 調變方式： BPSK  
- 子載波索引：  

$$
[7, 8, 9, 10, 11, 12, 13, 14,  
15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,  
27, 28, 29, 30, 31, 32, 34, 35, 36, 37, 38, 39,  
40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,  
52, 53, 54, 55, 56, 57, 58, 59]
$$

其中，  

$$
D_k \in \{-1, +1\}
$$  

為透過 BPSK 調變後的符號。  

## 接收訊號模型
接收訊號可表示為：  

$$
z[n] = x[n] * h[n] + n[n]
$$

其中：  
- $x[n]$ ：傳送訊號  
- $h[n]$ ：通道響應  
- $n[n]$ ：加性白雜訊 (AWGN)  

---
## (I) 實驗一：星座圖繪製

已知通道脈衝響應：  

$$
\tilde{h}[n] = [0.5 - 0.5j,\; 0,\; 0.15 + 0.12j,\; 0,\; 0,\; -0.1 + 0.05j]
$$  

並且正規化為：  

$$
h[n] = \frac{\tilde{h}[n]}{\text{norm}(\tilde{h}[n])}
$$  

### 實驗要求
- 模擬 **兩個 OFDM 符號週期**，並繪製輸出星座圖 (scatter plot)。  
- 在 **$M$ 點 DFT (FFT) 輸出**處觀察 $Y_k$。  
- 由於有 $2 \times 52 = 104$ 個子載波點，因此星座圖應包含 104 個點。  
- 雜訊條件設定為：  

$$
\frac{E_b}{N_0} = 20 \ \text{dB}
$$  

其中 $\frac{E_b}{N_0}$ 是在 FFT 輸出端量測的。  

### 等化器測試
- 請繪製通過 **頻域等化器 (frequency domain equalizer)** 後的星座圖。  
- 等化器採用 **零迫零 (Zero-Forcing) 準則**：  

$$
\hat{D}_k = \frac{Y_k}{H_k}
$$  

- 注意：需要將時域通道脈衝響應 $h[n]$ 轉換為頻域響應 $H_k$ (透過 FFT)。  
---
## Matlab ver
```matlab
%% ==================== OFDM System Simulation (BPSK) ==================== %%
clear all; close all; clc;

%% -------------------- 系統參數設定 -------------------- %%
M       = 64;                % FFT 點數 (Number of FFT points)
sc      = 52;                % 資料子載波數 (Number of data sub-carriers)
ofdmbit = 52;                % 每個 OFDM symbol 的 bits 數 (這裡用 BPSK，所以 bits = sub-carrier 數)
Nsymbol = 100;               % 模擬 OFDM 符號數 (Number of OFDM symbols)

%% -------------------- BPSK 調變 -------------------- %%
Data   = randi([0,1], ofdmbit*Nsymbol, 1);   % 產生隨機位元資料 (0/1)
dk     = 2*Data - 1;                         % BPSK 調變: 0→-1, 1→+1
dk_sym = reshape(dk, sc, Nsymbol);           % 串列轉並列 (Serial → Parallel)

%% -------------------- 子載波配置 -------------------- %%
% OFDM 有 64 個子載波，其中：
%   - 使用 52 個做資料傳送
%   - 中間 DC (第33個) 不用
%   - 兩側 guard band 也不用
Dk = zeros(M, Nsymbol);
Dk(7:32, :)  = dk_sym(1:26, :);   % 左半部 26 個資料子載波
Dk(34:59, :) = dk_sym(27:52, :);  % 右半部 26 個資料子載波

%% -------------------- IFFT (轉到時域) -------------------- %%
IDFT = [];
for i = 1:Nsymbol
    IDFT = [IDFT ifft(Dk(:, i), M)];
end

%% -------------------- 加入循環字首 (Cyclic Prefix, CP) -------------------- %%
D_cp = zeros(M+16, Nsymbol);
D_cp(1:16, :)   = IDFT(49:64, :);   % 複製最後 16 點到前面
D_cp(17:80, :)  = IDFT(1:64, :);    % 原本的 64 點

%% -------------------- 通道模型 (多路徑衰落) -------------------- %%
% 定義一個多徑通道 impulse response
h = [0.5-0.5j  0  0.15+0.12j  0  0  -0.1+0.05j];

%% -------------------- 串列化並通過通道 -------------------- %%
x_n = reshape(D_cp, 1, (M+16)*Nsymbol);   % 並列轉串列
h_n = h ./ norm(h);                       % 通道正規化
xh  = conv(x_n, h_n);                     % 傳輸: x * h

%% -------------------- 加入 AWGN 雜訊 -------------------- %%
Eb_N0 = 20;                              % Eb/No (dB)
Eb    = mean(abs(x_n).^2);               % 計算平均 bit 能量
Np    = Eb * 10^(-(Eb_N0/10));           % 噪聲功率
znn   = xh(1:end-5);                     % 調整長度 (扣掉通道記憶)
z_n   = awgn(znn, 20, 'measured');       % 加入 AWGN

%% -------------------- 移除循環字首 (CP Removal) -------------------- %%
receiver = reshape(z_n, 80, Nsymbol);
rn       = receiver(17:end, :);   % 去掉前 16 點，只留 64 點

%% -------------------- FFT (轉回頻域) -------------------- %%
DFT = [];
for i = 1:Nsymbol
    DFT = [DFT fft(rn(:, i), M)];
end

%% -------------------- 移除虛零子載波 (刪掉 DC 與 Guard) -------------------- %%
Yk64 = [DFT(7:32, :) ; DFT(34:59, :)];
Yk52 = [Yk64(:,1).' Yk64(:,2).'].';   % 串列化
scatterplot(Yk52);
title('Constellation (Before Equalization)');

%% -------------------- Equalizer (通道補償) -------------------- %%
Hn = fft(h, M);   % 通道頻率響應
Dk = [];
for i = 1:Nsymbol
    Dk = [Dk DFT(:,i) ./ Hn.'];   % Z(k)/H(k)
end

%% -------------------- 還原資料子載波 -------------------- %%
Dkreceiver = [Dk(7:32, :) ; Dk(34:59, :)];   % 去掉不用的載波
Dkr        = reshape(Dkreceiver, 1, 52*Nsymbol);

scatterplot(Dkr(1, 1:104));
title('Constellation (After Equalization)');

```
---

## (I) Result
<img src="2.jpg" alt="2" width="50%">
<img src="3.jpg" alt="3" width="50%">

---
## (II) 實驗二：BER 曲線繪製 (理想通道)

### 實驗要求
- 繪製 **位元錯誤率 (BER) 曲線**。  
- 誤碼率統計需在以下訊雜比條件下進行：  

$$
\frac{E_b}{N_0} = 0, \; 2, \; 4, \; 6, \; 8 \ \text{dB}
$$  

- 通道模型為 **理想通道**：  

$$
h[n] = [1]
$$  

---

## (III) 實驗三：BER 曲線繪製 (多路徑通道)

### 實驗要求
- 繪製 **位元錯誤率 (BER) 曲線**。  
- 誤碼率統計需在以下訊雜比條件下進行：  

$$
\frac{E_b}{N_0} = 0, \; 2, \; 4, \; 6, \; 8 \ \text{dB}
$$  

- 通道模型為 **多路徑通道**，其脈衝響應為：  

$$
h[n] = \frac{\tilde{h}[n]}{\text{norm}(\tilde{h}[n])}
$$  

在這種情況下，因為通道存在衰落與多路徑效應，模擬得到的 BER 曲線會與理論值有明顯差異。


