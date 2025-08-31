## (II)
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
## AWGN
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

### 🔹 程式逐行解釋

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

```matlab
% theoretical BER
BPSK_Pb(Eb_N0/2+1) = 0.5*erfc(sqrt(10^(Eb_N0/10)));
```

* BPSK 在 AWGN 通道下的理論 BER 公式：

  $P_b = Q\left(\sqrt{2E_b/N_0}\right) = \tfrac{1}{2} \, \text{erfc}\left(\sqrt{E_b/N_0}\right)$

* `10^(Eb_N0/10)`：把 dB 轉成線性。

* `erfc()`：互補誤差函數，用來計算 Q-function。

* 存到 `BPSK_Pb`，後面畫理論曲線。

---

```matlab
Np = Eb * 10^(-(Eb_N0/10));
```

* 根據公式 \$ N\_0 = E\_b / (E\_b/N\_0)\$，計算 **雜訊功率**。
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



