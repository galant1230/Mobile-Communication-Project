# 1. 🌊 Jakes Model 與 Gaussian Rayleigh 的對應與實作筆記

---

### 🧠 主要目的：

用 **Jakes model** 模擬時間變化的 **Rayleigh fading**（即通道係數 $h(t)$ 隨時間變化的隨機過程）。

### Jakes Model 簡介

* **用途**：用於模擬移動通訊中的 **時間相關 Rayleigh 衰落**。

* **核心想法**：

  * 把接收信號看成許多不同角度的平面波疊加。
  * 每個平面波因為移動造成 **多普勒效應**，頻率略有偏移。
  * 疊加後形成一個隨機但有時間相關性的信號。

* **數學模型**：

$$
T(t) = E_0 \sum_{m=1}^N C_m e^{j(\cos \alpha_m \, \omega_d t + \phi_m)}
$$

  * \$N\$：平面波數目（通常取 30 即可近似）
  * \$C\_m\$：權重 (\$C\_m^2 = 1/N\$)
  * \$\alpha\_m\$：到達角度，假設均勻分佈
  * \$\omega\_d = 2\pi f\_d\$：最大多普勒頻率
  * \$\phi\_m\$：隨機初始相位
    
---

## 🧪 為什麼 Jakes 模型的 Envelope 是 Rayleigh？

> 這個複數總和公式長得跟 Rayleigh 分布好像不太一樣，為什麼取絕對值之後就會是 Rayleigh？

讓我們一步步解析：

### 📌 Step 1：Jakes 模型公式

$$
T(t) = E_0 \sum_{m=1}^{N} C_m e^{j(\cos \alpha_m \, \omega_d t + \phi_m)}
$$

我們關心的是其 **envelope**：

$$
|T(t)| = \left| \sum_{m=1}^N C_m e^{j(\cos \alpha_m \omega_d t + \phi_m)} \right|
$$

---

### 📌 Step 2：轉換為實部與虛部的總和

利用歐拉公式：

$$
e^{j\theta} = \cos(\theta) + j \sin(\theta)
$$

所以總和變成：

$$
T(t) = \sum C_m \cos(\psi_m(t)) + j \sum C_m \sin(\psi_m(t))
$$

其中 $\psi_m(t) = \cos \alpha_m \, \omega_d t + \phi_m$

---

### 📌 Step 3：中心極限定理（CLT）應用

若：

* $\phi_m \sim \text{Uniform}(0,2\pi)$
* $\alpha_m$ 多方向分布 → 對應不同 Doppler
* $C_m = 1/\sqrt{N}$（或常數）

則：

* $\text{Re}(T(t)) \sim \mathcal{N}(0,\sigma^2)$
* $\text{Im}(T(t)) \sim \mathcal{N}(0,\sigma^2)$

因為是很多獨立小分量的總和，根據中心極限定理會趨近 Gaussian！

---

### 📌 Step 4：複數 Gaussian 取絕對值會變成 Rayleigh

$$
|T(t)| = \sqrt{\text{Re}(T(t))^2 + \text{Im}(T(t))^2} \sim \text{Rayleigh}(\sigma)
$$

這正是 Rayleigh 分布的定義：兩個獨立高斯隨機變數（均值 0）組合成的長度。
   
---

## 📌 原始 MATLAB 程式碼（Jakes model + Gaussian Rayleigh 對比）

```matlab
clc;clear;close all;
%------------------先給定參數--------------------%
T=1;                        % 錄時間長 T=1，若設太多，例如100則電腦無法負荷停量
delta_t=0.000001;           % 以下為總時間所佔的微度
N=30;
f=3.5*10^9;
c=3*10^8;
Eo=sqrt(2);
v=(60000)/3600;             % 先換成(m/s)較為常用的單位
landa=c/fc;
t=0:delta_t:T;              % 以下為所使用的時間範圍(1/m)
n=0:1:N-1;
m=0:1:N-1;
cn=msqrt(1/N);              % 以下為baseband representation的參數表示
w_m=(2*pi*v)/landa;
rfa_m=(2*pi*m)/N;

%------------------計算Jakes's model的h(t)--------------------%
for i=1:N                     % 用for來重Tisigama運算
    phi_m=2*pi*rand(1,1);     % 代表的是uniform distribution
    h(i,:)=cn*exp(j*(cos(rfa_m(i))*w_m*n'+phi_m));
end
h_t=Eo*sum(h);

%------------------計算Jakes's rayleigh fading pdf--------------------%
A=abs(h_t);                   % 取絕對值
[a_h,y]=hist(A,h_y);         % hist指令將整體資料依長度大小分為組，各組件中之數目作為繪圖之根據。
p_a_h=a_h./(sum(a_h)*(0.1)); % 每個bin寬是0.1 → 統一一個相相乘 [2 3].*2 [2 3 4]  本行是要算密度

%------------------再算 Gaussian rayleigh fading pdf--------------------%
var=1;
c=sqrt(var/2)*(randn(1,1000000)+j*randn(1,1000000)); % 產生一般的complex Gassion
a_c=abs(c);
y_c=[0:0.1:10];
[a_c_y,c_y]=hist(a_c,y_c);
p_a_c=a_c./(sum(a_c)*(0.1));

%------------------此處為計算平均值與變異數--------------------%
mean_c=mean(c);mean_h_t=mean(h_t);
var_c=std(c)^2;var_h_t=std(h_t).^2;
mean_A_c=mean(A);mean_A_h=mean(A_h);
var_A_c=std(A_c)^2;var_A_h=std(A_h).^2;

%------------------繪圖顯示--------------------%
figure(1);
plot(y,p_a_h,'c',y_c,p_a_c,'r--');     % 將Gausison Rayleigh使用紅色虛線表示
xlabel('h(t)/var');ylabel('P(h(t)/var)');
title('Normalized PDF of Gaussian rayleigh and Jakes rayleigh');
legend('Jakes rayleigh','Gaussian rayleigh')
axis([0 4.5 0 1])

figure(2);
hist(A,h_y);title('Histogram of Jakes rayleigh');
xlabel('frequency');
ylabel('Number of Hits');
axis([0 4.5 0 10^5])

figure(3);
hist(A_c,y_c);title('Histogram of Gaussian rayleigh');
xlabel('|y(t)|');ylabel('frequency');
axis([-0.5 4.5 0 10^5]);
```

---

### 📌 Part 1 : 參數

```matlab
clc; clear; close all;

%------------------先給定參數--------------------%
T = 1;                         % 總模擬時間 T（秒）
delta_t = 0.000001;            % 時間間隔 Δt（秒） → 影響時間解析度
N = 30;                        % 散射路徑數 N
f = 3.5*10^9;                  % 頻率 f（Hz）
c = 3*10^8;                    % 光速（m/s）
Eo = sqrt(2);                  % 正規化能量（後面乘在總和上）
v = 60000/3600;                % 移動速度（km/hr → m/s）

landa = c / f;                % 波長 λ = c / f

% 時間軸定義，從 0 到 T，以 delta_t 為間隔
% 會產生一個時間序列，用來觀察 fading 過程

t = 0:delta_t:T;              

% 基底向量 index（for multiple path sum）
n = 0:1:N-1;                  % 對應時間序列樣本 index
m = 0:1:N-1;                  % 用來產生不同方向與相位的 index

% baseband representation 的係數正規化常數
cn = sqrt(1/N);               % 每個散射路徑的功率平均化（均分能量）

% 計算最大 Doppler 頻率：w_m = 2πf_D
w_m = (2*pi*v) / landa;       % Doppler spread（與移動速度、波長有關）

% 將每條路徑的角度轉換成參數化角度 rfa_m
rfa_m = (2*pi*m) / N;         % 均勻分佈於 [0, 2π) 範圍內
```

---

### ✅ 對應到 Jakes 模型公式：

Jakes 模型完整表達式：

$$
T(t) = E_0 \sum_{m=1}^N C_m e^{j(\cos\alpha_m \omega_d t + \phi_m)}
$$

| MATLAB 變數  | 對應參數       | 說明                   |
| ---------- | ---------- | -------------------- |
| `Eo`       | $E_0$      | 輸出能量放大因子             |
| `cn`       | $C_m$      | 每條路徑的權重（通常平均）        |
| `rfa_m(i)` | $\alpha_m$ | 每條路徑的角度（經由 index 均分） |
| `w_m`      | $\omega_d$ | 最大 Doppler 頻率        |
| `phi_m`    | $\phi_m$   | 初始相位（隨機）             |
| `t`, `n`   | $t$        | 時間序列                 |

這段 code 就是完整對應上面公式中每一項的初始化步驟。

---
### 📌 Part 2 : 計算 Jakes 的 h(t) 

```markdown
%------------------計算 Jakes 的 h(t) --------------------%
for i = 1:N
    phi_m = 2*pi*rand(1,1);                       % 初始隨機相位 φ_m ∈ [0, 2π)
    h(i,:) = cn * exp(j*(cos(rfa_m(i))*w_m*n' + phi_m));
end
h_t = Eo * sum(h);                                % 將所有 N 條路徑相加後乘上能量因子
```

---

### ✅ 對應到 Jakes 模型公式：

回顧 Jakes 模型核心公式：

$$
T(t) = E_0 \sum_{m=1}^{N} C_m e^{j(\cos\alpha_m \, \omega_d t + \phi_m)}
$$

這段程式碼對應如下：

| MATLAB 程式碼                         | 對應公式項                                 | 說明                  |
| ---------------------------------- | ------------------------------------- | ------------------- |
| `phi_m = 2*pi*rand(1,1);`          | $\phi_m \sim \text{Uniform}(0, 2\pi)$ | 每條路徑隨機初始相位          |
| `cos(rfa_m(i)) * w_m * n' + phi_m` | $\cos\alpha_m \, \omega_d t + \phi_m$ | 多徑方向造成不同 Doppler 偏移 |
| `exp(j*(...))`                     | $e^{j(...)}$                          | 複數正弦波（旋轉向量）         |
| `h(i,:) = cn * exp(...)`           | $C_m e^{j(...)}$                      | 每條路徑帶有相同權重 C\_m     |
| `h_t = Eo * sum(h);`               | $E_0 \sum_{m=1}^{N} ...$              | 最後的通道係數是所有路徑加總      |

---

## Part 3 : Jakes model + Gaussian Rayleigh 對比

```matlab
%------------------計算 Jakes 的 Rayleigh fading PDF --------------------%
A = abs(h_t);                           % 對複數通道係數取絕對值 → 得到 fading 包絡
[a_h, y] = hist(A, h_y);               % 將 A 資料分成區間 (bin)，得到直方圖
p_a_h = a_h ./ (sum(a_h) * 0.1);       % 將直方圖正規化成機率密度函數（PDF），bin 寬為 0.1

%------------------計算 Gaussian Rayleigh fading PDF --------------------%
var = 1;                                % 設定變異數
c = sqrt(var/2) * (randn(1,1000000) + j*randn(1,1000000));
                                       % 產生複數 Gaussian 隨機變數，實虛部獨立

a_c = abs(c);                          % 取 magnitude 得到 Rayleigh 分布樣本

y_c = 0:0.1:10;                         % 分割範圍
[a_c_y, c_y] = hist(a_c, y_c);         % 建立直方圖
p_a_c = a_c_y ./ (sum(a_c_y) * 0.1);   % 正規化成機率密度函數（PDF）
```

---

### ✅ 對應說明：

#### Jakes 模型部分：

這裡從前面建構好的時間序列通道係數 `h_t`，
經由：

```matlab
A = abs(h_t)
```

取得其包絡 → **這就是 Rayleigh fading** 的樣本。

再使用 `hist()`：

* 對這些樣本作統計，計算落在每個 bin 區間的數量
* 除以總樣本數與 bin 寬，即可正規化為 PDF

---

#### Gaussian 模型部分：

這部分直接產生理論上的 Rayleigh 分布：

* 複數數列 `c = x + jy`，其中 x, y 為獨立高斯
* 取其 `abs(c)` 得到 Rayleigh 分布

這裡的 `1000000` 是樣本數，越大越接近理論曲線。

---

| 模型                | 建構方式                        | 分布類型              |
| ----------------- | --------------------------- | ----------------- |
| Jakes 模型          | 多徑複數正弦和，取 envelope          | Rayleigh-like（實測） |
| Gaussian Rayleigh | 實部虛部獨立 Gaussian，取 magnitude | Rayleigh（理論）      |

---
##  Part 4 : Figure
```matlab
%------------------此處為計算平均值與變異數--------------------%
mean_c = mean(c); % 複數 Gaussian 模型的均值（理論上應接近 0）
mean_h_t = mean(h_t); % Jakes 模型的均值（理論上應接近 0）


var_c = std(c)^2; % 複數 Gaussian 模型的變異數
var_h_t = std(h_t).^2; % Jakes 模型的變異數


mean_A_c = mean(a_c); % Gaussian Rayleigh 包絡的均值（理論上 = √(π/2)）
mean_A_h = mean(A); % Jakes Rayleigh 包絡的均值（模擬值）


var_A_c = std(a_c)^2; % Gaussian Rayleigh 包絡的變異數（理論上 = (4 - π)/2）
var_A_h = std(A).^2; % Jakes Rayleigh 包絡的變異數（模擬值）


%------------------繪圖顯示--------------------%
figure(1);
plot(y, p_a_h, 'c', y_c, p_a_c, 'r--'); % 將 Jakes Rayleigh 與 Gaussian Rayleigh 的 PDF 畫在同一張圖
xlabel('h(t)/var');
ylabel('P(h(t)/var)');
title('Normalized PDF of Gaussian rayleigh and Jakes rayleigh');
legend('Jakes rayleigh', 'Gaussian rayleigh');
axis([0 4.5 0 1]);


figure(2);
hist(A, h_y); % Jakes 模型的 Rayleigh 包絡直方圖
xlabel('frequency');
ylabel('Number of Hits');
title('Histogram of Jakes rayleigh');
axis([0 4.5 0 1e5]);


figure(3);
hist(a_c, y_c); % Gaussian 模型的 Rayleigh 包絡直方圖
xlabel('|y(t)|');
ylabel('frequency');
title('Histogram of Gaussian rayleigh');
axis([-0.5 4.5 0 1e5]);
```

---
## 🌊 Jakes Model 與 Gaussian Rayleigh 的對應與實作筆記

---

### 🧠 主要目的：

用 **Jakes model** 模擬時間變化的 **Rayleigh fading**（即通道係數 $h(t)$ 隨時間變化的隨機過程）。

---

## 📌 補充模型：另一種 Rayleigh Fading 模擬法（圖三）

以下是另一份 MATLAB 實作，使用幾何方式模擬散射波的相位差，形成 Rayleigh fading：

```matlab
clc; clear; close all;

%------------------ 先給定參數 ------------------%
n = 30;                           % number of scatterers
c = sqrt(1/n);                    % 每條路徑平均能量分配
Eo = 1;                           % 波的總能量
v = (60000)/3600;                 % 移動速度 (m/s)
fc = 3.5*10^9;                    % 頻率 (Hz)
landa = 3*10^8 / fc;              % 波長 λ

d = 1:n;
wn = 2 * pi * v / landa;          % 最大 Doppler 角頻率
phi = 2 * pi * (d - 1) / n;       % 均勻分佈的入射角

t = 0:0.000001:1;                 % 時間序列

theta = (2 * pi) * rand(1, n);    % 每條路徑的隨機相位 θ_k

xc = Eo * c * cos(wn * t);        % 初始化：實部
xs = zeros(1, length(xc));        % 初始化：虛部

%------------------ 計算信號總和 ------------------%
for k = 1:n
    xc = xc + c * cos(cos(phi(k)) * wn * t + theta(k));   % 加總每個路徑的實部
    xs = xs + c * sin(cos(phi(k)) * wn * t + theta(k));   % 加總每個路徑的虛部
end

ag = abs(xc + sqrt(-1) * xs);     % 得到 Rayleigh 包絡（取 magnitude）

%------------------ 繪圖 ------------------%
hist(ag, 300);
ylabel('Number of Hits');
xlabel('Magnitude');

figure(2);
semilogy(0.000001*(1:100000), ag(1:100000));
xlabel('Time');
ylabel('Magnitude (dB)');
```

---

### ✅ 對應數學公式：

這裡對應的基本模型公式為：

$$
T(t) = \sum_{k=1}^{N} A_k e^{j(\omega_d \cos(\phi_k) t + \theta_k)}
$$

拆成實部與虛部：

$$
\text{Re}[T(t)] = \sum A_k \cos(\omega_d \cos\phi_k t + \theta_k) \\
\text{Im}[T(t)] = \sum A_k \sin(\omega_d \cos\phi_k t + \theta_k)
$$

包絡為：

$$
|T(t)| = \sqrt{\text{Re}^2 + \text{Im}^2} \sim \text{Rayleigh}
$$

---

### 📌 對應說明表格：

| MATLAB 變數        | 數學符號       | 說明                   |   |                      |
| ---------------- | ---------- | -------------------- | - | -------------------- |
| `n`              | $N$        | 散射路徑數                |   |                      |
| `phi(k)`         | $\phi_k$   | 第 k 條路徑的入射角          |   |                      |
| `theta(k)`       | $\theta_k$ | 初始相位（隨機）             |   |                      |
| `wn`             | $\omega_d$ | 最大 Doppler shift     |   |                      |
| `cos(phi(k))*wn` | $\omega_k$ | 每條路徑的 Doppler offset |   |                      |
| `ag = abs(...)`  | (          | T(t)                 | ) | 總包絡 → 應為 Rayleigh 分布 |

---

### 🧪 與前面 Jakes Model 相比：

| 項目     | Jakes Model             | 幾何散射模型（本段）         |
| ------ | ----------------------- | ------------------ |
| 路徑相位建模 | 統一相位與方向參數               | 使用角度 `phi(k)` 明確建模 |
| 實虛部構成  | 直接複數旋轉向量累加              | 拆實部與虛部，顯式加總        |
| 調變方式   | $\cos(\alpha) \omega t$ | 同上                 |
| 結果包絡   | Rayleigh 分布（模擬驗證）       | Rayleigh 分布（模擬驗證）  |

---

## 📌 3. 模型：Rician Fading Channel 模擬（有 LOS 成分）

以下 MATLAB 程式碼為模擬具有直視路徑 (Line-of-Sight, LOS) 的 **Rician Fading Channel**：

```matlab
clc; clear; close all;

%------------------先給定參數--------------------%
n = 30;                              % 散射路徑數
c = sqrt(1/n);                       % 每條路徑平均能量分配
km = 2;                              % K-factor（線性）
kf = 10*log10(km);                   % K-factor（dB）
Eo = 1;                              % 總能量
A = sqrt(2 * kf * Eo^2);             % LOS 分量的振幅
v = 60000 / 3600;                    % 移動速度 (m/s)
fc = 3.5*10^9;                       % 頻率 (Hz)
landa = 3e8 / fc;                    % 波長

d = 1:n;
wn = 2 * pi * v / landa;            % 最大 Doppler shift
phi = 2 * pi * (d - 1) / n;          % 均勻分佈角度

t = 0:0.000001:1;
theta = 2*pi*rand(1, n);            % 每條路徑隨機相位

xc = Eo * c * cos(wn * t);          % 初始化實部
xs = zeros(1, length(xc));          % 初始化虛部
xm = zeros(1, length(xc));          % 初始化 LOS 分量

%------------------加總多徑分量--------------------%
for k = 1:n
    xc = xc + c * cos(cos(phi(k)) * wn * t + theta(k));  % 散射實部
    xs = xs + c * sin(cos(phi(k)) * wn * t + theta(k));  % 散射虛部
end

%------------------加上直視 LOS 成分--------------------%
xm = A * cos(wn * t) + 1i * A * sin(wn * t);  % 產生 A * exp(jωt) 的複數直視分量

ag = abs(xm + xc + 1i * xs);       % 總和後取 magnitude

%------------------繪圖顯示--------------------%
hist(ag, 300);
ylabel('Number of Hits');
xlabel('Magnitude');

figure(2);
semilogy(0.000001*(1:100000), ag(1:100000));
xlabel('Time');
ylabel('Magnitude');
```

---

### ✅ 對應公式說明

Rician Fading 的通道模型可以寫成：

$$
T(t) = A e^{j\omega t} + \sum_{k=1}^{N} a_k e^{j(\omega_k t + \phi_k)}
$$

* 第一項：直視 LOS 分量
* 第二項：N 條多徑散射分量

**取 magnitude** 後得到的是 **Rician 分布**。

---

### 📌 對應表格

| MATLAB 變數       | 對應公式參數                   | 說明                |   |              |
| --------------- | ------------------------ | ----------------- | - | ------------ |
| `A`             | LOS 振幅                   | 強度由 K-factor 決定   |   |              |
| `xm`            | LOS 分量 $A e^{j\omega t}$ | 使用 cos + j sin 組成 |   |              |
| `xc + j*xs`     | 多徑散射分量總和                 | N 條散射方向相加         |   |              |
| `ag = abs(...)` | (                        | T(t)              | ) | 得到 Rician 包絡 |

---

### 📊 與 Rayleigh 比較：

| 模型類型            | 是否有 LOS？ | 分布類型        | 控制參數            |
| --------------- | -------- | ----------- | --------------- |
| Rayleigh Fading | 無 LOS    | Rayleigh 分布 | 無               |
| Rician Fading   | 有 LOS    | Rician 分布   | K-factor（強 LOS） |

---
# 4. Ber Vs Ebno Analysis
```matlab
% -------------------- 先處理參數 ------------------------
clear all; close all;
n=32;                   % 資料數
Ts=1e-3;                % sampling time
T=2*1e-3;               % time
landa=0.15;             % wavelength (lambda)
v=60000/3600;           % velocity (m/s)
N=8;                    % number of scatterers

% -------------------- 隨機產生 binary source ------------------------
x=rand(1,n);            % 原始 bit stream
d=sign(x-0.5);          % 轉為 +1/-1 binary 資料

% -------------------- 產生 Jakes 模型 h(t) ------------------------
Cm=sqrt(1/N);
wn=2*pi*v/landa;        % 角都卜勒頻率
h_t=zeros(1,T/Ts);
for m=0:N-1
    am=2*pi*m/N;
    q=rand(1,1); % random 相位
    h_tmp=Cm*exp(j*(cos(am)*wn*(Ts:Ts:T)+q));
    h_t=h_t+h_tmp;
end
h=h_t./max(abs(h_t));   % normalize

% -------------------- 產生 Rayleigh 通道 (高斯分佈) ------------------------
ray_tmp=randn(1,T/Ts)+j*randn(1,T/Ts);
r=ray_tmp/max(abs(ray_tmp));

% -------------------- 模擬傳輸 ------------------------
xt=r.*d;                % Rayleigh 模型乘上資料
xh=h.*d;                % Jakes 模型乘上資料

% -------------------- 模擬通道 + AWGN ------------------------
for snr=1:6
    xh=awgn(xh,snr,'measured');
    xr=awgn(xt,snr,'measured');

    % receiver 端做 equalization
    yh=conj(h).*xh;           % 解調 (Jakes)
    yr=conj(r).*xr;           % 解調 (Rayleigh)

    % 判斷符號
    yhd=sign(real(conj(h).*xh));
    yrd=sign(real(conj(r).*xr));

    % 計算 BER
    [~,ratio_h]= symerr(d,yhd);
    [~,ratio_r]= symerr(d,yrd);

    ber_h(n)=ratio_h;
    ber_r(n)=ratio_r;
end

% -------------------- 畫圖 ------------------------
figure
subplot(2,1,1);
x_axe=1:6;
semilogy(x_axe,ber_h,'r.-');
title('BER (Jakes model)');
xlabel('SNR (dB)'); ylabel('BER');

subplot(2,1,2);
semilogy(x_axe,ber_r,'b.-');
title('BER (Rayleigh model)');
xlabel('SNR (dB)'); ylabel('BER');

% 比較圖
figure;
semilogy(x_axe,ber_h,'r-o',x_axe,ber_r,'b-o');
legend('jakes','rayleigh');
xlabel('SNR(dB)'); ylabel('BER');
grid on;
```

---

## ：BER vs Eb/No 模擬

### 📌 模型設定區段（對應理論公式）

```matlab
Ts=1e-3;         % 取樣時間 ∆t
T=2*1e-3;        % 總模擬時間 T
N=8;             % 散射體數目（路徑數）
v=60000/3600;    % 使用者移動速度 v (m/s)
fc=3.5*10^9;     % 頻率 fc
landa=3e8/fc;    % 波長 λ = c/f
```

### 🔹 對應理論：

**都卜勒頻率：**

$$
f_D = \frac{v}{c}f_c \quad \Rightarrow \quad \omega_D = 2\pi f_D = \frac{2\pi v}{\lambda}
$$

這邊的 `wn = 2*pi*v/landa` 對應的就是都卜勒角頻率

---

### 📌 通道產生 (Jakes 模型)

```matlab
for m=0:N-1
    am=2*pi*m/N;
    q=rand(1,1);   % 隨機相位
    h_tmp=Cm*exp(j*(cos(am)*wn*(Ts:Ts:T)+q));
    h_t=h_t+h_tmp;
end
```

**對應公式：**

Jakes 模型的 channel impulse response (無 LOS)：

$$
h(t) = \sum_{m=1}^{N} C_m e^{j(\omega_D \cos\alpha_m t + \phi_m)}
$$

其中：

* \$C\_m = \frac{1}{\sqrt{N}}\$
* \$\alpha\_m = \frac{2\pi m}{N}\$
* \$\phi\_m \sim \text{Unif}(0,2\pi)\$

---

### 📌 正規化處理

```matlab
h=h./max(abs(h_t));
```

📈 為了讓所有的 fading 都有統一平均力量

---

### 📌 建立高斯隨機通道 (Rayleigh 分佈)

```matlab
ray_tmp=randn(1,T/Ts)+j*randn(1,T/Ts);
r=ray_tmp/max(abs(ray_tmp));
```

**對應理論：**

Rayleigh fading 模型定義：

$$
h(t) = X(t) + jY(t),\quad X,Y \sim \mathcal{N}(0,\sigma^2)
$$

其 envelope \$|h(t)|\$ 會呈 Rayleigh 分佈

---

### 📌 BER 模擬

```matlab
d=sign(rand(1,n)-0.5);    % binary source: ±1
x=r.*d;                   % Rayleigh
x=h.*d;                   % Jakes
```

模擬輸入資料經通道處理

---

### 📌 加上噪聲 + 解讀

```matlab
xh=awgn(xt,snr,'measured');
```

對應理論公式：

$$
y(t) = h(t)x(t) + n(t)
$$

受 noise 影響

```matlab
yh=conj(h).*xh;
yhd=sign(real(yh));
```

對應二項法的 channel equalization

---

### 📌 错誤率計算

```matlab
[number_h,ratio_h] = symerr(d,yhd);
[number_r,ratio_r] = symerr(d,yrd);
```

BER 公式：

$$
BER = \frac{\text{number of error bits}}{\text{total transmitted bits}}
$$

---

### 📊 圖表呈現

```matlab
semilogy(x_axe,ber_h,'r-o'); % Jakes
semilogy(x_axe,ber_r,'b-o'); % Rayleigh
legend('jakes','rayleigh');
```

semilogy 用於 y-axis 為 log-scale

---

### ✅ 小結分析

| 項目            | 對應模型          | 影響說明                    |
| ------------- | ------------- | ----------------------- |
| `h_t`         | Jakes         | 模擬 multipath fading     |
| `ray_tmp`     | Rayleigh      | 理想 Gaussian 通道          |
| `x = h.*d`    | channel input | 通道影響輸入                  |
| `awgn()`      | noise         | 加入 AWGN                 |
| `conj(h).*xh` | equalization  | 使用 channel estimate 做解碼 |
| `symerr()`    | BER           | 算計誤碼率                   |

---
# 5. Ber Kfactor Bpsk·

```matlab
clc;clear; 
%----------------------------- 先給定參數 -----------------------------%
T=3;                            % 總模擬時間
delta_t=0.0001;                 % 每個時間點的間距
N=100;                          % 散射體數量（注意太大會造成效能問題）
fc=3.5*10^9;                    % 載波頻率
c=3*10^8;                       % 光速
Eo=sqrt(10);                    % 能量參數（根號10）
v=(60000)/3600;                 % 使用者速度 (m/s)
landa=c/fc;                     % 波長
cm=sqrt(1/N);                   % 每條路徑的增益（根號1/N）
w_n=(2*pi*v)/landa;             % 都卜勒角頻率
T=0:delta_t:T;                  % 時間向量
i=0:1:N-1;                      % 散射體索引
rfa_m=(2*pi*i/N);               % 每個路徑的散射角
%------------------------- 計算Jakes模型的h(t) -------------------------%
for i=1:N
    phi_m=2*pi*rand(1,1);                                  % 每條路徑一個隨機相位
    h1(i,:)=cm*exp(1j*(cos(rfa_m(i))*w_n*T+phi_m));       % 第一組Jakes波形
end

for i=1:N
    phi_m=2*pi*rand(1,1);                                  % 第二組Jakes波形
    h2(i,:)=cm*exp(1j*(cos(rfa_m(i))*w_n*T+phi_m));
end

h_t1=Eo*sum(h1);      % 第一組 Jakes 模型總和（強LOS情況）
h_t2=Eo*sum(h2);      % 第二組 Jakes 模型總和（NLOS情況）
%-------------------------- 產生c[k] ------------------------------%
c_k1=h_t1;            % c[k] 第一通道（LOS為主，Rician）
c_k2=h_t2;            % c[k] 第二通道（類Rayleigh）

%-------------------------- 產生d[k] ------------------------------%
bit_num=length(c_k1);           % 資料總長度
d_k=source(bit_num);            % 產生 ±1 的 BPSK 資料

%--------------------- 計算BPSK錯誤率 ---------------------------%
for SNR=0.1:0.1:20
    %----- 模擬通道增益（取樣均值與變異數） -----%
    mean_c_k1=mean(c_k1);
    var_c_k1=std(c_k1)^2;
    Eb=mean(abs(c_k1).^2);
    No=Eb/SNR;
    var_w=No/2;
    Eb_Gauss=mean(abs(c_k1).^2);
    No_Gauss=Eb_Gauss/SNR;
    var_w_Gauss=No_Gauss/2;

    %----- 加入AWGN雜訊 -----%
    w_k1=sqrt(var_w)*randn(1,length(d_k));
    w_k2=sqrt(var_w)*randn(1,length(d_k));
    w_k1_Gauss=sqrt(var_w_Gauss)*randn(1,length(d_k));
    w_k2_Gauss=sqrt(var_w_Gauss)*randn(1,length(d_k));

    %----- 通道增益乘上資料，加上雜訊 -----%
    x_k1=c_k1.*d_k+w_k1;
    x_k2=c_k2.*d_k+w_k2;
    x_k1_Gauss=c_k1.*d_k+w_k1_Gauss;
    x_k2_Gauss=c_k2.*d_k+w_k2_Gauss;

    %------------------ 模擬接收機處理 ------------------%
    % Maximal Ratio Combining (MRC)
    decision_k=real(conj(c_k1).*x_k1 + conj(c_k2).*x_k2);
    decision_k_Gauss=real(conj(c_k1).*x_k1_Gauss + conj(c_k2).*x_k2_Gauss);

    %------------------ 判斷錯誤位元數 ------------------%
    bit_err=0; bit_Gauss_err=0;
    for i=1:length(c_k1)
        if decision_k(i)<0
            d_k_est(i)=-1;
        else
            d_k_est(i)=1;
        end

        if d_k_est(i)~=d_k(i)
            bit_err=bit_err+1;
        end

        if decision_k_Gauss(i)<0
            d_k_Gauss_esti(i)=-1;
        else
            d_k_Gauss_esti(i)=1;
        end

        if d_k_Gauss_esti(i)~=d_k(i)
            bit_Gauss_err=bit_Gauss_err+1;
        end
    end

    Pb(SNR)=bit_err/bit_num;
    Pb_Gauss(SNR)=bit_Gauss_err/bit_num;
end

%----------------------- 畫圖顯示 -----------------------%
SNR=0.1:0.1:20;
semilogy(SNR,Pb,'b.-',SNR,Pb_Gauss,'r.-');
xlabel('Eb/No');
ylabel('BER');
title('BPSK 在 K=5 的 Jakes & Gaussian Rayleigh 通道下的性能比較');
legend('Jakes Rayleigh (K=5)','Gaussian Rayleigh (K=5)');
grid on;
```

---
## 📘 BPSK 傳輸於 Rayleigh (K=5) 通道 — 錯誤率模擬詳解

### 📌 模擬背景與目標

* 模擬 BPSK 在 Rayleigh Fading 通道下的 BER 表現
* 比較 **兩種模型通道**：

  * 使用 **Jakes 模型產生的 Rayleigh 通道**
  * 使用 **理想 Gaussian 分布的 Rayleigh 通道**
* 模擬中的 **K-factor = 5**：來自通道中直視徑功率與散射分量功率的比值

---

### 📌 K-factor = 5 是怎麼設定的？

K-factor 是透過程式碼這一段實現的：

```matlab
Eo = sqrt(10);   % 這裡 Eo^2 = A^2 = 10
Cm = sqrt(1/N);  % 每條路徑功率，總功率=1
h_t1 = Eo * sum(h1);
h_t2 = Eo * sum(h2);
```

對應公式：

$$
K = \frac{A^2}{2\sigma^2} = \frac{10}{2 \cdot 1} = 5
$$

其中：

* \$A^2 = 10\$：直視徑功率（Eo 的平方）
* \$2\sigma^2 = 2\$：兩個 Rayleigh 分支總平均功率

---

### 🔧 程式碼區塊解析

#### 🔹 1. SNR 迴圈（模擬不同的 SNR）

```matlab
for SNR = 0.1:0.1:20
```

---

#### 🔹 2. 通道統計特性：平均功率與雜訊變異數

```matlab
Eb = mean(abs(c_k1).^2);
No = Eb/SNR;
var_w = No / 2;
```

Gaussian 通道對應：

```matlab
Eb_Gauss = mean(abs(c_k1).^2);
No_Gauss = Eb_Gauss / SNR;
var_w_Gauss = No_Gauss / 2;
```

---

#### 🔹 3. AWGN 雜訊加入

```matlab
w_k1 = sqrt(var_w) * randn(1,length(d_k));
w_k1_Gauss = sqrt(var_w_Gauss) * randn(1,length(d_k));
```

---

#### 🔹 4. 傳輸信號（BPSK）經通道與雜訊處理

```matlab
x_k1 = c_k1 .* d_k + w_k1;
x_k1_Gauss = c_k1 .* d_k + w_k1_Gauss;
```

---

#### 🔹 5. Maximal Ratio Combining (MRC)

```matlab
decision_k = real(conj(c_k1).*x_k1 + conj(c_k2).*x_k2);
decision_k_Gauss = real(conj(c_k1).*x_k1_Gauss + conj(c_k2).*x_k2_Gauss);
```

### 💡 MRC 是什麼？

> **最大比率合併（MRC, Maximal Ratio Combining）** 是一種分集技術，
> 可將多個接收支路的訊號依據其通道品質加權後合併，以最大化 SNR。

對應公式：

$$
\hat{d} = \sum_{i=1}^L h_i^* x_i
$$

* \$h\_i\$：第 \$i\$ 條支路的通道增益
* \$x\_i\$：接收訊號
* \$h\_i^\*\$：複數共軛，對應 matched filter

📌 優點：

* 對雜訊弱但通道好者給予較高權重
* 合併後總接收 SNR 提升，有效減少錯誤率

---

#### 🔹 6. 解調與錯誤判斷

```matlab
if decision_k(i) < 0
    d_k_est(i) = -1;
else
    d_k_est(i) = 1;
```

然後與真實 BPSK bit 比對：

```matlab
if d_k_est(i) ~= d_k(i)
    bit_err = bit_err + 1;
end
```

---

#### 🔹 7. 計算錯誤率

```matlab
Pb(SNR) = bit_err / bit_num;
Pb_Gauss(SNR) = bit_Gauss_err / bit_num;
```

---

### 📊 畫圖呈現

```matlab
semilogy(SNR,Pb,'b.-',SNR,Pb_Gauss,'r.-');
legend('Jakes Rayleigh (K=5)','Gaussian Rayleigh (K=5)');
```

* y 軸取 log 顯示 BER
* 比較兩種通道在不同 SNR 下的 BER 曲線

---

### ✅ 總結分析

| 項目          | 說明                                      |
| ----------- | --------------------------------------- |
| 模型          | Rayleigh fading（Jakes 模型 vs Gaussian）   |
| K-factor 設定 | \$K = 5\$，代表有一條較強的 LOS 分量               |
| 技術          | BPSK 傳輸 + MRC 接收                        |
| 評估指標        | Bit Error Rate（BER） vs SNR              |
| 結果預期        | Jakes 模型因建構方式 BER 較差；Gaussian 模型 BER 較好 |

