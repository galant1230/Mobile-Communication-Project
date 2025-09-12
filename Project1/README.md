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
