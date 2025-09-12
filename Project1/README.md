## 🌊 Jakes Model 與 Gaussian Rayleigh 的對應與實作筆記

---

### 🧠 主要目的：

用 **Jakes model** 模擬時間變化的 **Rayleigh fading**（即通道係數 $h(t)$ 隨時間變化的隨機過程）。

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

## 📌 Part 1: Jakes model 的 MATLAB 實作（對應圖一）

```matlab
for i=1:N
    phi_m=2*pi*rand(1,1);              % 隨機相位：模擬散射角的差異
    h(i,:)=m*exp(j*(cos(rfa_m(i))*w_m*n'+phi_m)); 
end
```

### ✅ 對應模型說明：

Jakes 模型假設：
通道係數是很多路徑 $h(t) = \sum a_k e^{j\omega_k t + \phi_k}$ 的疊加，
其中：

* 每個路徑來自不同方向 → 會有不同相位與頻率偏移。
* 頻率偏移： $\omega_k = 2\pi f_D \cos(\theta_k)$ ，模擬多徑的 Doppler shift。
* $\phi_k \sim \text{Uniform}(0,2\pi)$ ，模擬多徑隨機性。

這裡就等價於：

```matlab
cos(rfa_m(i))*w_m*n' + phi_m
```

→ 就是每條路徑的 Doppler + 相位。

### 🧪 結果：

產生的 $h(t)$ 是一條 **複數時間序列**，符合 Rayleigh 隨機變化的特性。

---

## 📌 Part 2: Gaussian Rayleigh 模型（圖二）

```matlab
c = sqrt(var/2)*(randn(1,1000000)+j*randn(1,1000000));
```

### ✅ 對應模型說明：

這是另一種「理論建模法」來模擬 Rayleigh：

* 產生 **實部與虛部皆為獨立高斯分佈（mean=0, var/2）**
* 再取 magnitude：$|c|$ 就會形成 **Rayleigh 分布**

這就是：

> Rayleigh fading = 實部 + 虛部為高斯（中心極限定理）

與 Jakes Model 的物理建模（來自路徑和頻移）不同，
這是統計建模（直接產生分布）。

---

## 📌 Part 3: Rayleigh Channel Simulation（圖三）

```matlab
xc = Eo * cos(wn*t + theta(k));
xs = Eo * sin(wn*t + theta(k));
```

這一段也是 Jakes model 的變體版本。
這裡針對每條路徑做：

* 一條餘弦項（實部）：$\cos(\omega_k t + \theta_k)$
* 一條正弦項（虛部）：$\sin(\omega_k t + \theta_k)$

之後將多路徑相加再取絕對值 `abs()` → 得到 fading magnitude

### ✅ 與 Rayleigh 對應：

這種模型在物理上模擬了「同時有多個不同角度入射的波」，
其總和的包絡值（envelope）符合 Rayleigh 分布。

---

## 📊 總結比較表：

| 模型                    | 特點                   | 實作方式                                                 | 對應分布              |
| --------------------- | -------------------- | ---------------------------------------------------- | ----------------- |
| **Jakes Model**       | 從物理角度模擬（Doppler, 多徑） | $\sum e^{j(\omega_k t + \phi_k)}$                    | Rayleigh Envelope |
| **Gaussian Rayleigh** | 統計角度建模               | $\mathcal{N}(0,\sigma^2) + j\mathcal{N}(0,\sigma^2)$ | Rayleigh Envelope |

---

## 🧾 延伸筆記建議：

* 可以比較 histogram（PDF）是否符合 Rayleigh。
* 可以觀察 time domain 中 fading 的快慢（依 Doppler spread）。
* 若增加角度數目 N，則更趨近理論 Rayleigh。

---

✅ 有需要我幫你做出 PDF / DOCX / 繪圖比較也可以繼續說。

