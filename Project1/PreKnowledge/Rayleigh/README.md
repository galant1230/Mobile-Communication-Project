# Rayleigh 分佈從高斯（NLOS 多路徑）推導

## 0. 物理假設（為什麼是 Rayleigh？）
在 **NLOS** 的環境中，接收端看到的等效通道（或雜訊）可以視為許多散射路徑的向量疊加：
- 路徑很多且彼此獨立；
- 每條路徑的到達相位近似 **均勻分佈**（ $[0,2\pi)$ ）；
- 沒有「特別強」的主導直視徑（**沒有 LOS**），功率由眾多路徑分攤。

在此條件下，**同相/正交**分量（I/Q）會因中心極限定理趨近 **零均值高斯**；其包絡（到原點的距離）自然呈 **Rayleigh**。

---

## 1. 數學模型（從 NLOS 多路徑到 I/Q 高斯）
令等效複數通道（或複數雜訊）在某一個固定時間截面為

$$
h \=\ \sum_{k=1}^{L} a_k e^{j\phi_k} \=\ X + jY,
$$

其中 $a_k \ge 0$ 為第 $k$ 條路徑的振幅， $\phi_k \sim \text{Unif}(0,2\pi)$ 且相互獨立， $L$ 很大。

把實部與虛部分開：

$$
X \=\ \sum_{k=1}^{L} a_k \cos\phi_k, 
\qquad
Y \=\ \sum_{k=1}^{L} a_k \sin\phi_k.
$$

因為 $\phi_k$ 均勻且與 $a_k$ 獨立，有

$$
\mathbb{E}[\cos\phi_k] = \mathbb{E}[\sin\phi_k] = 0,\qquad
\mathbb{E}[\cos\phi_k\sin\phi_k]=0,
$$

且

$$
\mathbb{E}[\cos^2\phi_k] = \mathbb{E}[\sin^2\phi_k] = \frac{1}{2}.
$$

---
### 為什麼 \$\mathbb{E}\[\cos^2\phi\_k] = \mathbb{E}\[\sin^2\phi\_k] = \tfrac{1}{2}\$

假設 \$\phi\_k \sim \text{Unif}(0,2\pi)\$。

#### 1. 使用三角恆等式

$$
\cos^2\phi = \tfrac{1}{2}(1+\cos 2\phi), \qquad
\sin^2\phi = \tfrac{1}{2}(1-\cos 2\phi).
$$

#### 2. 對均勻分佈取期望

$$
\mathbb{E}[\cos^2\phi]
= \frac{1}{2\pi}\int_0^{2\pi}\cos^2\phi \, d\phi,
$$

$$
\mathbb{E}[\sin^2\phi]
= \frac{1}{2\pi}\int_0^{2\pi}\sin^2\phi \, d\phi.
$$

代入恆等式：

$$
\mathbb{E}[\cos^2\phi]
= \frac{1}{2\pi}\int_0^{2\pi}\tfrac{1}{2}(1+\cos 2\phi)\, d\phi.
$$

拆成兩項：

$$
= \tfrac{1}{2}\cdot \frac{1}{2\pi}\int_0^{2\pi}1\, d\phi + 
 \tfrac{1}{2}\cdot \frac{1}{2\pi}\int_0^{2\pi}\cos 2\phi\, d\phi.
$$

因為

$$
\int_0^{2\pi}\cos 2\phi \, d\phi = 0,
$$

所以：

$$
\mathbb{E}[\cos^2\phi] = \tfrac{1}{2}.
$$

同理：

$$
\mathbb{E}[\sin^2\phi] = \tfrac{1}{2}.
$$


#### 3. 結果

因此：

$$
\mathbb{E}[\cos^2\phi_k] = \mathbb{E}[\sin^2\phi_k] = \tfrac{1}{2}.
$$

👉 直觀上，因為 \$\phi\$ 在 $\[0,2\pi)\$ 上均勻分佈，平均來看， \$\cos^2\$ 和 \$\sin^2\$ 分別佔一半的能量。

---

因此

$$
\mathbb{E}[X]=\mathbb{E}[Y]=0,
$$

$$
\mathrm{Var}(X)=\sum_{k=1}^{L}\mathbb{E}[a_k^2]\cdot \mathbb{E}[\cos^2\phi_k]
= \tfrac{1}{2}\sum_{k=1}^{L}\mathbb{E}[a_k^2],
$$

$$
\mathrm{Var}(Y)=\tfrac{1}{2}\sum_{k=1}^{L}\mathbb{E}[a_k^2],
\qquad
\mathrm{Cov}(X,Y)=0.
$$

當 $L$ 很大（許多獨立小貢獻），由**中心極限定理**可得

$$
X \sim \mathcal N(0,\sigma^2),\qquad Y \sim \mathcal N(0,\sigma^2),\qquad
\sigma^2 \=\ \frac{1}{2}\sum_{k=1}^{L}\mathbb{E}[a_k^2],
$$

且 $X$ 與 $Y$ 獨立同分佈。這稱為**圓對稱複高斯**（CSCG）隨機變數。

---

## 2. 由 I/Q 高斯推導包絡的 Rayleigh 分佈（座標轉換）
定義包絡與相位：

$$
R = \sqrt{X^2+Y^2},\qquad \Theta=\mathrm{atan2}(Y,X).
$$

因 $X,Y$ 獨立且同為 $\mathcal N(0,\sigma^2)$，其聯合密度為

$$
f_{X,Y}(x,y)=\frac{1}{2\pi\sigma^2}\exp\\left(-\frac{x^2+y^2}{2\sigma^2}\right).
$$

---
### 為什麼叫「二維圓對稱高斯」？

我們的聯合分佈是：

$$
f_{X,Y}(x,y) = \frac{1}{2\pi\sigma^2} \exp\!\left(-\frac{x^2+y^2}{2\sigma^2}\right).
$$

#### 觀察重點

* 指數項只跟 \$x^2+y^2\$ 有關，也就是「到原點的距離平方」：

$$
r^2 = x^2+y^2.
$$

* 所以只要半徑 \$r\$ 一樣，密度值就一樣。
* 等機率密度的曲線滿足：

$$
x^2+y^2=c \quad (\text{常數}).
$$

* 這正是一個「以原點為圓心的圓」。

#### 結論

因此它是 **圓對稱 (circularly symmetric)** 的二維高斯分佈。

### ASCII 圓環示意

```text
          y
          ↑
      *********
    **         **
   *             *
   *      ●      *+----------------→ x   ← 等密度圓 (x²+y² = c₁)
   *             *
    **         **
      *********

    ***************
  **               **
 *                   *
 *        ●          *   ← 更大一圈等密度圓 (x²+y² = c₂ > c₁)
 *                   *
  **               **
    ***************

```

👉 每一圈都是 \$x^2+y^2=c\$ 的等密度線，半徑越大，機率密度越小。

---

改成極座標（ $x=r\cos\theta,\, y=r\sin\theta$， $r\ge0, \theta\in(-\pi,\pi]$ ），Jacobian 為

$$
\left|\frac{\partial(x,y)}{\partial(r,\theta)}\right| = r.
$$

因此

$$
f_{R,\Theta}(r,\theta)
= f_{X,Y}(x,y)\cdot r
= \frac{1}{2\pi\sigma^2}\exp\\left(-\frac{r^2}{2\sigma^2}\right)\ r.
$$

邊際化掉 $\Theta$：

$$
f_R(r) = \int_{-\pi}^{\pi} f_{R,\Theta}(r,\theta)\, d\theta
= \frac{r}{\sigma^2}\exp\\left(-\frac{r^2}{2\sigma^2}\right),\qquad r\ge 0.
$$

這就是 **Rayleigh 分佈** 的 PDF：

$$
\boxed{\ f_R(r) \=\ \dfrac{r}{\sigma^2}\ e^{-\,r^2/(2\sigma^2)}\quad r\ge0 \ }.
$$

其 CDF、均值、方差分別為

$$
F_R(r)=1-\exp\\left(-\frac{r^2}{2\sigma^2}\right),\quad r\ge0,
$$

$$
\mathbb{E}[R]=\sigma\sqrt{\tfrac{\pi}{2}},\qquad
\mathrm{Var}(R)=\left(2-\tfrac{\pi}{2}\right)\sigma^2.
$$

---

## 3. 另一種極簡觀點： $\chi^2$ / 指數
由於

$$
\frac{X^2}{\sigma^2} + \frac{Y^2}{\sigma^2} \;\sim\; \chi^2_{(2)},
$$

等價於一個參數為 $\lambda=\tfrac{1}{2}$ 的指數分佈，因此

$$
\frac{R^2}{\sigma^2} \sim \chi^2_{(2)}
\quad\Longrightarrow\quad
R^2 \sim \text{Exp}\!\left(\frac{1}{2\sigma^2}\right).
$$

把 $Z=R^2$ 的指數密度 $f_Z(z)=\tfrac{1}{2\sigma^2}e^{-z/(2\sigma^2)}$（ $z\ge0$ ）做變數回代
（ $z=r^2,\, dz=2r\,dr$ ）即可到達同一個 Rayleigh PDF。

---

## 4. 為何強調「NLOS」？
- **NLOS** 多路徑 ⇒ 許多隨機相位的微小向量相加 ⇒ $X,Y$ 零均值、等方差、獨立高斯 ⇒ **Rayleigh 包絡**。
- **若存在強 LOS**（或主導徑）， $h$ 會有**非零均值**， $R$ 的分佈就變成 **Rician**：
  
### Rician 分佈：從有均值的複高斯推導

### 問題設定

有 LOS（或主導徑）時，複通道可寫成

$$
h = X + jY = \mu + n,\qquad n\sim \mathcal{CN}(0,\sigma^2),
$$

其中 \$\mu\in\mathbb{C}\$ 為**固定**（非隨機）的 LOS 成分，\$n\$ 為零均值複高斯雜訊，且

$$
X\sim\mathcal N(\mu_x,\sigma^2),\qquad Y\sim\mathcal N(\mu_y,\sigma^2),\qquad X\perp Y,
$$

令 \$s\equiv|\mu|=\sqrt{\mu\_x^2+\mu\_y^2}\$，並記 \$\mu\_x=s\cos\theta\_0,\ \mu\_y=s\sin\theta\_0\$。
我們欲推導包絡 \$R=|h|=\sqrt{X^2+Y^2}\$ 的分佈。

---

### 第一步：二維高斯（有均值）的聯合密度

$$
f_{X,Y}(x,y) = \frac{1}{2\pi\sigma^2}\,\exp\!\left(-\frac{(x-\mu_x)^2+(y-\mu_y)^2}{2\sigma^2}\right).
$$

將平方和展開、並以 \$r,\theta\$ 表示：

$$
(x-\mu_x)^2+(y-\mu_y)^2
= x^2+y^2 + \mu_x^2+\mu_y^2 - 2(x\mu_x+y\mu_y)
= r^2 + s^2 - 2rs\cos(\theta-\theta_0),
$$

其中 \$x=r\cos\theta,\ y=r\sin\theta\$。

---

### 第二步：極座標變數變換

Jacobian 為 \$\left|\partial(x,y)/\partial(r,\theta)\right|=r\$，故

$$
\begin{aligned}
f_{R,\Theta}(r,\theta)
&= f_{X,Y}(x,y)\cdot r\\
&= \frac{1}{2\pi\sigma^2}\,\exp\!\left(-\frac{r^2+s^2-2rs\cos(\theta-\theta_0)}{2\sigma^2}\right)\,r.
\end{aligned}
$$

---

### 第三步：對角度積分，得到 \$R\$ 的邊際分佈

因為只關心 \$R\$，對 \$\theta\in(-\pi,\pi]\$ 積分：

$$
\begin{aligned}
f_R(r)
&= \int_{-\pi}^{\pi} f_{R,\Theta}(r,\theta)\,d\theta\\
&= \frac{r}{2\pi\sigma^2} e^{-(r^2+s^2)/(2\sigma^2)}
\int_{-\pi}^{\pi} \exp\!\left( \frac{rs}{\sigma^2}\cos(\theta-\theta_0) \right) d\theta.
\end{aligned}
$$

利用 \$\int\_{-\pi}^{\pi} e^{\kappa\cos u},du = 2\pi I\_0(\kappa)\$ 與平移不影響積分（ \$u=\theta-\theta\_0\$ ），得

$$
\boxed{\; f_R(r) = \frac{r}{\sigma^2}\\exp\!\left(-\frac{r^2+s^2}{2\sigma^2}\right) I_0\\left(\frac{rs}{\sigma^2}\right),\quad r\ge 0.\;}
$$

這就是 **Rician 分佈** 的 PDF。

> 註：當 \$s\to 0\$（無 LOS）時， \$I\_0(0)=1\$ ，即退化為 Rayleigh： \$f\_R(r)=\frac{r}{\sigma^2}e^{-r^2/(2\sigma^2)}\$ 。

---

### CDF（以 Marcum Q 函數表示）

Rician 的 CDF 可優雅地寫為

$$
F_R(r) = \Pr(R\le r) = 1 - Q_1\!\left(\frac{s}{\sigma},\frac{r}{\sigma}\right),\qquad r\ge 0,
$$

其中 \$Q\_1(\cdot,\cdot)\$ 為一階 Marcum \$Q\$ 函數。

---

### 物理/幾何意義（直觀）

* \$n\$ 在平面上形成以原點為中心的高斯「隨機雲」。
* 常數向量 \$\mu\$ 把這個雲**整體平移**到以 \$\mu\$ 為中心。
* 從原點量到點 \$(X,Y)\$ 的距離 \$R\$ 就形成「偏心的幅度分佈」—這正是 Rician。

---

### 參數與常見特例

*  \$s=|\mu|\$ ：LOS 強度； \$\sigma^2\$ ：每個 I/Q 分量的變異數。
* 常見以 \$K\$ 因子表徵 LOS 與散射能量比： \$K\equiv \dfrac{s^2}{2\sigma^2}\$ .（ \$K\to 0\$：Rayleigh； \$K\to \infty\$ ：近似常數幅度）


---

## 5. 一句話總結
**Rayleigh 分佈** 來自 **NLOS 多路徑** 下許多隨機相位路徑向量的疊加：

I/Q 由 **中心極限定理** 成為零均值獨立高斯，故其包絡（到原點距離）必然是 **Rayleigh**。

---
## NLOS 與 LOS 下的通道分佈差異

### 1. NLOS（純隨機散射、多徑）

* 每條路徑的相位 \$\phi\_k \sim \text{Unif}(0,2\pi)\$，彼此獨立。
* 很多小路徑相加後，實部 \$X\$ 和虛部 \$Y\$ 會因 **中心極限定理** 變成零均值高斯。
* 對稱性保證：

$$
\mathbb{E}[\cos\phi_k] = \mathbb{E}[\sin\phi_k] = 0.
$$

* 所以整體通道係數為

$$
h = X + jY,
$$

  其均值為 \$0\$。
* 對應的包絡分佈就是 **Rayleigh 分佈**。

---

### 2. 有 LOS（Line-of-Sight）

* 當有一條直視徑存在時，情況不同：

  * 直視徑有「固定的振幅與相位」，不像散射徑那樣隨機。
  * 可以看作是一個常數向量 \$\mu\$。

* 通道可寫成：

$$
h = \underbrace{\mu}_{\text{LOS 固定成分}} + \underbrace{(\text{很多隨機小路徑})}_{\sim \mathcal{CN}(0,\sigma^2)}.
$$

* 隨機小路徑仍然平均為零，但 LOS 成分 \$\mu\$ 會把整體均值偏移到 \$\mu\$：

$$
\mathbb{E}[h] = \mu \ne 0.
$$

* 此時包絡分佈對應 **Rician 分佈**。
