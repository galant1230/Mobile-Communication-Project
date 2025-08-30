# 📡 行動通訊課程專案 (Mobile Communication Projects)

本專案包含兩個行動通訊課程中的模擬實驗，重點在於 **通道建模** 與 **OFDM 系統模擬**。

---

## 🔹 Project 1 — Rayleigh/Rician 通道模擬與 BPSK 傳輸效能分析

- 採用 **Jake’s Model** 產生瑞利衰落 (Rayleigh fading) 通道，並比較模擬結果與理論 PDF。  
- 嘗試不同參數 (路徑數 N、速度 v、時間長度 T) 對通道行為的影響。  
- 擴展至 **Rician fading** (不同 K 因子)。  
- 在此通道下模擬 **BPSK 調變** 的 BER 表現，並與高斯隨機通道模型比較。  
- 探討 **MRC 多重接收分集** (L=2,3) 下的效能改善。  

📈 **主要成果**：  
- Fading 分佈直方圖 vs 理論 PDF  
- 不同 K 值下的 Rician channel 模擬  
- BER 曲線比較 (Rayleigh / Rician / Gaussian / 分集)  

---

## 🔹 Project 2 — OFDM 系統模擬與等化實驗

- 實作一個 **OFDM 收發機模擬** (M=64 FFT, 52 個子載波, CP=16)。  
- 採用 **BPSK 調變**，在多路徑衰落通道下進行傳輸。  
- 分析等化前後的星座圖 (scatter plot)，比較 **頻率選擇性衰落對 OFDM 的影響**。  
- 使用 **Zero-Forcing 頻域等化器**還原訊號。  
- 繪製不同 \(E_b/N_0\) 下的 **BER 曲線**，比較單一路徑通道與多路徑通道。  

📈 **主要成果**：  
- Constellation 圖 (等化前後對比)  
- BER 曲線在不同通道條件下的比較  

---

👉 總結：
- **Project 1**：專注於 **無線通道建模與 BPSK 性能分析**  
- **Project 2**：專注於 **OFDM 系統設計、通道效應與等化補償**  

---
