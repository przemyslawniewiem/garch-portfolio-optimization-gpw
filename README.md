# Financial Econometrics: MVP vs. GARCH Portfolio Optimization

## 📊 Project Overview
[cite_start]This project performs a comprehensive comparative analysis of **Minimum Variance Portfolio (MVP)** strategies for five major instruments on the Warsaw Stock Exchange (GPW): **PZU, BNP, CDR, ORANGE, and ORLEN**[cite: 5]. [cite_start]The study covers the period 2019-2025[cite: 6].

[cite_start]The main goal was to compare the classical static Markowitz approach with dynamic models of multivariate volatility[cite: 6].

## 🛠 Key Features & Methodology
- [cite_start]**Static MVP:** Optimization using a 500-day rolling window with weekly rebalancing[cite: 31, 32].
- [cite_start]**Dynamic Volatility Models:** Implementation of advanced GARCH-class models[cite: 81]:
  - [cite_start]**DCC-GARCH** (Dynamic Conditional Correlation) with Normal and Student-t distributions[cite: 83].
  - [cite_start]**BEKK-GARCH** to capture volatility spillovers between assets[cite: 85].
- [cite_start]**Econometric Testing:** - **Jarque-Bera** test for normality of returns[cite: 7, 91].
  - [cite_start]**Diebold-Mariano** test to compare variance prediction accuracy[cite: 7, 96].
  - [cite_start]**Ljung-Box** test for residual autocorrelation[cite: 78, 79].

## 📈 Results Highlights
- [cite_start]The **MVP strategy** successfully reduced total risk (Standard Deviation) by nearly **50%** compared to the most volatile assets (e.g., CDR)[cite: 40].
- [cite_start]**DCC-GARCH (Normal)** provided the best risk minimization (lowest VaR at 1.94%)[cite: 107].
- [cite_start]All strategies demonstrated robustness during market shocks, including the **2020 pandemic**[cite: 54, 150].


## 📂 Project Structure
- `portfolio_optimization_garch.R`: Main R script containing the full analytical pipeline.
- `/data`: CSV datasets for the analyzed GPW instruments.
- [cite_start]`/docs`: Full technical report (in Polish) with detailed visualizations[cite: 1].

## 🚀 How to Run
1. Clone the repository.
2. Ensure you have the following R libraries: `quantmod`, `PerformanceAnalytics`, `rugarch`, `rmgarch`, `mgarchBEKK`.
3. Run the script `portfolio_optimization_garch.R`.
