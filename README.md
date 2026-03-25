# Financial Econometrics: MVP vs. GARCH Portfolio Optimization

## 📊 Project Overview
[cite_start]This project presents a comprehensive analysis of the **Minimum Variance Portfolio (MVP)** for five major GPW instruments: **PZU, BNP, CDR, ORANGE, and ORLEN**[cite: 5, 18]. [cite_start]The study covers the period 2019-2025 and focuses on comparing the classical static Markowitz approach with dynamic models of multivariate volatility[cite: 6].

[cite_start]The models are verified based on efficiency measures (Sharpe Ratio), extreme risk (VaR, ES), and econometric tests such as Jarque-Bera and Diebold-Mariano[cite: 7].

## 📈 Strategy Performance
The following visualization compares the cumulative returns of the static MVP against dynamic GARCH-class models.

<p align="center">
  <img src="img/equity_curves.png" width="90%" alt="Equity Curves">
</p>
<p align="center"><em>Figure 1: Comparison of cumulative returns for different MVP models (2019-2025).</em></p>

## 🛠 Key Features & Methodology
* [cite_start]**Static MVP:** Optimization using a 500-day rolling window with weekly rebalancing[cite: 31, 32].
* **Dynamic Volatility Models:** Implementation of advanced GARCH-class models:
    * [cite_start]**DCC-GARCH:** Estimated with Normal and Student-t distributions[cite: 83].
    * [cite_start]**BEKK-GARCH:** Capturing volatility spillovers between assets[cite: 85].
* **Econometric Testing:**
    * [cite_start]**Jarque-Bera:** Confirmed strong non-normality and "fat tails" in returns[cite: 91, 94].
    * [cite_start]**Diebold-Mariano:** Comparing accuracy of variance predictions[cite: 95, 96].
    * [cite_start]**VAR(1) Filtration:** Used to eliminate residual autocorrelation[cite: 76, 78].

## ⚠️ Risk Dynamics
[cite_start]The models demonstrate high responsiveness to market shocks, especially during the 2020 pandemic[cite: 150].

<p align="center">
  <img src="img/rolling_risk.png" width="90%" alt="Rolling Risk Analysis">
</p>
<p align="center"><em>Figure 2: 60-day rolling annualized standard deviation. Models adapt effectively to volatility spikes.</em></p>

## 🏆 Results Highlights
* [cite_start]**Risk Reduction:** The MVP strategy successfully reduced total risk (Standard Deviation) by nearly **50%** compared to high-volatility assets like CDR[cite: 40].
* [cite_start]**Best Performer:** **DCC-GARCH (Normal)** achieved the best risk minimization with the lowest VaR at **1.94%**[cite: 107].
* [cite_start]**Efficiency:** While GARCH models provided the lowest risk metrics, the **Static MVP** maintained the highest Sharpe Ratio (0.4437) due to its economic efficiency[cite: 105, 108].

## 📂 Project Structure
* `portfolio_optimization_garch.R`: Main R script with the full analytical pipeline.
* `/data`: CSV datasets for the analyzed GPW instruments.
* `/docs`: Full technical report (in Polish) with detailed statistical interpretations.
* `/img`: Visualizations used in this documentation.

## 🚀 How to Run
1. Clone the repository.
2. Install required R libraries: `quantmod`, `PerformanceAnalytics`, `rugarch`, `rmgarch`, `mgarchBEKK`.
3. Run the script `portfolio_optimization_garch.R`.
