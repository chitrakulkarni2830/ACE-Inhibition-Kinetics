# ACE Inhibition Kinetics Analysis Dashboard

This project provides a professional diagnostic tool for the kinetic characterization of **Angiotensin-Converting Enzyme (ACE)**. Using computational modeling, it identifies the inhibition mechanism of test compounds through non-linear regression and double-reciprocal transformations.

## 🧪 Scientific Insights
The analysis focuses on identifying how an inhibitor interacts with the enzyme's active site. Based on the experimental data provided:

* **Inhibition Mechanism**: **Competitive**.
* **Diagnostic Proof**: The Lineweaver-Burk plot shows lines intersecting on the $y$-axis ($1/V_{max}$), confirming that $V_{max}$ is unaffected while the apparent $K_m$ increases.
* **Inhibitor Constant ($K_i$)**: **5.82 µM**.
* **Affinity Shift**: The inhibitor increased the $K_m$ from **22.67 µM** to **61.65 µM**.



## 💻 Technical Features
* **Non-Linear Regression**: Built with `scipy.optimize.curve_fit` to provide higher accuracy for $V_{max}$ and $K_m$ than traditional linear-only methods.
* **Cohesive Visualization**: Implements a 3-panel dashboard using `matplotlib.gridspec`, combining raw saturation curves, linear diagnostics, and numerical summaries.
* **Validation**: Automated calculation of $R^2$ values (Goodness of Fit) for all conditions.



## 🚀 Installation & Usage
1. **Clone the repository**:
   ```bash
   git clone [https://github.com/chitrakulkarni2830/ACE-Inhibition-Kinetics.git](https://github.com/chitrakulkarni2830/ACE-Inhibition-Kinetics.git)
   