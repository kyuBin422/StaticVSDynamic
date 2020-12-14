Description
---------------------

Simulate the Brown point in the dynamic model and evaluate with the static model.

### Implementation Notes

- Programming language: **MATLAB** and **C++**
  - Implement in the MATLAB 2019b
  - Microsft C++ compiler
  - IDE Microsoft Visual Studio 2019 Community 
- Intersection Point Algorithm
  - Compute the intersection point with line segement and polygon

Package Contents
---------------------

### Script Code

The script code is located in parent dictionary. there are **four** main MATLAB script:
1. [`HMsimPnL`](HMsimPnL.m) script simulate the dynamic model by Brownian motion and output the optimal q and theta.
  - this formula for boundary construction is described in [[1]](#references).
  - **Usage:**
    `HMsimPnL` `<X0>` `<p>` `<OutputBool>` `<s0>` `<gamma>` `<kappa>` `<c>` `<F>` `<R>` `<G>` `<b>`
   - **Where:**
        `<X0>`  &ndash; float ; initial belief for theta <sub>0</sub> =X0
        `<p>`  &ndash; bool ; decide whether output figures, 1 is true to output
        `<OutputBool>`  &ndash; bool ; decide output result
        `<s0>`  &ndash; float ; variance under q=0
        `<gamma>`  &ndash; float ; risk aversion
        `<kappa>`  &ndash; float ; initial wealth
        `<c>`  &ndash; float ; marginal cost of information
        `<F>`  &ndash; float ; fixed cost of information
        `<R>`  &ndash; float ; lending and borrowing payoff multiplier
        `<G>`  &ndash; float ; goverment variable
        `<b>`  &ndash; float ; another goverment variable 
        
2. [`staticVSDynamic`](staticVSDynamic.m) script simulate the performance of dynamic model and static model
3. [`measurePDF`](measurePDF.m) estimate the probability density function of optimal q and thetaq by normal, gamma, log normal, weibull and GH distrbution
4. [`calculateStoppingIndex`](calculateStoppingIndex.cpp) compute the stopping index by C++ language

References
---------------------
1. Keppo, Jussi and Tan, Hong Ming and Zhou, Chao, Smart City Investments under Dynamic Information Acquisition (June    23, 2020). Available at SSRN:[doi:10.2139/ssrn.3141043] (http://dx.doi.org/10.2139/ssrn.3141043)
