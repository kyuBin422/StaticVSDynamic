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
    `MATLAB HMsimPnL` `<X0>` `<p>` `<OutputBool>` `<s0>` `<gamma>` `<kappa>` `<c>` `<F>` `<R>` `<G>` `<b>`
    **where:**
        +`<X0>`  &ndash; float ; initial belief for theta <sub>0</sub> =X0
        +`<p>`  &ndash; bool ; decide whether output figures, 1 is true to output


References
---------------------
1. Keppo, Jussi, Hong Ming Tan, and Chao Zhou. "Smart City Investments." Available at SSRN 3141043 (2019).

2. Keppo, Jussi and Tan, Hong Ming and Zhou, Chao, Smart City Investments under Dynamic Information Acquisition (June    23, 2020). Available at SSRN:[doi:10.2139/ssrn.3141043] (http://dx.doi.org/10.2139/ssrn.3141043)
