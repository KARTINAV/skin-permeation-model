# Dancik Skin Permeation Model

This repository contains the modernized, high-fidelity Python port of the Dancik et al. (2012) mathematical skin permeation model. It simulates the transient and steady-state absorption of chemicals through the skin (Stratum Corneum, Viable Epidermis, and Dermis compartments).

The legacy Java solver has been completely refactored into a high-performance Python application utilizing SciPy's banded matrix solvers for the implicit Crank-Nicolson integration. The physical and numerical logic natively matches the original literature benchmarks, including complex state-transitions for finite evaporating doses.

## 🌟 Key Features

*   **Robust Transient Solver:** Fully reproduces the original Crank-Nicolson implicit numerical integration.
*   **Exact Mass Balance:** Maintains strict conservative mass balance across all simulation phases (Saturated, Finite, and Depleted).
*   **Built-in Pharmacokinetic Database:** Includes validated properties for Caffeine, Nicotine, and Testosterone (In Vivo & In Vitro).
*   **Modern Web Dashboard:** An interactive, premium UI built with Streamlit and Plotly for real-time visualization of Absorptive Flux, Cumulative Systemic Absorption, and Evaporative Flux.
*   **Data Export:** One-click generation of `.csv` time-series data for downstream statistical analysis.

## 📁 Repository Structure

*   `skin_perm/`: Core module containing the mathematical physics engine.
    *   `solver.py`: The Crank-Nicolson implicit integration logic.
    *   `physics.py`: Partition coefficient (`Kv`), permeability (`Psc`), and thermodynamic (`SW`, `Ro`) calculations.
    *   `app.py`: The Streamlit Web Dashboard application.
*   `standard_benchmark/`: Automated benchmark test scripts used to guarantee 1:1 mathematical parity with the original Java output.
*   `Decompiled_f_4c3l_CDC2020/`: The legacy Java source files retained strictly for historical reference.

## 🚀 Deployment & Usage

The application relies on a modern Python scientific stack (NumPy, SciPy, Pandas, Plotly, and Streamlit).

### 1. Install Requirements
Make sure you have the necessary dependencies installed:
```bash
pip install numpy scipy pandas plotly streamlit
```

### 2. Run the Dashboard locally
To launch the interactive web dashboard on your local machine:
```bash
python -m streamlit run skin_perm/app.py
```
The application will be hosted locally at `http://localhost:8501`.

### 3. Deploying to Production
Because the dashboard is built with Streamlit, it is exceptionally easy to deploy:
*   **Streamlit Community Cloud:** You can directly link this repository to [share.streamlit.io](https://share.streamlit.io) for free, instant global hosting.
*   **Docker:** The application can be easily containerized and deployed to AWS, Google Cloud Run, or Azure App Service.
