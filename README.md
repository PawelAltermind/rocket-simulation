# 🚀 Rocket Engine Internal Ballistics Simulation

Numerical simulation of internal ballistics for solid and hybrid rocket motors.

The project compares different propulsion systems under a fixed design constraint:

> **Total impulse ≈ 3000 Ns**

---

## 🎯 Project Goal

The main objective is to analyze how different propellants influence:

- thrust profile
- burn time
- nozzle geometry requirements
- total impulse
- stability of combustion

Four propulsion configurations are simulated:

- KNSB (solid)
- APCP (solid)
- N₂O / Paraffin (hybrid)
- N₂O / HTPB (hybrid)

---

## ⚙️ Model Overview

The simulation is based on simplified internal ballistics equations:

- Saint Robert’s law (solid motors)
- Marxman model (hybrid motors)
- time-dependent chamber pressure solver
- evolving grain geometry

---

## 📊 Key Results

### Thrust profiles

![Thrust profile](figures/wykres_idealny_Ft.png)

### Chamber pressure

![Pressure](figures/fizyczny_profil_cisnienia.png)

### Hybrid performance comparison

![Hybrid](figures/profil_ciagu_3kN_realistyczny.png)

---

## 📦 Project Structure
