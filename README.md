# Rocket Internal Ballistics Simulation

## Overview

This project simulates internal ballistics of rocket engines and compares different propulsion systems under a fixed total impulse constraint (~3000 Ns).

## Engine types

- Solid: KNSB, APCP  
- Hybrid: N2O + Paraffin, N2O + HTPB  

## Goal

To analyze how different propellants affect:
- thrust profile
- burn time
- nozzle sizing
- stability of combustion

## Key feature

The simulation enforces a **constant impulse target (~3000 Ns)** to allow fair comparison between propulsion systems.

## Outputs

- Thrust vs time curves
- Chamber pressure evolution
- Nozzle throat sizing
- Total impulse calculation
