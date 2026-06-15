import math
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

# ==============================================================================
# 1. STAŁE I FUNKCJE POMOCNICZE
# ==============================================================================
G0 = 9.80665
P_ATM = 101325.0

def oblicz_gamma_funkcja(gamma):
    return np.sqrt(gamma * (2.0 / (gamma + 1.0))**((gamma + 1.0) / (gamma - 1.0)))

# ==============================================================================
# 2. FUNKCJE SYMULACYJNE
# ==============================================================================
def symulacja_paliwa_stale(isp, m_paliwa, pc_design_bar, nazwa, n_exp, c_star, rho, a_si, gamma_spalin):
    P_design = pc_design_bar * 1e5
    r_outer = 0.065 / 2
    r_init = 0.012 
    L_total = m_paliwa / (rho * np.pi * (r_outer**2 - r_init**2))
    num_grains = 4
    L_grain = L_total / num_grains
    A_b_init = num_grains * (2 * np.pi * r_init * L_grain + 2 * np.pi * (r_outer**2 - r_init**2))
    A_t = (A_b_init * rho * a_si * (P_design**n_exp) * c_star) / P_design
    
    Gamma = oblicz_gamma_funkcja(gamma_spalin)
    t, dt = 0.0, 0.0005
    Pc = P_ATM
    r_core = r_init
    L_g = L_grain
    V_c = 0.0002 + (num_grains * np.pi * r_core**2 * L_g)
    czasy, ciagi = [], []
    
    while Pc > P_ATM + 10 or (r_core < r_outer and L_g > 0):
        web_burned = r_core - r_init
        sliver_factor = 1.0 if web_burned < 0.94 * (r_outer - r_init) and L_g > 0.06 * L_grain else np.clip((r_outer - r_core) / (0.06 * (r_outer - r_init)), 0, 1) * np.clip(L_g / (0.06 * L_grain), 0, 1)

        if r_core < r_outer and L_g > 0:
            A_b = num_grains * (2 * np.pi * r_core * L_g + 2 * np.pi * (r_outer**2 - r_core**2)) * sliver_factor
            rdot = a_si * (Pc**n_exp)
            ign_factor = 1.0 - np.exp(-t / 0.04)
            mdot_gen = rho * A_b * rdot * ign_factor
            dVc_dt = A_b * rdot
        else:
            mdot_gen, dVc_dt, rdot = 0.0, 0.0, 0.0

        mdot_ext = (Pc * A_t) / c_star
        RT = (Gamma * c_star)**2
        dPc_dt = (RT / V_c) * (mdot_gen - mdot_ext) - (Pc / V_c) * dVc_dt

        Pc += dPc_dt * dt
        if Pc < P_ATM:
            Pc = P_ATM

        r_core += rdot * dt
        L_g -= 2 * rdot * dt
        V_c += dVc_dt * dt
        t += dt

        F = mdot_ext * isp * G0 if Pc > P_ATM else 0.0

        czasy.append(t)
        ciagi.append(F)

    impuls = np.trapezoid(ciagi, czasy)
    return np.array(czasy), np.array(ciagi), A_t, L_total, m_paliwa, r_init * 2, impuls


def symulacja_hybryda(isp_max, m_tot, pc_design_bar, nazwa, of_opt, rho, a_marxman, n_marxman, c_star_max, gamma_spalin):
    P_design = pc_design_bar * 1e5
    m_f_tot = m_tot / (1 + of_opt)
    m_ox_tot = m_tot - m_f_tot

    t_burn_nominal = 4.0
    mdot_ox_nominal = m_ox_tot / t_burn_nominal
    mdot_f_design = mdot_ox_nominal / of_opt

    r_init = 0.015 
    G_ox_init = mdot_ox_nominal / (np.pi * r_init**2)
    rdot_init = a_marxman * (G_ox_init**n_marxman)

    L = mdot_f_design / (2 * np.pi * r_init * rdot_init * rho)
    r_outer = np.sqrt(r_init**2 + m_f_tot / (rho * np.pi * L))

    A_t = ((mdot_ox_nominal + mdot_f_design) * c_star_max) / P_design

    Gamma = oblicz_gamma_funkcja(gamma_spalin)
    t, dt = 0.0, 0.0005
    Pc = P_ATM
    r_port = r_init
    m_ox_used = 0.0
    V_c = 0.0003 + (np.pi * r_port**2 * L)

    czasy, ciagi = [], []

    while Pc > P_ATM + 10 or m_ox_used < m_ox_tot:

        if m_ox_used < m_ox_tot:
            mdot_ox = mdot_ox_nominal * (1.0 - np.exp(-t / 0.06))
            m_ox_used += mdot_ox * dt
        else:
            mdot_ox = 0.0

        if r_port < r_outer:
            A_port = np.pi * r_port**2
            G_ox = mdot_ox / A_port if mdot_ox > 0 else 0
            rdot = a_marxman * (G_ox**n_marxman) if G_ox > 0 else 0
            sliver_factor = np.clip((r_outer - r_port) / (0.02 * (r_outer - r_init)), 0, 1)
            A_b = 2 * np.pi * r_port * L * sliver_factor
            mdot_f = A_b * rdot * rho
            dVc_dt = A_b * rdot
        else:
            rdot, mdot_f, dVc_dt = 0.0, 0.0, 0.0

        mdot_gen = mdot_ox + mdot_f

        if mdot_f > 0:
            of_inst = mdot_ox / mdot_f
            deviation = (of_inst - of_opt) / of_opt
            efektywnosc = np.clip(1.0 - 0.25 * abs(deviation), 0.5, 1.0)
        else:
            efektywnosc = 0.4

        c_star_inst = c_star_max * efektywnosc
        isp_inst = isp_max * efektywnosc

        mdot_ext = (Pc * A_t) / c_star_inst
        RT = (Gamma * c_star_inst)**2
        dPc_dt = (RT / V_c) * (mdot_gen - mdot_ext) - (Pc / V_c) * dVc_dt

        Pc += dPc_dt * dt
        if Pc < P_ATM:
            Pc = P_ATM

        r_port += rdot * dt
        V_c += dVc_dt * dt
        t += dt

        F = mdot_ext * isp_inst * G0 if Pc > P_ATM else 0.0

        czasy.append(t)
        ciagi.append(F)

    impuls = np.trapezoid(ciagi, czasy)
    return np.array(czasy), np.array(ciagi), A_t, L, m_tot, r_init * 2, impuls


# ==============================================================================
# 3. DANE
# ==============================================================================
DATA_PALIW = {
    "KNSB": {"isp": 130, "m_paliwa": 2.446, "pc_bar": 45, "n_exp": 0.32, "C_star": 890, "rho": 1841, "a_si": 6.47e-5, "gamma": 1.13, "type": "solid"},
    "APCP": {"isp": 220, "m_paliwa": 1.456, "pc_bar": 50, "n_exp": 0.42, "C_star": 1550, "rho": 1750, "a_si": 2.15e-5, "gamma": 1.20, "type": "solid"},
    "N2O_Parafina": {"isp_max": 215, "m_tot": 1.699, "pc_bar": 45, "of_opt": 4.8, "C_star_max": 1600, "rho": 900, "a_marxman": 8.0e-5, "n_marxman": 0.62, "gamma": 1.24, "type": "hybrid"},
    "N2O_HTPB": {"isp_max": 225, "m_tot": 1.529, "pc_bar": 50, "of_opt": 6.0, "C_star_max": 1530, "rho": 920, "a_marxman": 3.0e-5, "n_marxman": 0.62, "gamma": 1.25, "type": "hybrid"}
}


# ==============================================================================
# 4. RUN + CSV + PLOT
# ==============================================================================
def uruchom_symulacje():
    os.makedirs("results", exist_ok=True)

    plt.figure(figsize=(11, 6))

    for nazwa, d in DATA_PALIW.items():

        if d["type"] == "solid":
            t, f, at_val, L_val, m_val, d_in_val, impuls = symulacja_paliwa_stale(
                d["isp"], d["m_paliwa"], d["pc_bar"], nazwa,
                d["n_exp"], d["C_star"], d["rho"], d["a_si"], d["gamma"]
            )
        else:
            t, f, at_val, L_val, m_val, d_in_val, impuls = symulacja_hybryda(
                d["isp_max"], d["m_tot"], d["pc_bar"], nazwa,
                d["of_opt"], d["rho"], d["a_marxman"], d["n_marxman"],
                d["C_star_max"], d["gamma"]
            )

        d_sym_mm = (2 * math.sqrt(at_val / math.pi)) * 1000

        print(f"--- {nazwa} ---")
        print(f"Śr. wew. początkowa: {d_in_val*1000:.1f} mm")
        print(f"Długość ziarna:      {L_val*1000:.1f} mm")
        print(f"Masa paliwa:         {m_val:.3f} kg")
        print(f"Impuls całkowity:    {impuls:.2f} Ns")
        print(f"Wymagane A_t:        {at_val*1000000:.2f} mm^2 -> Średnica dyszy: {d_sym_mm:.2f} mm\n")

        plt.plot(t, f, label=f"{nazwa}")

        with open(f"results/{nazwa}_thrust.csv", "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["time", "thrust"])
            writer.writerows(zip(t, f))

    plt.title("Fizyczny profil ciągu silników")
    plt.xlabel("Czas [s]")
    plt.ylabel("Ciąg [N]")
    plt.legend()
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    uruchom_symulacje()
