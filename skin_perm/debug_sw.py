#!/usr/bin/env python3
"""Debug water solubility calculation."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from skin_perm.chem_data import ChemDat
from skin_perm.skin_data import *
from skin_perm.physics import *
from skin_perm.constants import *

c = ChemDat()
c.set_mw(194.191); c.set_logkow(-0.07); c.set_melting_point(236.0)
c.set_density(1.23); c.set_density_temperature(17.78)
# Java input: 2.17 g/100mL = 21700 mg/L
c.set_water_solubility(21700.0)  # mg/L
c.set_water_solubility_temperature(25.0)

env = EnvOpt()
ro = Ro(c, env)

# Step through SW calculation manually
mp = c.melting_point  # 236
sc_temp = env.temperature  # 32
mw = c.mw  # 194.191
logk = c.logP  # -0.07
sol_temp = c.water_solubility_temperature  # 25

print("=== SW MANUAL CALCULATION ===")
print(f"mp={mp}, sc_temp={sc_temp}, sol_temp={sol_temp}")
print(f"mp > sc_temp: {mp > sc_temp}")
print(f"chem.water_solubility = {c.water_solubility} mg/L")

sg_exp = c.water_solubility  # 21700 mg/L
smol_exp = sg_exp / mw
print(f"sg_exp = {sg_exp}, smol_exp = {smol_exp}")

# Conditions
print(f"\nsc_temp < mp: {sc_temp < mp}")
print(f"mp < sol_temp: {mp < sol_temp}")
print(f"sol_temp < mp < sc_temp: {sol_temp < mp < sc_temp}")

if sc_temp < mp and mp < sol_temp:
    t1 = mp
    print(f"Branch 1: t1 = mp = {t1}")
else:
    t1 = sol_temp
    print(f"Branch 2: t1 = sol_temp = {t1}")

if sc_temp < mp:
    t2 = sc_temp
    print(f"Branch A: t2 = sc_temp = {t2}")
elif sol_temp < mp < sc_temp:
    t2 = mp
    print(f"Branch B: t2 = mp = {t2}")
else:
    t2 = sol_temp
    print(f"Branch C: t2 = sol_temp = {t2}")

smol_2 = smol_exp * 10**(0.01*(t2 - t1))
sg_2 = smol_2 * mw
print(f"\nt1={t1}, t2={t2}")
print(f"10^(0.01*(t2-t1)) = 10^(0.01*{t2-t1}) = {10**(0.01*(t2-t1)):.6f}")
print(f"smol_2 = {smol_2:.6f}")
print(f"sg_2 = {sg_2:.4f} mg/L")
print(f"sg_2/1000 = {sg_2/1000:.4f} g/cm³ vs rho = {ro.eff_density:.4f}")
if sg_2/1000 < ro.eff_density:
    final_sw = sg_2
else:
    final_sw = ro.eff_density * 1000
    print("CAPPED at density!")
print(f"\nFinal SW = {final_sw:.4f} mg/L")

# What Java gets
# Java Csat = 0.109, Java Ksc = 4.28
# → Java SW = 0.109/4.28 = 0.02547 g/cm³ = 25470 mg/L ≈ 25.5 g/L
# But Java shows "Water solubility: 21.7 g/L at 25.0°C" 
# That's the INPUT sw, not the temperature-adjusted one
# So Java's temperature-adjusted SW should be around 25.5 g/L ≈ 25500 mg/L
print(f"\nJava Csat/Ksc = {0.109/4.28:.4f} g/cm³ = {0.109/4.28*1e6:.1f} mg/L")
print(f"This is {0.109/4.28*1e3:.2f} g/L")

# The issue: our code stores WS in mg/L but Java uses different units
# Java input "2.17 g/100mL" = 21.7 g/L = 21700 mg/L
# After T adjustment: 21700 → 25470 at 32°C? No that's too high
# Actually: 10^(0.01*(32-25)) = 10^0.07 = 1.174
# 21700 * 1.174 = 25478 mg/L = 25.5 g/L → Csat = 4.28 * 25.5e-3 = 0.109 ✓
print(f"\nDirect calc: 21700 * 10^(0.07) = {21700*10**0.07:.1f} mg/L")
print(f"= {21700*10**0.07/1000:.2f} g/L")
print(f"Ksc * SW = 4.28 * {21700*10**0.07/1e6:.6f} = {4.28*21700*10**0.07/1e6:.6f}")
