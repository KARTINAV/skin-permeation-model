#!/usr/bin/env python3
"""Debug f_u calculation."""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from skin_perm.constants import *

# Caffeine: logP=-0.07, pka_HA=14.0, pka_BH+=0.6
logP = -0.07
pka_HA = 14.0
pka_BH = 0.6
has_HA = True
has_BH = True
two_pka = True

# ph_vals: [0, 7, 14, dermis_ph, vehicle_ph, sc_ph]
ph_vals = [0.0, 7.0, 14.0, 7.4, 7.4, 5.0]

f1 = [1.0/(1.0 + 10**(ph - pka_HA)) for ph in ph_vals]
f2 = [1.0/(1.0 + 10**(pka_BH - ph)) for ph in ph_vals]

print("f1 (HA dissociation):")
for i,ph in enumerate(ph_vals):
    print(f"  pH={ph:.1f}: f1={f1[i]:.8f}")
print("f2 (BH+ protonation):")
for i,ph in enumerate(ph_vals):
    print(f"  pH={ph:.1f}: f2={f2[i]:.8f}")

if two_pka:
    Net = [1.0 - (f1[i]+f2[i]) + 2*f1[i]*f2[i] for i in range(6)]
else:
    Net = [1.0] * 6

print("\nNet (non-ionic fraction):")
for i,ph in enumerate(ph_vals):
    print(f"  pH={ph:.1f}: Net={Net[i]:.8f}")

logD = [logP + math.log10(max(Net[i], 1e-30)) for i in range(6)]
print("\nlogD:")
for i,ph in enumerate(ph_vals):
    print(f"  pH={ph:.1f}: logD={logD[i]:.4f}")

# T1 = binding terms
T1 = [0]*5
T1[0] = 0.7936 * math.exp(logP) + 0.2239
T1[1] = 0.5578 * math.exp(logD[3]) + 0.0188
T1[2] = T1[1]
T1[3] = 0.3127 * math.exp(logP) + 0.5121
T1[4] = T1[3]

print("\nT1 (binding terms):")
for i in range(5):
    print(f"  T1[{i}] = {T1[i]:.6f}")

PBR = [t/(1+t) for t in T1]
f_u = [max(1 - p, 0.005) if p < 0.995 else 0.005 for p in PBR]
print("\nPBR, f_u:")
for i in range(5):
    print(f"  PBR[{i}]={PBR[i]:.6f}  f_u[{i}]={f_u[i]:.6f}")

sd3 = logD[0] < logD[1] and logD[1] > logD[2]
print(f"\nsd3 (logD[0]<logD[1] and logD[1]>logD[2]): {sd3}")
print(f"logD[0]={logD[0]:.4f} >= logD[2]={logD[2]:.4f}: {logD[0] >= logD[2]}")

# Selection logic:
if not sd3 and logD[0] >= logD[2]:
    sel_fu = f_u[3]
    print(f"\nSelected: f_u[3] = {sel_fu:.6f}")
else:
    sel_fu = f_u[1]
    print(f"\nSelected: f_u[1] = {sel_fu:.6f}")

# Java output: fu = 0.55
# Let's check what logD Java might compute differently
# For caffeine: Net[0](pH=0) = ?, Net[1](pH=7) = ?, Net[2](pH=14) = ?
# With pka_HA=14, pka_BH=0.6:
#   At pH=0: f1 = 1/(1+10^(0-14)) = 1.0, f2 = 1/(1+10^(0.6-0)) = 0.2009
#   Net[0] = 1 - (1.0+0.2009) + 2*1.0*0.2009 = 1 - 1.2009 + 0.4018 = 0.2009
#   logD[0] = -0.07 + log10(0.2009) = -0.07 - 0.6968 = -0.767
#   At pH=7: f1 = 1/(1+10^(7-14)) = 1/(1+1e-7) ≈ 1.0, f2 = 1/(1+10^(0.6-7)) = ≈ 1.0
#   Net[1] = 1 - 2 + 2 = 1.0, logD[1] = -0.07
#   At pH=14: f1 = 1/(1+10^(14-14)) = 0.5, f2 = 1/(1+10^(0.6-14)) ≈ 1.0
#   Net[2] = 1 - (0.5+1) + 2*0.5*1 = 1 - 1.5 + 1 = 0.5
#   logD[2] = -0.07 + log10(0.5) = -0.07 - 0.301 = -0.371

# Since logD[0]=-0.767 < logD[1]=-0.07: True
# And logD[1]=-0.07 > logD[2]=-0.371: True
# So sd3 = True → should use f_u[1]
# But in Java, f_u[1] corresponds to T1[1] = 0.5578*exp(logD[3])
# logD[3](pH=7.4): Net[3] ≈ 1.0 → logD[3] ≈ -0.07
# T1[1] = 0.5578*exp(-0.07) + 0.0188 = 0.5578*0.9324 + 0.0188 = 0.5201 + 0.0188 = 0.539
# PBR[1] = 0.539/1.539 = 0.3503, f_u[1] = 0.6497

# But Java gives 0.55! There might be a different n_double_bonds or n_ring issue
# Let me check: in-vivo screenshot shows Double Bonds = 4
# Does this affect logP or Kow? No, it affects Va (molar volume)
# Wait - the issue is that I'm using n_double_bonds=4 in formula but
# the Va computation uses it differently
