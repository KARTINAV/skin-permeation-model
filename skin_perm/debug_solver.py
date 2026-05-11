#!/usr/bin/env python3
"""Check the DE profile and bottom flux."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import lil_matrix

Dsc=2.43e-11; Ded=3.36e-7; Dde=3.36e-7
Ksc=42.4; Ked=3.61; Kde=3.61
hsc=1.337e-3; hed=0.01; hde=0.2
nsc=10; ned=10; nde=10; N=30
dxsc=hsc/nsc; dxed=hed/ned; dxde=hde/nde
Csat=1.167
c1=1/(1+Ded*Ked*dxsc/Dsc/Ksc/dxed)
c2=1/(1+Dsc*Ksc*dxed/Ded/Ked/dxsc)
c3=1/(1+Dde*Kde*dxed/Ded/Ked/dxde)
c4=1/(1+Ded*Ked*dxde/Dde/Kde/dxed)

def rhs(t,y):
    dydt=np.zeros(N)
    for i in range(1,nsc-1): dydt[i]=Dsc/dxsc**2*(y[i-1]-2*y[i]+y[i+1])
    dydt[0]=(Csat-y[0])*1e8
    dydt[nsc-1]=Dsc/dxsc**2*(y[nsc-2]+(-1-2*(1-c1))*y[nsc-1]+2*c2*y[nsc])
    dydt[nsc]=Ded/dxed**2*(2*c1*y[nsc-1]+(-1-2*(1-c2))*y[nsc]+y[nsc+1])
    for i in range(nsc+1,nsc+ned-1): dydt[i]=Ded/dxed**2*(y[i-1]-2*y[i]+y[i+1])
    dydt[nsc+ned-1]=Ded/dxed**2*(y[nsc+ned-2]+(-1-2*(1-c3))*y[nsc+ned-1]+2*c4*y[nsc+ned])
    dydt[nsc+ned]=Dde/dxde**2*(2*c3*y[nsc+ned-1]+(-1-2*(1-c4))*y[nsc+ned]+y[nsc+ned+1])
    for i in range(nsc+ned+1,N-1): dydt[i]=Dde/dxde**2*(y[i-1]-2*y[i]+y[i+1])
    dydt[N-1]=Dde/dxde**2*(y[N-2]-2*y[N-1])
    return dydt

Y0=np.zeros(N); Y0[0]=Csat
jac=lil_matrix((N,N),dtype=int)
for i in range(N):
    jac[i,i]=1
    if i>0: jac[i,i-1]=1
    if i<N-1: jac[i,i+1]=1

sol=solve_ivp(rhs,[0,500*3600],Y0,method='BDF',t_eval=[50*3600,200*3600,500*3600],
              rtol=1e-6,atol=1e-10,jac_sparsity=jac)
print(f"Success: {sol.success}")

for ki in range(len(sol.t)):
    t_h = sol.t[ki]/3600
    y = sol.y[:,ki]
    print(f"\n=== t={t_h:.0f}h ===")
    print(f"SC:  Y[0]={y[0]:.4e}  Y[{nsc-1}]={y[nsc-1]:.4e}")
    print(f"VE:  Y[{nsc}]={y[nsc]:.4e}  Y[{nsc+ned-1}]={y[nsc+ned-1]:.4e}")
    print(f"DE:  Y[{nsc+ned}]={y[nsc+ned]:.4e}  Y[{N-1}]={y[N-1]:.4e}")
    # Bottom flux: J = D * dC/dx at x=hde, with C_ghost=0
    # Standard: J = D * Y[N-1] / (dxde/2)  (one-sided with ghost=0)
    Jbot = Dde * y[N-1] / (dxde/2) * 3600  # cm²/s * g/cm³ / cm * s/h = g/cm²/h
    # Scale by Kde/Ksc to convert from SC-normalized
    print(f"Bottom flux (raw)    = {Jbot:.6e} g/cm2/h")
    # Total absorbed flux (in water-equivalent) 
    Qsc = sum(y[:nsc])*dxsc
    Qed = sum(y[nsc:nsc+ned])*dxed*Ked/Ksc
    Qde = sum(y[nsc+ned:N])*dxde*Kde/Ksc
    print(f"Qsc={Qsc:.6e}  Qed={Qed:.6e}  Qde={Qde:.6e}")
    print(f"Qmem_total={Qsc+Qed+Qde:.6e} g/cm2")
