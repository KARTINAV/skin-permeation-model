#!/usr/bin/env python3
"""Compare our steady-state against Java reference for caffeine in-vivo."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from skin_perm.chem_data import ChemDat
from skin_perm.skin_data import *
from skin_perm.physics import *
from skin_perm.m2hp import M2H, M2P
from skin_perm.constants import *
import math

def make_caffeine():
    c = ChemDat()
    c.chem_name = "Caffeine"
    c.set_formula("C8H10N4O2")
    c.set_mw(194.191)
    c.set_logkow(-0.07)
    c.set_melting_point(236.0)
    c.set_boiling_point(178.0)
    c.set_vapour_pressure(0.0000015, scale=TORR_TO_PA)  # 1.5E-6 torr → Pa
    c.set_density(1.23)
    c.set_density_temperature(17.78)
    # Java stores in g/L. Input was 2.17 g/100mL = 21.7 g/L
    c.set_water_solubility(21.7)  # g/L (Java internal unit)
    c.set_water_solubility_temperature(25.0)
    c.grain_class = GRAIN_OTHER
    c.pka_HA = 14.0; c._is_pka_HA = True
    c.pka_BHplus = 0.6; c._is_pka_BHplus = True
    c.n_double_bonds = 4; c.n_triple_bonds = 0; c.n_ring_systems = 2
    return c

def test_case(c, mode_name, in_vivo_vitro, java_ref):
    """Run one test case and compare to Java reference."""
    skin = SkinProp()
    skin.set_hydration(PARTIALLY_HYDRATED)
    skin.set_in_vitro_vivo(in_vivo_vitro)
    env = EnvOpt()
    dose = DosageDat(); dose.set_permeant_amount(500.0); dose.set_area(1.0)
    op = OutputParameters(); op.max_duration = 1000.0
    nvv = VehicleDat(0); vv = VehicleDat(1)

    fu = FuFnom(c, skin, nvv)
    ro = Ro(c, env)
    m2h = M2H(c); m2p = M2P(c)
    psc = Psc(skin, m2p, m2h)
    sw = SW(c, env, ro)
    sys_a = Sys(skin, dose)
    kv = Kv(c, ro, nvv, dose, env, sw, fu, psc)
    po = ProgramOpt(sys_a, env, skin, dose, op, vv)
    sc = SC(c, psc, sw, ro, skin, fu, sys_a, po, kv)
    ev = Evap(c, env, skin, ro, sc, po)

    print(f"\n{'='*60}")
    print(f"  {mode_name}")
    print(f"{'='*60}")
    print(f"{'Parameter':<25} {'Java':>12} {'Ours':>12} {'Match':>6}")
    print("-" * 60)

    all_ok = True
    def cmp(name, java, ours, tol=0.03):
        nonlocal all_ok
        match = abs(java - ours) / max(abs(java), 1e-30) < tol
        all_ok &= match
        print(f"{name:<25} {java:>12.4e} {ours:>12.4e} {'✓' if match else '✗'}")

    cmp("Dsc0 [cm²/s]", java_ref['Dsc0'], sc.sc_Dsc0)
    cmp("Ksc/w", java_ref['Ksc'], sc.est_ksc_w)
    cmp("Ded [cm²/s]", java_ref['Ded'], sc.ve_ded)
    cmp("Ked/w", java_ref['Ked'], sc.ve_kedw)
    cmp("Dde [cm²/s]", java_ref['Dde'], sc.de_dde)
    cmp("Kde", java_ref['Kde'], sc.de_kde)
    cmp("Csat [g/cm³]", java_ref['Csat'], sc.best_csat)
    cmp("Msat [µg/cm²]", java_ref['Msat'], sc.sc_msat)
    cmp("kevap [cm/h]", java_ref['kevap'], ev.kevap)
    cmp("kp_eff [cm/h]", java_ref['kp'], sc.effective_kp)
    cmp("tau_sc [h]", java_ref['tau'], ev.hsc2_dsc0)
    cmp("rho_eff", java_ref['rho'], ro.eff_density)
    cmp("f_u", java_ref['fu'], fu.cf_u_vt)
    cmp("f_non", java_ref['fnon'], fu.cf_non_vt)
    cmp("SW [g/L]", java_ref['sw_gL'], sw.solubility)

    # Additional: thickness, deposition
    print(f"\n  SC: {po.sc_thickness:.3f} µm  VE: {po.ve_thickness:.1f} µm  DE: {po.de_thickness:.1f} µm")
    print(f"  f_dep: {po.deposition_fraction:.4f}")
    print(f"  DE clearance: {sc.de_ked_clear:.4e} s⁻¹")

    return all_ok

# Java reference: In-Vivo
ref_iv = dict(Dsc0=1.76e-12, Ksc=4.28, Ded=1.80e-6, Ked=0.75,
              Dde=1.80e-6, Kde=0.75, Csat=0.109, Msat=3.645,
              kevap=3.370e-9, kp=2.024e-5, tau=282.34, rho=1.216,
              fu=0.55, fnon=1.00, sw_gL=21.7)  # java output says 21.7 g/L

# Java reference: In-Vitro
ref_ivt = dict(Dsc0=1.76e-12, Ksc=4.28, Ded=1.80e-6, Ked=0.75,
               Dde=1.80e-6, Kde=0.75, Csat=0.109, Msat=3.645,
               kevap=3.370e-9, kp=2.025e-5, tau=282.27, rho=1.216,
               fu=0.55, fnon=1.00, sw_gL=21.7)

c = make_caffeine()
ok1 = test_case(c, "CAFFEINE IN-VIVO", IN_VIVO, ref_iv)
ok2 = test_case(c, "CAFFEINE IN-VITRO", IN_VITRO, ref_ivt)

print(f"\n{'='*60}")
if ok1 and ok2:
    print("✅ ALL STEADY-STATE PARAMETERS MATCH JAVA")
else:
    print("❌ SOME PARAMETERS DON'T MATCH")
print(f"{'='*60}")
