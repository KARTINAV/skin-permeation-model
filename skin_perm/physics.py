"""
Psc.java, SW.java, Ro.java, Fu_Fnom.java, Kv.java, Evap.java,
ProgramOpt.java, Sys.java, SC.java → Python.
All the intermediate computation classes.
"""
import math
from .constants import *
from .chem_data import ChemDat
from .skin_data import SkinProp, EnvOpt, DosageDat, VehicleDat, OutputParameters
from .m2hp import M2H, M2P


class Psc:
    """Select correct SC permeability expression based on (sigma,r) region."""
    def __init__(self, skin: SkinProp, m2p: M2P, m2h: M2H):
        ls_fh, lr_fh = m2h.log10_s, m2h.log10_r
        cat_fh = self._categorize(ls_fh, lr_fh)
        ls_ph, lr_ph = m2p.log10_s, m2p.log10_r
        cat_ph = self._categorize(ls_ph, lr_ph)

        ksc_fh = m2h.ksc_w
        ksc_ph = m2p.ksc_w
        self.ksc_final = ksc_fh if skin.hydration == FULLY_HYDRATED else ksc_ph
        self.ksc_scaling_factor = ksc_ph / ksc_fh if ksc_fh != 0 else 1.0

        psc_fh = self._select_psc_fh(cat_fh, m2h)
        psc_ph = self._select_psc_ph(cat_ph, m2p)
        self.psc_scaling_factor = psc_ph / psc_fh if psc_fh != 0 else 1.0

        if skin.hydration == FULLY_HYDRATED:
            self.psc_final = psc_fh
        else:
            self.psc_final = psc_ph

    @staticmethod
    def _categorize(ls, lr):
        e4 = ls < -3.5; e5 = -3.5 <= ls <= 1.5; e6 = 1.5 < ls < 4.5; e7 = ls >= 4.5
        e10 = lr < -0.5; e11 = -0.5 <= lr <= 1.0; e12 = 1.0 < lr < 4.5
        e13 = 2.7 < lr < 4.5; e14 = lr >= 4.5
        cat = 1
        if e5 and e12: cat = 6
        elif e6 and e12: cat = 7
        elif e7 and e14: cat = 5
        elif e10: cat = 4
        elif e11 or (e12 and e7): cat = 3
        elif e4 and e13: cat = 2
        return cat

    @staticmethod
    def _select_psc_fh(cat, m):
        return {2: m.pscw_analytic_hour, 3: m.pcompscw, 4: m.pcompscw_a3b_lor,
                5: m.pcompscw_a3b_hir, 6: m.psc_hour, 7: m.psc_hour_hole}.get(cat, m.psc_hour)

    @staticmethod
    def _select_psc_ph(cat, m):
        return {2: m.pscw_analytic_hour, 3: m.pcompscw, 4: m.pcompscw_a3b_lor,
                5: m.pcompscw_a3b_hir, 6: m.psc_hour, 7: m.psc_hour_hole}.get(cat, m.psc_hour)


class Ro:
    """Effective density at skin temperature (Ro.java)."""
    def __init__(self, chem: ChemDat, env: EnvOpt):
        va_bp = chem.get_va()
        ts = env.temperature + 273.15
        mw_a = chem.mw
        bp = (chem.boiling_point + 273.15) if chem._is_bp else 0.0
        exp_map = {GRAIN_OTHER: 0.31, GRAIN_ALCOHOL: 0.25, GRAIN_HYDROCARBON: 0.29}
        exponent = exp_map.get(chem.grain_class, 0.31)

        aws = (va_bp != MINUS_BIG and bp > 0.0 and not chem.formula.has_other())
        cp_ts = mw_a * (3.0 - 2.0*ts/bp)**exponent / va_bp if (aws and bp > 0) else 0.0

        is_ro_ts = chem._is_density and chem._is_density_temp
        ro_ts = 0.0
        if is_ro_ts and bp > 0:
            u_dens_temp = chem.density_temperature + 273.15
            denom = 3.0*bp - 2.0*u_dens_temp
            if denom != 0:
                ro_ts = chem.density * ((3.0*bp - 2.0*ts) / denom)**exponent

        if is_ro_ts: self.eff_density = ro_ts
        elif aws: self.eff_density = cp_ts
        elif chem._is_density: self.eff_density = chem.density
        else: self.eff_density = 1.0


class SW:
    """Water solubility at skin temperature (SW.java)."""
    def __init__(self, chem: ChemDat, env: EnvOpt, ro: Ro):
        mp = chem.melting_point
        sc_temp = env.temperature
        mw = chem.mw
        logk = chem.logP
        sol_temp = chem.water_solubility_temperature if chem._is_ws_temp else 0.0

        mol_25 = 10**(0.5 - logk - 0.01*(mp - 25.0))
        mol_32_a = mol_25 * 10**0.07
        mol_32_b = 10**(0.5 - logk)

        g_32_a = mol_32_a * mw if (mol_32_a*mw/1000 < ro.eff_density) else ro.eff_density*1000
        g_32_b = mol_32_b * mw if (mol_32_b*mw/1000 < ro.eff_density) else ro.eff_density*1000
        csw = g_32_b if mp < sc_temp else g_32_a

        if chem._is_ws and chem._is_ws_temp:
            sg_exp = chem.water_solubility
            smol_exp = sg_exp / mw
            if sc_temp < mp and mp < sol_temp: t1 = mp
            else: t1 = sol_temp
            if sc_temp < mp: t2 = sc_temp
            elif sol_temp < mp < sc_temp: t2 = mp
            else: t2 = sol_temp
            smol_2 = smol_exp * 10**(0.01*(t2 - t1))
            sg_2 = smol_2 * mw
            self.solubility = sg_2 if sg_2/1000 < ro.eff_density else ro.eff_density*1000
        else:
            self.solubility = csw


class FuFnom:
    """Ionization fractions and protein binding (Fu_Fnom.java)."""
    def __init__(self, chem: ChemDat, skin: SkinProp, nvv: VehicleDat):
        has_HA = chem._is_pka_HA
        has_BH = chem._is_pka_BHplus
        two_pka = has_HA and has_BH
        no_pka = (not has_HA) and (not has_BH)

        ph_vals = [0.0, 7.0, 14.0, skin.dermis_ph, nvv.vehicle_ph, skin.sc_ph]
        f1 = [1.0/(1.0 + 10**(ph - chem.pka_HA)) for ph in ph_vals] if has_HA else [-1.0]*6
        f2 = [1.0/(1.0 + 10**(chem.pka_BHplus - ph)) for ph in ph_vals] if has_BH else [0.0]*6

        if two_pka:
            Net = [1.0 - (f1[i]+f2[i]) + 2*f1[i]*f2[i] for i in range(6)]
        elif has_HA:
            Net = [f1[i] for i in range(6)]
        elif has_BH:
            Net = [f2[i] for i in range(6)]
        else:
            Net = [1.0]*6

        logD = [chem.logP + math.log10(max(Net[i], 1e-30)) for i in range(6)]

        T1 = [0]*5
        T1[0] = 0.7936 * math.exp(chem.logP) + 0.2239
        T1[1] = 0.5578 * math.exp(logD[3]) + 0.0188
        T1[2] = T1[1]
        T1[3] = 0.3127 * math.exp(chem.logP) + 0.5121
        T1[4] = T1[3]

        PBR = [t/(1+t) for t in T1]
        f_u = [max(1 - p, 0.005) if p < 0.995 else 0.005 for p in PBR]

        # Java sd[] flags:
        sd3 = logD[0] < logD[1] and logD[1] > logD[2]  # zwitterion peak
        sd1 = (not sd3) and (logD[0] == logD[2])
        sd2 = (not sd3) and (logD[0] < logD[2])

        self.cf_non_vt = chem.f_non_vt if chem._is_f_non_vt else Net[3]
        if chem._is_f_u_vt:
            self.cf_u_vt = chem.f_u_vt
        elif (not sd1) and (not sd2):
            # Java: if (!sd[1] && !sd[2]) → use f_u[3]
            self.cf_u_vt = f_u[3]
        else:
            self.cf_u_vt = f_u[1]

        self.cf_non_veh = chem.f_non_veh if chem._is_f_non_veh else Net[4]
        self.cf_non_sc = chem.f_non_sc if chem._is_f_non_sc else Net[5]
        self.cf_u_sc = chem.f_u_sc


class Sys:
    """Systemic compartment (Sys.java). All metabolism/clearance rates default 0."""
    def __init__(self, skin: SkinProp, dp: DosageDat):
        self.ksc_met = 0.0; self.kve_met = 0.0; self.kde_met = 0.0
        self.kblood_met = 0.0; self.kblood_fat = 0.0; self.kfat_blood = 0.0
        self.kelim = 0.0; self.vblood = 0.0; self.vfat = 0.0


class Kv:
    """Vehicle/SC partition coefficient (Kv.java)."""
    def __init__(self, chem: ChemDat, ro: Ro, nvv: VehicleDat, dose: DosageDat,
                 env: EnvOpt, sw: SW, ff: FuFnom, psc: Psc):
        density = ro.eff_density
        mp = chem.melting_point
        mw = chem.mw
        temperature = env.temperature
        best_sol = sw.solubility / 1000.0

        f_non_veh = chem.f_non_veh if chem._is_f_non_veh else ff.cf_non_veh
        nvv_sol = (chem.non_volatile_vehicle_solubility / 1000.0
                   if (nvv.is_setup() and chem._is_nvv_sol) else 0.0)
        nvv_aa = (dose.nvv_amount_applied
                  if (nvv.is_setup() and dose.is_nvv_aa) else 0.0)

        best_ku = nvv_sol / best_sol if nvv_sol != 0 else 0.0
        if nvv.is_water:
            k_v_w = 1.0 / f_non_veh if f_non_veh != 0 else 1e10
        elif nvv_aa > 0:
            k_v_w = best_ku
        else:
            k_v_w = density / best_sol if best_sol != 0 else 1e10
        self.k_v_w = k_v_w
        self.nvv_amount_applied = nvv_aa

        # Vehicle volume and Vsat
        nvv_density = nvv.density if (nvv.is_setup() and nvv._is_density) else 0.0
        hv = 10000*nvv_aa/1000/nvv_density if (nvv_aa > 0 and nvv_density > 0) else 0.0
        self.v_sat = nvv_sol * 1e6 * (hv/10000) if hv > 0 else 0.0

        # Ideal saturation solubility in octanol
        r_gas = 1.98739632
        r_t = r_gas * (temperature + 273.15)
        logP = chem.logP
        oct_density = 6.9e-4 * temperature + 0.813
        oct_mw = 130.22
        k_oct = 10**logP
        ln_s_oct = -13.5 * (mp - temperature) / r_t if mp > temperature else 0.0
        s_oct = math.exp(ln_s_oct)
        self.ist_s_oct = ((oct_density*(1-s_oct) + density*s_oct) * s_oct * mw /
                          (mw*s_oct + (1-s_oct)*oct_mw))


class ProgramOpt:
    """Simulation control parameters (ProgramOpt.java)."""
    def __init__(self, sys_obj: Sys, env: EnvOpt, skin: SkinProp, dp: DosageDat,
                 op: OutputParameters, vv: VehicleDat):
        self.R = 0.082056
        self.area = dp.area
        self.headspace_volume = 0.0
        self.method = 3
        self.num_sublayers = 120
        self.init_step_size = 1e-6
        self.minimum_step_size = 1e-8
        self.maximum_step_size = op.max_step_size
        self.dy_max = 0.01
        self.output_flag = op.output_flag
        self.max_output_step_size = op.min_output_step_size
        self.wash = dp.is_terminal_wash

        self.sc_thickness = skin.sc_thickness if skin.is_sc else 0.0
        self.ve_thickness = skin.ve_thickness if skin.is_ve else 0.0
        self.de_thickness = skin.dermis_thickness if skin.is_dermis else 0.0

        if skin.is_sc:
            self.deposition_fraction = 0.1 if vv.is_present else 3.0 / self.num_sublayers
        else:
            self.deposition_fraction = 0.0
        self.deposition = self.deposition_fraction * self.sc_thickness


class SC:
    """SC / VE / DE properties — the heart of the steady-state model (SC.java)."""
    def __init__(self, chem: ChemDat, psc_obj: Psc, sw_obj: SW, ro: Ro,
                 skin: SkinProp, ff: FuFnom, sys_obj: Sys, po: ProgramOpt, kv: Kv):
        density = ro.eff_density
        mw = chem.mw
        logP = chem.logP
        k_oct = 10**logP
        potts_guy = 10**(-2.7437 + 0.71*logP - 0.0061*mw)

        # ── firstSC ──
        jmn = psc_obj.psc_final
        exp_kp = chem.permeability_coefficient if chem._is_kp else 0.0
        psc_scale = psc_obj.psc_scaling_factor

        # ── calculateDE ──
        de_kfree, de_falbumin, de_flipid = 0.6, 0.32, 0.001
        de_Q, de_S, de_m, de_Db = 0.0022, 100.0, 100.0, 1e-7
        de_pcap = de_m * potts_guy / 3600.0
        de_fu = chem.f_u_vt if chem._is_f_u_vt else ff.cf_u_vt
        de_fnon = chem.f_non_vt if chem._is_f_non_vt else ff.cf_non_vt
        de_binding = 1 - de_falbumin + de_falbumin/de_fu + de_flipid*de_fnon*k_oct
        de_dfree = 10**(-4.15 - 0.655*math.log10(mw))
        self.de_dde = (de_dfree + de_falbumin*de_Db*(1-de_fu)/de_fu) / de_binding
        self.de_kde = de_kfree * de_binding / de_fnon
        self.de_cfree_cmax = (1 - de_falbumin) / de_binding
        de_kfree_bis = 1.0 / (1.0/(de_pcap*de_S) + 1.0/de_Q)
        self.de_ked_clear = de_kfree_bis / de_binding if skin.in_vivo_vitro == IN_VIVO else 0.0
        de_ked_met = sys_obj.kde_met
        de_hde = po.de_thickness / 10000.0
        if self.de_ked_clear > 0:
            de_perm_clear = max(3600*self.de_kde*math.sqrt(self.de_ked_clear*self.de_dde),
                                3600*self.de_dde*self.de_kde/de_hde)
        else:
            de_perm_clear = 3600*self.de_dde*self.de_kde/de_hde
        if self.de_ked_clear == 0:
            de_pde = de_perm_clear
        elif de_ked_met == 0:
            de_pde = de_perm_clear
        else:
            de_pde = (3600*self.de_kde*math.sqrt(self.de_dde*(self.de_ked_clear+de_ked_met))
                      * self.de_ked_clear / (self.de_ked_clear+de_ked_met))
        self.de_omega = 1.0 / de_pde

        # ── firstVE ──
        self.ve_ded = self.de_dde
        self.ve_kedw = self.de_kde
        self.ve_ked_clear = 0.0
        ve_ked_met = sys_obj.kve_met
        ve_hed = po.ve_thickness / 10000.0
        if ve_ked_met > 0:
            ve_perm = (3600*self.ve_kedw*math.sqrt(ve_ked_met*self.ve_ded) /
                       math.sinh(ve_hed*math.sqrt(ve_ked_met/self.ve_ded)))
        else:
            ve_perm = 3600*self.ve_ded*self.ve_kedw / ve_hed
        self.ve_omega = 1.0 / ve_perm

        # ── secondSC ──
        if exp_kp != 0:
            sc_psc_h = 1.0 / (1.0/exp_kp - self.ve_omega - self.de_omega)
            sc_psc_ph = sc_psc_h * psc_scale
            if skin.hydration == FULLY_HYDRATED: sc_sel = sc_psc_h
            else: sc_sel = sc_psc_ph
        else:
            sc_sel = jmn
        sel_ksc_w = psc_obj.ksc_final
        best_sw = sw_obj.solubility / 1000.0
        sc_csat = sel_ksc_w * best_sw
        self.best_csat = min(sc_csat, density)
        self.sc_msat = self.best_csat * (po.deposition / 10000.0) * 1e6
        self.est_ksc_w = self.best_csat / best_sw if best_sw != 0 else 0
        sc_hsc = po.sc_thickness / 10000.0
        self.sc_base_D = sc_hsc * sc_sel / self.est_ksc_w / 3600 if self.est_ksc_w != 0 else 0
        self.sc_ksc_clear = 0.0
        sc_ksc_met = sys_obj.ksc_met
        if sc_ksc_met > 0:
            sc_perm_met = (3600*self.est_ksc_w*math.sqrt(sc_ksc_met*self.sc_base_D) /
                           math.sinh(sc_hsc*math.sqrt(sc_ksc_met/self.sc_base_D)))
        else:
            sc_perm_met = 3600*self.sc_base_D*self.est_ksc_w / sc_hsc
        self.sc_Dsc0 = self.sc_base_D
        self.sc_resistance = 1.0 / sc_perm_met

        # ── finalVE (compute effective kp) ──
        total_res = self.sc_resistance + self.ve_omega + self.de_omega
        self.effective_kp = 1.0 / total_res
        self.sc_dsat_d0 = MINUS_BIG
        self.sc_ctrans = MINUS_BIG
        self.sc_m = MINUS_BIG


class Evap:
    """Evaporation rate constant (Evap.java)."""
    def __init__(self, chem: ChemDat, env: EnvOpt, skin: SkinProp,
                 ro: Ro, sc: SC, po: ProgramOpt):
        pvp = chem.vapour_pressure / TORR_TO_PA  # convert Pa → torr
        ts_K = env.temperature + 273.15
        u_wind = env.wind_speed
        self.w_wind = 44.0 * u_wind**0.78
        mw = chem.mw
        kgp = 6320.0 * u_wind**0.78 / mw**(1/3)
        self.hsc2_dsc0 = (skin.sc_thickness/10000)**2 / sc.sc_Dsc0 / 3600 if sc.sc_Dsc0 > 0 else 0
        prog_R = po.R
        kevap_ro = kgp * pvp * mw / 0.76 / prog_R / (ts_K * 1e6)
        self.kevap = kevap_ro / ro.eff_density if ro.eff_density > 0 else 0
