"""
Orchestrator: wires all physics classes together like OutputDat.SetupSimulation(),
then calls the transient solver. This is the single entry point.
"""
from .chem_data import ChemDat
from .skin_data import SkinProp, EnvOpt, DosageDat, VehicleDat, OutputParameters
from .m2hp import M2H, M2P
from .physics import Psc, Ro, SW, FuFnom, Sys, Kv, ProgramOpt, SC, Evap
from .solver import run_simulation, SimResult


def setup_and_run(chem: ChemDat,
                  skin: SkinProp = None,
                  env: EnvOpt = None,
                  dose: DosageDat = None,
                  nvv: VehicleDat = None,
                  vv: VehicleDat = None,
                  out_par: OutputParameters = None) -> SimResult:
    """
    Full simulation pipeline — mirrors OutputDat.SetupSimulation() + StartSimulation().
    Returns SimResult with time-series and summary statistics.
    """
    skin = skin or SkinProp()
    env = env or EnvOpt()
    dose = dose or DosageDat()
    nvv = nvv or VehicleDat(volatile_type=0)
    vv = vv or VehicleDat(volatile_type=1)
    out_par = out_par or OutputParameters()

    # ── Build object graph (same order as Java) ──
    fu_fnom = FuFnom(chem, skin, nvv)
    ro = Ro(chem, env)
    m2h = M2H(chem)
    m2p = M2P(chem)
    psc = Psc(skin, m2p, m2h)
    sw = SW(chem, env, ro)
    sys_a = Sys(skin, dose)
    kv = Kv(chem, ro, nvv, dose, env, sw, fu_fnom, psc)
    po = ProgramOpt(sys_a, env, skin, dose, out_par, vv)
    sc = SC(chem, psc, sw, ro, skin, fu_fnom, sys_a, po, kv)
    evap = Evap(chem, env, skin, ro, sc, po)

    # ── Package parameters for solver ──
    hsc_cm = po.sc_thickness / 10000.0
    hed_cm = po.ve_thickness / 10000.0
    hde_cm = po.de_thickness / 10000.0

    trans_params = dict(
        Dose=dose.permeant_amount_applied,
        rho=ro.eff_density,
        DoseV=kv.nvv_amount_applied,
        rhoV=nvv.density if nvv._is_density else 0.0,
        Kv=kv.k_v_w,
        Msat=sc.sc_msat,
        kevap=evap.kevap,
        Dsc0=sc.sc_Dsc0,
        F=po.deposition_fraction,
        Ksc=sc.est_ksc_w,
        Ded=sc.ve_ded,
        Ked=sc.ve_kedw,
        Dde=sc.de_dde,
        Kde=sc.de_kde,
        hsc=hsc_cm,
        hed=hed_cm,
        hde=hde_cm,
        finish_time=out_par.max_duration,
        kloss_sc=sys_a.ksc_met,
        kloss_ed=sys_a.kve_met,
        kloss_de=sys_a.kde_met,
    )
    sim_opts = dict(
        N=po.num_sublayers,
        method=po.method,
        h0=po.init_step_size,
        hmin=po.minimum_step_size,
        hmax=po.maximum_step_size,
        dymax=po.dy_max,
        min_output_step=out_par.min_output_step_size,
    )

    result = run_simulation(trans_params, sim_opts)

    # Attach intermediate values for display
    result.kp_eff = sc.effective_kp
    result.Dsc0 = sc.sc_Dsc0
    result.Ksc_w = sc.est_ksc_w
    result.Msat = sc.sc_msat
    result.kevap = evap.kevap
    result.density = ro.eff_density
    result.sw = sw.solubility
    result.sc_resistance_pct = 100 * sc.sc_resistance / (sc.sc_resistance + sc.ve_omega + sc.de_omega)
    result.ve_resistance_pct = 100 * sc.ve_omega / (sc.sc_resistance + sc.ve_omega + sc.de_omega)
    result.de_resistance_pct = 100 * sc.de_omega / (sc.sc_resistance + sc.ve_omega + sc.de_omega)
    return result


def get_steady_state_kp(chem: ChemDat, skin: SkinProp = None,
                        env: EnvOpt = None) -> dict:
    """Quick steady-state kp calculation without running transient solver."""
    skin = skin or SkinProp()
    env = env or EnvOpt()
    dose = DosageDat()
    dose.set_permeant_amount(100.0)
    nvv = VehicleDat(volatile_type=0)
    vv = VehicleDat(volatile_type=1)
    out_par = OutputParameters()

    fu_fnom = FuFnom(chem, skin, nvv)
    ro = Ro(chem, env)
    m2h = M2H(chem)
    m2p = M2P(chem)
    psc_obj = Psc(skin, m2p, m2h)
    sw = SW(chem, env, ro)
    sys_a = Sys(skin, dose)
    kv = Kv(chem, ro, nvv, dose, env, sw, fu_fnom, psc_obj)
    po = ProgramOpt(sys_a, env, skin, dose, out_par, vv)
    sc = SC(chem, psc_obj, sw, ro, skin, fu_fnom, sys_a, po, kv)

    total_R = sc.sc_resistance + sc.ve_omega + sc.de_omega
    return {
        'kp_eff': 1.0 / total_R,
        'Dsc0': sc.sc_Dsc0,
        'Ksc_w': sc.est_ksc_w,
        'Csat': sc.best_csat,
        'Msat': sc.sc_msat,
        'density': ro.eff_density,
        'sw_mg_L': sw.solubility,
        'sc_resist_pct': 100 * sc.sc_resistance / total_R,
        've_resist_pct': 100 * sc.ve_omega / total_R,
        'de_resist_pct': 100 * sc.de_omega / total_R,
        'psc_final': psc_obj.psc_final,
        'ksc_final': psc_obj.ksc_final,
        'log10_s_fh': m2h.log10_s,
        'log10_r_fh': m2h.log10_r,
    }
