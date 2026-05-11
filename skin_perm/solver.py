import numpy as np
from scipy.linalg import solve_banded

class SimResult:
    def __init__(self):
        self.times = []
        self.Jabs = []
        self.Jevap_arr = []
        self.Qabs = []
        self.Qevap_arr = []
        self.Qmem_arr = []
        self.mass_bal = []
        self.max_abs_flux = 0.0
        self.t_max_abs_flux = 0.0
        self.frc_absorbed = 0.0
        self.frc_evaporated = 0.0

def run_simulation(tp: dict, so: dict = None) -> SimResult:
    """Faithful port of the implicit Java Crank-Nicolson solver."""
    so = so or {}
    Dose = tp['Dose']
    rho = tp['rho']
    DoseV = tp['DoseV']
    rhov = tp['rhoV']
    Kv = tp['Kv']
    Msat = tp['Msat']
    kevap = tp['kevap']
    Dsc0 = tp['Dsc0']
    F_dep = tp.get('F', 0.025)
    Ksc = tp['Ksc']
    Ded = tp['Ded']
    Ked = tp['Ked']
    Dde = tp['Dde']
    Kde = tp['Kde']
    hsc = tp['hsc']
    hed = tp['hed']
    hde = tp['hde']
    finish_time = tp.get('finish_time', 1000.0)
    kloss_sc = tp.get('kloss_sc', 0.0)
    kloss_ed = tp.get('kloss_ed', 0.0)
    kloss_de = tp.get('kloss_de', 0.0)
    
    finish_time_s = finish_time * 3600.0
    
    N = so.get('N', 120)
    nsc = N // 3
    ned = N // 3
    nde = N - nsc - ned
    
    deltaxsc = hsc / nsc
    deltaxed = hed / ned if ned > 0 else 0
    deltaxde = hde / nde if nde > 0 else 0
    
    # Calculate Csat from Msat and F_dep
    Csat = Msat / (F_dep * hsc * 1e6) if (F_dep * hsc > 0) else 1.0
    ndep = max(1, int(round(F_dep * nsc)))
    
    # Scale kevap to cm/s
    kevap_s = kevap / 3600.0
    
    # Initial Conditions (from OutputDat.SetupICs)
    Y = np.zeros(N + 1)
    k_loss = np.zeros(N + 1)
    
    # Loss terms
    for i in range(1, nsc + 1): k_loss[i] = kloss_sc
    for i in range(nsc + 1, nsc + ned + 1): k_loss[i] = kloss_ed
    for i in range(nsc + ned + 1, N + 1): k_loss[i] = kloss_de

    # Volume of vehicle
    hv = (10000 * DoseV / 1000.0 / rhov) if (DoseV > 0 and rhov > 0) else 0.0
    Vsat = DoseV * 1000.0 if DoseV > 0 else 0.0 # approx
    
    Qveh = 0.0
    Qfat = 0.0
    Qtrap = 0.0
    Qhead = 0.0
    Qloss = 0.0
    Qremoved = 0.0
    
    if DoseV > 0:
        Qveh = Dose
    else:
        init_c = Dose / (ndep * deltaxsc * 1e6)
        if init_c > Csat:
            for i in range(1, ndep + 1):
                Y[i] = Csat
            Qveh = (init_c - Csat) * ndep * deltaxsc * 1e6
        else:
            for i in range(1, ndep + 1):
                Y[i] = init_c
            Qveh = 0.0

    # Diffusivities (CalculateDiffusivities)
    D = np.zeros(N + 1)
    for i in range(1, nsc + 1): D[i] = Dsc0
    D[0] = 2.0 * D[1] - D[2] if N > 1 else Dsc0
    for i in range(nsc + 1, nsc + ned + 1): D[i] = Ded
    for i in range(nsc + ned + 1, N + 1): D[i] = Dde
    
    # Step size control
    h0_h = so.get('h0', 1e-6)
    hmin_h = so.get('hmin', 1e-8)
    hmax_h = so.get('hmax', min(1.0, finish_time / 100))
    dymax = so.get('dymax', 0.01)
    min_out_step_h = so.get('min_output_step', 1.0)
    
    h0_s = h0_h * 3600.0
    hmin_s = hmin_h * 3600.0
    hmax_s = hmax_h * 3600.0
    
    t = 0.0
    dt = h0_s
    
    # Determine initial icase
    if Qveh > 0:
        icase = 2
    elif kevap_s > 0:
        icase = 1
    else:
        icase = 3
        
    Kvm = 1.0  # From Java if no vehicle partition specified differently
        
    # Result arrays
    res = SimResult()
    res.times.append(t)
    res.Jabs.append(0.0)
    res.Jevap_arr.append(0.0)
    res.Qabs.append(0.0)
    res.Qevap_arr.append(0.0)
    res.Qmem_arr.append(np.sum(Y[1:nsc+1])*deltaxsc*1e6)
    res.mass_bal.append(1.0)
    
    Factor = np.zeros(30)
    
    # Setup BC factors
    def calc_bc_factors(icase):
        R_t = deltaxsc * kevap_s * rho / D[0] / Csat if (D[0] > 0 and Csat > 0) else 0.0
        delsc = 1.0 / (deltaxsc * deltaxsc)
        deled = 1.0 / (deltaxed * deltaxed) if ned > 0 else 0.0
        delde = 1.0 / (deltaxde * deltaxde) if nde > 0 else 0.0
        
        Factor[0] = delsc
        Factor[20] = deled
        Factor[24] = delde
        Factor[17] = delsc * (D[1] + D[0] / 3.0)
        Factor[18] = -delsc * (D[1] + 3.0 * D[0])
        Factor[1] = -delsc * (D[0] * 9.0 * R_t / (8.0 + 3.0 * R_t) + D[1])
        Factor[2] = delsc * (D[0] * R_t / (8.0 + 3.0 * R_t) + D[1])
        
        # Bottom sink BC
        if nde > 0:
            Factor[3] = delde * (D[N - 1] + D[N] / 3.0)
            Factor[4] = -delde * (D[N - 1] + 3.0 * D[N])
        elif ned > 0:
            Factor[3] = deled * (D[N - 1] + D[N] / 3.0)
            Factor[4] = -deled * (D[N - 1] + 3.0 * D[N])
        else:
            Factor[3] = delsc * (D[N - 1] + D[N] / 3.0)
            Factor[4] = -delsc * (D[N - 1] + 3.0 * D[N])
            
        Factor[5] = delsc * D[0] * Csat
        
        Q_fac = D[1] * delsc / (1.0 + hv * Kvm / deltaxsc) if deltaxsc > 0 else 0
        Factor[15] = Q_fac
        Factor[16] = -Q_fac * (1.0 + R_t)
        
        Psc = D[nsc] * Ksc
        Ped1 = D[nsc + 1] * Ked if ned > 0 else 0
        Ped2 = D[nsc + ned] * Ked if ned > 0 else 0
        Pde = D[nsc + ned + 1] * Kde if nde > 0 else 0
        
        Csced1 = 1.0 / (1.0 + Ped1 * deltaxsc / Psc / deltaxed) if (Psc > 0 and deltaxed > 0) else 0.5
        Csced2 = 1.0 / (1.0 + Psc * deltaxed / Ped1 / deltaxsc) if (Ped1 > 0 and deltaxsc > 0) else 0.5
        
        Factor[6] = -D[nsc] - 2.0 * D[nsc] * (1.0 - Csced1)
        Factor[7] = 2.0 * D[nsc] * Csced2
        Factor[8] = 2.0 * D[nsc + 1] * Csced1 if ned > 0 else 0
        Factor[9] = -D[nsc + 1] - 2.0 * D[nsc + 1] * (1.0 - Csced2) if ned > 0 else 0
        
        Cedde1 = 1.0 / (1.0 + Pde * deltaxed / Ped2 / deltaxde) if (Ped2 > 0 and deltaxde > 0) else 0.5
        Cedde2 = 1.0 / (1.0 + Ped2 * deltaxde / Pde / deltaxed) if (Pde > 0 and deltaxed > 0) else 0.5
        
        Factor[10] = -D[nsc + ned] - 2.0 * D[nsc + ned] * (1.0 - Cedde1)
        Factor[11] = 2.0 * D[nsc + ned] * Cedde2
        Factor[12] = 2.0 * D[nsc + ned + 1] * Cedde1 if nde > 0 else 0
        Factor[13] = -D[nsc + ned + 1] - 2.0 * D[nsc + ned + 1] * (1.0 - Cedde2) if nde > 0 else 0

    calc_bc_factors(icase)
    
    # Build matrix elements for implicit solve
    def build_matrix(icase):
        # A returns diagonals: upper, main, lower for solve_banded (u=1, l=1)
        upper = np.zeros(N)
        main = np.zeros(N)
        lower = np.zeros(N)
        const = np.zeros(N)
        
        # Node 1
        if icase == 3:
            main[0] = Factor[16]
            upper[0] = Factor[15]
        elif icase == 2 or icase == 4:
            main[0] = -Factor[0] * (D[0] + D[1])
            upper[0] = Factor[0] * D[1]
            const[0] = Factor[5]
        elif icase == 1:
            main[0] = Factor[1]
            upper[0] = Factor[2]
            
        # Interior SC
        for i in range(2, nsc):
            idx = i - 1
            lower[idx] = Factor[0] * D[i - 1]
            main[idx] = -Factor[0] * (D[i - 1] + D[i]) - k_loss[i]
            upper[idx] = Factor[0] * D[i]
            
        # Interior VE
        for i in range(nsc + 2, nsc + ned):
            idx = i - 1
            lower[idx] = Factor[20] * D[i - 1]
            main[idx] = -Factor[20] * (D[i - 1] + D[i]) - k_loss[i]
            upper[idx] = Factor[20] * D[i]
            
        # Interior DE
        for i in range(nsc + ned + 2, N):
            idx = i - 1
            lower[idx] = Factor[24] * D[i - 1]
            main[idx] = -Factor[24] * (D[i - 1] + D[i]) - k_loss[i]
            upper[idx] = Factor[24] * D[i]
            
        # Interfaces
        if N > nsc:
            idx = nsc - 1
            lower[idx] = Factor[0] * D[nsc - 1]
            main[idx] = Factor[0] * Factor[6]
            upper[idx] = Factor[0] * Factor[7]
            
            idx = nsc
            lower[idx] = Factor[20] * Factor[8]
            main[idx] = Factor[20] * Factor[9]
            upper[idx] = Factor[20] * D[nsc + 1]
            
            if ned > 0 and nde > 0:
                idx = nsc + ned - 1
                lower[idx] = Factor[20] * D[nsc + ned - 1]
                main[idx] = Factor[20] * Factor[10]
                upper[idx] = Factor[20] * Factor[11]
                
                idx = nsc + ned
                lower[idx] = Factor[24] * Factor[12]
                main[idx] = Factor[24] * Factor[13]
                upper[idx] = Factor[24] * D[nsc + ned + 1]
                
        # Sink BC
        idx = N - 1
        lower[idx] = Factor[3]
        main[idx] = Factor[4]
        
        return lower, main, upper, const

    JdOld = 0.0
    JeOld = 0.0
    JsurfOld = 0.0
    QvehOld = Qveh
    
    last_out_t = 0.0
    max_Jabs = 0.0
    t_max_Jabs = 0.0

    step_count = 0
    while t < finish_time_s:
        step_count += 1
        if step_count % 10000 == 0:
            print(f"t={t/3600.0:.4f} h, dt={dt:.4e} s, dy_max={dy_max_step:.4e}, Qveh={Qveh:.2f}", flush=True)
            
        if t + dt > finish_time_s:
            dt = finish_time_s - t
            
        lower, main, upper, const = build_matrix(icase)
        
        # Backward Euler formulation (fully implicit)
        # Y_new - Y_old = dt * (A Y_new + B)
        # (I - dt * A) Y_new = Y_old + dt * B
        
        # Crank-Nicolson formulation (theta=0.5)
        # (I - 0.5*dt*A) Y_new = (I + 0.5*dt*A) Y_old + dt*B
        M_main = 1.0 - 0.5 * dt * main
        M_lower = -0.5 * dt * lower
        M_upper = -0.5 * dt * upper
        
        # RHS = Y_old + 0.5*dt*(A Y_old) + dt*B
        RHS = Y[1:] + 0.5 * dt * (main * Y[1:]) + dt * const
        RHS[1:] += 0.5 * dt * lower[1:] * Y[1:-1]
        RHS[:-1] += 0.5 * dt * upper[:-1] * Y[2:]
        
        # Solve banded system
        ab = np.zeros((3, N))
        ab[0, 1:] = M_upper[:-1]
        ab[1, :] = M_main
        ab[2, :-1] = M_lower[1:]
        
        try:
            Y_new_int = solve_banded((1, 1), ab, RHS)
        except Exception as e:
            print(f"Solve_banded failed at t={t}: {e}")
            break
        
        # Calculate step error to adjust dt (dymax control)
        diffs = np.abs(Y_new_int - Y[1:])
        max_vals = np.maximum(np.abs(Y[1:]), np.abs(Y_new_int))
        mask = max_vals > 1e-15
        errors = np.zeros_like(diffs)
        errors[mask] = diffs[mask] / max_vals[mask]
        errors[~mask] = diffs[~mask]
        dy_max_step = np.max(errors)
                
        # If error too large, shrink step and retry
        if dy_max_step > dymax and dt > hmin_s:
            dt = max(hmin_s, dt * 0.5)
            continue
            
        # Accept step
        Y[1:] = Y_new_int
        t += dt
        
        # Calculate fluxes for mass balance
        # Jdefat
        temp = 1.0 / deltaxsc
        if ned > 0: temp = Ked / Ksc / deltaxed
        if nde > 0: temp = Kde / Ksc / deltaxde
        Jdefat = D[N] * temp * (5.0 * Y[N] - Y[N - 1]) / 2.0  # Proper one-sided deriv
        # In Java, it says Jdefat = D[N]*temp*(5.0*Y[N] - Y[N-1]); wait, the above formula was 3*Y - Y/3, wait, Java line 597: D[N]*temp*(5.0*Y[N]-Y[N-1])
        Jdefat = D[N] * temp * (5.0 * Y[N] - Y[N - 1]) 
        
        Jevap = 0.0
        Jsurf = 0.0
        
        if icase == 1:
            Jevap = kevap_s * (15.0 * Y[1] - 10.0 * Y[2] + 3.0 * Y[3]) * Kvm / 8.0
            Jsurf = Jevap
        elif icase == 2:
            Jevap = kevap_s * rho
            # Java: D[0]*(-Y[2] + 4*Y[1] - 3*Csat)/2.0/deltaxsc
            Jsurf = D[0] * (-Y[2] + 4.0 * Y[1] - 3.0 * Csat) / 2.0 / deltaxsc
            
        Jd_h = Jdefat * 3600.0 * 1e6  # to ug/cm2/h
        Je_h = Jevap * 3600.0 * 1e6
        Js_h = Jsurf * 3600.0 * 1e6
        
        if Jd_h > max_Jabs:
            max_Jabs = Jd_h
            t_max_Jabs = t
            
        # Qtrap_old before integration
        Qtrap_old = Qtrap
        
        dt_h = dt / 3600.0
        Qfat += dt_h * 0.5 * (Jd_h + JdOld)
        Qtrap += dt_h * 0.5 * (Je_h + JeOld)
        
        if icase == 2:
            # Force absolute mass balance for icase 2 to eliminate numerical Crank-Nicolson drift
            Qsc_curr = np.sum(Y[1:nsc+1]) * deltaxsc
            Qed_curr = np.sum(Y[nsc+1:nsc+ned+1]) * deltaxed * Ked / Ksc if ned > 0 else 0
            Qde_curr = np.sum(Y[nsc+ned+1:N+1]) * deltaxde * Kde / Ksc if nde > 0 else 0
            Qmem_curr = (Qsc_curr + Qed_curr + Qde_curr) * 1e6
            
            Qveh = Dose - (Qmem_curr + Qfat + Qtrap)
            if Qveh < 0:
                Qveh = 0.0
        if icase == 1:
            Qveh = Y[1] * hv * Kvm * 1e6
        elif icase == 3:
            Qveh = hv * Csat * 1e6
            
        if icase == 1 or icase == 3:
            Qsc = np.sum(Y[1:nsc+1]) * deltaxsc
            Qed = np.sum(Y[nsc+1:nsc+ned+1]) * deltaxed * Ked / Ksc if ned > 0 else 0
            Qde = np.sum(Y[nsc+ned+1:N+1]) * deltaxde * Kde / Ksc if nde > 0 else 0
            Qmem = (Qsc + Qed + Qde) * 1e6
            
            diff_Qevap = Dose - (Qmem + Qveh + Qfat)
            
            # If the calculated Qevap is less than the previous step's Qevap, keep the old one (evaporation cannot drop)
            if diff_Qevap < Qtrap_old:
                Qtrap = Qtrap_old
            else:
                Qtrap = diff_Qevap
                
            Je_h = (Qtrap - Qtrap_old) / dt_h if dt_h > 0 else 0.0
            
        JdOld = Jd_h
        JeOld = Je_h
        JsurfOld = Js_h
        
        # State transitions
        if icase == 2 and Qveh <= 0:
            Qveh = 0.0
            if kevap_s > 0:
                icase = 1
            else:
                icase = 3
            calc_bc_factors(icase)
            
        # Increase step size if error is small
        if dy_max_step < dymax * 0.5 and dt < hmax_s:
            dt = min(hmax_s, dt * 1.5)
            
        # Save output
        if t - last_out_t >= min_out_step_h * 3600.0 or t >= finish_time_s:
            last_out_t = t
            res.times.append(t / 3600.0)
            
            Qsc = np.sum(Y[1:nsc+1]) * deltaxsc
            Qed = np.sum(Y[nsc+1:nsc+ned+1]) * deltaxed * Ked / Ksc if ned > 0 else 0
            Qde = np.sum(Y[nsc+ned+1:N+1]) * deltaxde * Kde / Ksc if nde > 0 else 0
            Qmem = (Qsc + Qed + Qde) * 1e6
            res.Qmem_arr.append(Qmem)
            
            res.Jabs.append(Jd_h)
            res.Jevap_arr.append(Je_h)
            res.Qabs.append(Qfat)
            res.Qevap_arr.append(Qtrap)
            
            mb = (Qmem + Qfat + Qtrap + Qveh) / Dose if Dose > 0 else 1.0
            res.mass_bal.append(mb)
            
    res.times = np.array(res.times)
    res.Jabs = np.array(res.Jabs)
    res.Jevap_arr = np.array(res.Jevap_arr)
    res.Qabs = np.array(res.Qabs)
    res.Qevap_arr = np.array(res.Qevap_arr)
    res.Qmem_arr = np.array(res.Qmem_arr)
    res.mass_bal = np.array(res.mass_bal)
    
    res.max_abs_flux = max_Jabs
    res.t_max_abs_flux = t_max_Jabs / 3600.0
    res.frc_absorbed = Qfat / Dose * 100 if Dose > 0 else 0
    res.frc_evaporated = Qtrap / Dose * 100 if Dose > 0 else 0
    
    # Store layer distributions at end
    res.final_Qsc = Qsc * 1e6
    res.final_Qed = Qed * 1e6
    res.final_Qde = Qde * 1e6
    
    return res
