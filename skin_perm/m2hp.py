"""
M2H.java, M2P.java → Python.  The micro-scale SC permeability models.
Fully-hydrated (M2H) and partially-hydrated (M2P) brick-and-mortar.
Line-by-line from Java.
"""
import math
from .chem_data import ChemDat


class M2H:
    """Fully-hydrated SC permeability (M2H.java)."""

    # Model-6 polynomial coefficients (25 terms, 4th-order bivariate in log10_s, log10_r)
    _mod6a = [-4.48116, 0.823855, -0.200615, 0.0874482, -0.00969133,
              0.00388257, 7.26031e-4, -0.0163362, 0.0115681, -0.00217086,
              8.34108e-5, -9.56469e-4, -7.16288e-4, 0.00195218, -7.07896e-4,
              -9.4815e-4, -8.09313e-4, 0.0037365, -0.00169451, 9.24718e-5,
              -2.09213e-4, -1.42974e-4, 8.35781e-4, -4.33507e-4, 4.3709e-5]
    _mod6_hole = [-8.9307385, 8.9434534, -4.8497856, 1.0799239, -0.0806654,
                  8.2525913, -14.7285415, 8.1655564, -1.6616479, 0.111266,
                  -4.3162569, 7.4620289, -3.9072178, 0.7232071, -0.0427657,
                  0.8782842, -1.4677591, 0.7222751, -0.1204092, 0.0059583,
                  -0.0622331, 0.1005588, -0.0463148, 0.0067971, -2.461e-4]
    _mod6e = [0,0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 4,4,4,4,4]
    _mod6f = [0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4]

    def __init__(self, chem: ChemDat):
        va = chem.get_va()
        hsch = 0.0043365  # SC thickness in cm (fully hydrated)
        mw = chem.mw
        logk_ow = chem.logP
        k_ow = 10.0 ** logk_ow

        klip_w = 0.43 * k_ow ** 0.81
        dlip = 1.24e-7 * (100.0 / mw) ** 2.43 + 2.34e-9
        log10ktrans = -0.725 - 0.792 * mw ** (1/3)
        ktrans = 10.0 ** log10ktrans

        if va < 445.2:
            daq = 1.92e-4 / va ** 0.6
            a_s = 0.145 * va ** 0.6
        else:
            daq = 3.78e-5 * va ** (1/3)
            a_s = 0.735 * va ** (1/3)

        lam = a_s / 35.0
        ff1 = 0.1928 * (1 + lam) ** 2
        k_cor_w_free = 1.0 - ff1
        dcor_free = daq * (1 - ff1) * (0.9999 - 1.2762*lam + 0.0718*lam**2 + 0.1195*lam**3)

        sigma = klip_w * dlip / k_cor_w_free / dcor_free
        r = ktrans / 0.132536 / dlip
        self.log10_s = math.log10(sigma)
        self.log10_r = math.log10(r)

        self.ksc_w = 0.014 * k_ow**0.81 + 0.782 + 1.381 * k_ow**0.27

        # Model-6 Psc
        log10_psc = self._model6(self._mod6a)
        psc_sec = 10**log10_psc * klip_w * dlip / hsch
        self.psc_hour = psc_sec * 3600.0

        # Model-6 hole pathway
        log10_psc_hole = self._model6(self._mod6_hole)
        psc_sec_hole = 10**log10_psc_hole * klip_w * dlip / hsch
        self.psc_hour_hole = psc_sec_hole * 3600.0

        # Composite Psc (B43)
        K5, K6, K7, K8, K9 = -4.49523, 0.854964, -0.120539, -0.0153812, 0.00339236
        log_phat = K5 + self.log10_r*(K6 + self.log10_r*(K7 + self.log10_r*(K8 + self.log10_r*K9)))
        self.pcompscw = 3600 * 10**log_phat * klip_w * dlip / hsch

        # Low-r asymptote (B47)
        phat_lor = 3.803e-5 * r
        self.pcompscw_a3b_lor = 3600 * phat_lor * klip_w * dlip / hsch

        # High-r asymptote (B50)
        self.pcompscw_a3b_hir = 2.88 * klip_w * dlip / hsch

        # Analytic (B56)
        p_scw_an = 1.0 / (0.9685 * sigma + 182700.0 / r)
        self.pscw_analytic_hour = 3600 * p_scw_an * klip_w * dlip / hsch

    def _model6(self, coeffs):
        s, r = self.log10_s, self.log10_r
        return sum(coeffs[i] * s**self._mod6e[i] * r**self._mod6f[i] for i in range(25))


class M2P:
    """Partially-hydrated SC permeability (M2P.java)."""

    _mod6a = [-4.95725, 0.838263, -0.207133, 0.0837613, -0.00861439,
              0.00187015, -1.05288e-4, -0.00849198, 0.00666713, -0.00119412,
              7.45813e-4, -1.8624e-4, -0.00360234, 0.00343394, -7.68073e-4,
              1.82904e-4, -8.28876e-5, -5.26007e-4, 7.1352e-4, -2.13017e-4,
              2.83838e-5, -1.98892e-5, -6.05354e-6, 4.84676e-5, -2.1566e-5]
    _mod6_hole = [-4.3591261, -0.899076, 1.4332295, -0.475971, 0.0500664,
                  -1.7486129, 4.6689678, -4.0774937, 1.3157635, -0.1333001,
                  1.7949242, -4.2658909, 3.3771508, -1.0093969, 0.0964929,
                  -0.5711852, 1.2859666, -0.9628191, 0.2729855, -0.0251096,
                  0.0552116, -0.1206319, 0.0873323, -0.0239349, 0.0021466]
    _mod6e = [0,0,0,0,0, 1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 4,4,4,4,4]
    _mod6f = [0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4, 0,1,2,3,4]

    def __init__(self, chem: ChemDat):
        va = chem.get_va()
        hschp = 0.0013365  # SC thickness in cm (partially hydrated)
        mw = chem.mw
        logk_ow = chem.logP
        k_ow = 10.0 ** logk_ow

        hlat = 3.0; htrans = 3.0
        klip_w = 0.43 * k_ow ** 0.81
        dlip = (1.24e-7 * (100.0/mw)**2.43 + 2.34e-9) / hlat
        log10ktrans = -0.725 - 0.792 * mw**(1/3) - math.log10(htrans)
        ktrans = 10.0 ** log10ktrans

        if va < 445.2:
            daq = 1.92e-4 / va**0.6
            a_s = 0.145 * va**0.6
        else:
            daq = 3.78e-5 * va**(1/3)
            a_s = 0.735 * va**(1/3)

        lam = a_s / 35.0
        ff1 = 0.6044 * (1 + lam)**2
        k_cor_w_free = 1.0 - ff1
        dcor_free = daq * (1-ff1) * (1.0001 - 2.4497*lam + 1.141*lam**2 + 0.5432*lam**3)

        sigma = klip_w * dlip / k_cor_w_free / dcor_free
        r = ktrans / 0.141916 / dlip
        self.log10_s = math.log10(sigma)
        self.log10_r = math.log10(r)

        self.ksc_w = 0.04 * k_ow**0.81 + 0.359 + 4.057 * k_ow**0.27

        log10_psc = self._model6(self._mod6a)
        psc_sec = 10**log10_psc * klip_w * dlip / hschp
        self.psc_hour = psc_sec * 3600.0

        log10_psc_hole = self._model6(self._mod6_hole)
        psc_sec_hole = 10**log10_psc_hole * klip_w * dlip / hschp
        self.psc_hour_hole = psc_sec_hole * 3600.0

        K5, K6, K7, K8, K9 = -4.97447, 0.86578, -0.113774, -0.016918, 0.00344476
        log_phat = K5 + self.log10_r*(K6 + self.log10_r*(K7 + self.log10_r*(K8 + self.log10_r*K9)))
        self.pcompscw = 3600 * 10**log_phat * klip_w * dlip / hschp

        phat_lor = 1.242e-5 * r
        self.pcompscw_a3b_lor = 3600 * phat_lor * klip_w * dlip / hschp
        self.pcompscw_a3b_hir = 1.21356 * klip_w * dlip / hschp

        p_scw_an = 1.0 / (0.8979 * sigma + 553600.0 / r)
        self.pscw_analytic_hour = 3600 * p_scw_an * klip_w * dlip / hschp

    def _model6(self, coeffs):
        s, r = self.log10_s, self.log10_r
        return sum(coeffs[i] * s**self._mod6e[i] * r**self._mod6f[i] for i in range(25))
