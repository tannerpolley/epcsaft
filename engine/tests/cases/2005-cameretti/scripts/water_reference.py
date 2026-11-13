"""Water saturation reference: Wagner and Pruss (1993), J. Phys. Chem. Ref. Data 22, 783, Eqs. (2.5) and (2.6),
the IAPWS auxiliary equations for the saturation pressure and the saturated-liquid density."""
import math

TC, PC, RHOC = 647.096, 22.064e6, 322.0  # K, Pa, kg/m3
A = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
B = (1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -6.74694450e5)


def saturation_pressure(T):
    tau = 1 - T / TC
    terms = (tau, tau**1.5, tau**3, tau**3.5, tau**4, tau**7.5)
    return PC * math.exp(TC / T * sum(a * t for a, t in zip(A, terms)))


def saturated_liquid_density(T):
    tau = 1 - T / TC
    exponents = (1 / 3, 2 / 3, 5 / 3, 16 / 3, 43 / 3, 110 / 3)
    return RHOC * (1 + sum(b * tau**e for b, e in zip(B, exponents)))


if __name__ == "__main__":
    # Published anchors: 3.1699 kPa and 997.05 kg/m3 at 298.15 K; 101.325 kPa at 373.124 K.
    assert abs(saturation_pressure(298.15) / 3169.9 - 1) < 2e-4
    assert abs(saturated_liquid_density(298.15) / 997.05 - 1) < 1e-4
    assert abs(saturation_pressure(373.124) / 101325 - 1) < 2e-4
