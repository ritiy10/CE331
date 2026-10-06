import numpy as np

# 1. Physical constants and hit condition derivation
v_b = 250.0  # m/s
v_s = 340.0  # m/s
t_react = 0.010  # seconds

D_crit = (t_react * v_b * v_s) / (v_s - v_b)  # 85/9 m

# 2. Field observations
RS_val = 30.000
sigma_RS_single = 0.002 + 2e-6 * RS_val  # 2mm + 2ppm
sigma_RS = sigma_RS_single / np.sqrt(2)

deg2rad = np.pi / 180.0
arcsec2rad = np.pi / (180.0 * 3600.0)

A_G_R = (69 + 25 / 60 + 53.75 / 3600) * deg2rad
A_B_R = (58 + 31 / 60 + 16.90 / 3600) * deg2rad
A_G_S = (55 + 46 / 60 + 26.47 / 3600) * deg2rad
A_B_S = (71 + 20 / 60 + 15.35 / 3600) * deg2rad

v_G_R = (2 + 49 / 60 + 43.50 / 3600) * deg2rad
v_B_R = (5 + 23 / 60 + 59.69 / 3600) * deg2rad
v_G_S = (3 + 39 / 60 + 42.88 / 3600) * deg2rad
v_B_S = (7 + 10 / 60 + 55.56 / 3600) * deg2rad

# Forward calculations
gamma1 = np.pi - A_G_R - A_G_S
RG = RS_val * np.sin(A_G_S) / np.sin(gamma1)
SG = RS_val * np.sin(A_G_R) / np.sin(gamma1)

gamma2 = np.pi - A_B_R - A_B_S
RB = RS_val * np.sin(A_B_S) / np.sin(gamma2)
SB = RS_val * np.sin(A_B_R) / np.sin(gamma2)

dH_R = RB * np.tan(v_B_R) - RG * np.tan(v_G_R)
dH_S = SB * np.tan(v_B_S) - SG * np.tan(v_G_S)
dH = (dH_R + dH_S) / 2.0

BG_R = np.sqrt(
    RB**2 + RG**2 - 2 * RB * RG * np.cos(abs(A_B_R - A_G_R))
)
BG_S = np.sqrt(
    SB**2 + SG**2 - 2 * SB * SG * np.cos(abs(A_B_S - A_G_S))
)
BG = (BG_R + BG_S) / 2.0

D = np.sqrt(BG**2 + dH**2)

# Constraint: D + 2*sigma_D <= D_crit
max_sigma_D = (D_crit - D) / 2.0


def compute_sigma_D(sigma_theta):
    X0 = np.array(
        [
            RS_val,
            A_G_R,
            A_B_R,
            A_G_S,
            A_B_S,
            v_G_R,
            v_B_R,
            v_G_S,
            v_B_S,
        ]
    )
    cov_X = np.diag([sigma_RS**2] + [sigma_theta**2] * 8)

    def f(X):
        rs, a_gr, a_br, a_gs, a_bs, v_gr, v_br, v_gs, v_bs = X
        g1 = np.pi - a_gr - a_gs
        rg = rs * np.sin(a_gs) / np.sin(g1)

        g2 = np.pi - a_br - a_bs
        rb = rs * np.sin(a_bs) / np.sin(g2)
        sb = rs * np.sin(a_br) / np.sin(g2)

        dh_r = rb * np.tan(v_br) - rg * np.tan(v_gr)
        dh_s = sb * np.tan(v_bs) - (rs * np.sin(a_gr) / np.sin(g1)) * np.tan(v_gs)
        dh = (dh_r + dh_s) / 2.0

        bg_r = np.sqrt(rb**2 + rg**2 - 2 * rb * rg * np.cos(abs(a_br - a_gr)))
        bg_s = np.sqrt(sb**2 + (rs * np.sin(a_gr) / np.sin(g1))**2 - 2 * sb * (rs * np.sin(a_gr) / np.sin(g1)) * np.cos(abs(a_bs - a_gs)))
        bg = (bg_r + bg_s) / 2.0

        return np.sqrt(bg**2 + dh**2)

    eps = 1e-8
    J = np.zeros(len(X0))
    for i in range(len(X0)):
        X_plus = X0.copy()
        X_plus[i] += eps
        X_minus = X0.copy()
        X_minus[i] -= eps
        J[i] = (f(X_plus) - f(X_minus)) / (2 * eps)

    return np.sqrt(J @ cov_X @ J)


# Root finding
low, high = 0.0, 100.0 * arcsec2rad
for _ in range(100):
    mid = (low + high) / 2.0
    if compute_sigma_D(mid) > max_sigma_D:
        high = mid
    else:
        low = mid

sigma_theta_req_arcsec = mid / arcsec2rad
print(f"Minimum required angular accuracy: {sigma_theta_req_arcsec:.2f} arcseconds")
