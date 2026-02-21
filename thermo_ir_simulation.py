"""
1D Transient Heat Conduction + IR Camera Simulation + Emissivity Correction
一维瞬态导热 + 红外相机测量模拟 + 发射率修正

Teaching-oriented numerical simulation.
"""

import numpy as np
import matplotlib.pyplot as plt


# PARAMETERS (easy to modify)

# Geometry
L1 = 0.05          # Length of material A [m]
L2 = 0.05          # Length of material B [m]
N = 51             # Number of grid points (fewer = faster, 51 is enough for teaching)

# Material A: k [W/(m·K)], rho [kg/m³], cp [J/(kg·K)]
k_A, rho_A, cp_A = 50.0, 7800.0, 500.0
# Material B
k_B, rho_B, cp_B = 0.2, 1000.0, 2000.0

# Boundary and initial
Thot = 400.0       # Hot end temperature [K]
Tcold = 300.0      # Cold end temperature [K]
T0 = 350.0         # Initial uniform temperature [K]

# Time
t_end = 100.0      # Simulation end time [s]
num_time_steps = 1500   # Fewer steps = faster (dt must still satisfy stability)
plot_times = [5.0, 20.0, 50.0, 100.0]  # Times at which to plot profiles

# Emissivity
epsilon_real_A = 0.85   # Real emissivity material A
epsilon_real_B = 0.95   # Real emissivity material B
epsilon_set = 0.90      # Emissivity assumed by camera (wrong)

# Physics
sigma_sb = 5.670374419e-8   # Stefan–Boltzmann constant [W/(m²·K⁴)]


def setup_grid(L1, L2, N):
    """
    Create 1D spatial grid for two-layer domain [0, L1+L2].
    建立一维空间网格 [0, L1+L2]。
    """
    L_total = L1 + L2
    x = np.linspace(0, L_total, N)
    dx = x[1] - x[0]
    return x, dx, L_total


def assign_material_properties(x, L1, L2,
                               k_A, rho_A, cp_A, k_B, rho_B, cp_B,
                               epsilon_real_A, epsilon_real_B):
    """
    Assign k, rho*cp, and epsilon_real at each grid point.
    Interface: properties change at x = L1.
    在每个网格点赋予 k、rho*cp 和 epsilon_real；界面在 x=L1。
    """
    N = len(x)
    k = np.where(x <= L1, k_A, k_B)
    rho_cp = np.where(x <= L1, rho_A * cp_A, rho_B * cp_B)
    alpha = k / rho_cp   # Thermal diffusivity 热扩散率
    epsilon_real = np.where(x <= L1, epsilon_real_A, epsilon_real_B)
    return k, rho_cp, alpha, epsilon_real


def solve_heat_conduction(x, dx, k, rho_cp, Thot, Tcold, T0, t_end, num_time_steps):
    """
    Solve 1D transient heat conduction using explicit FDM.
    Equation: rho*cp*dT/dt = d/dx(k*dT/dx).
    Interface handled by harmonic-mean thermal conductance.
    显式有限差分求解一维瞬态导热，界面处 k 用调和平均。
    """
    N = len(x)
    dt = t_end / num_time_steps

    # Harmonic mean of k at cell faces (for flux) 界面导热系数调和平均
    k_half = np.zeros(N + 1)
    k_half[0], k_half[-1] = k[0], k[-1]
    for i in range(1, N):
        if k[i] != 0 and k[i - 1] != 0:
            k_half[i] = 2 * k[i] * k[i - 1] / (k[i] + k[i - 1])
        else:
            k_half[i] = 0.5 * (k[i] + k[i - 1])

    T = np.full(N, T0, dtype=float)
    T[0], T[-1] = Thot, Tcold   # Boundary conditions 边界条件

    # Time history for selected indices (all points, all steps for colormap)
    T_history = np.zeros((N, num_time_steps + 1))
    T_history[:, 0] = T.copy()

    # Stability: dt <= min(dx^2/(2*alpha)); use 0.5 factor for safety
    alpha = k / rho_cp
    dt_stable = 0.4 * dx**2 / np.max(alpha)
    if dt > dt_stable:
        raise ValueError(
            f"Unstable: dt={dt:.2e} > dt_stable={dt_stable:.2e}. "
            "Increase num_time_steps or decrease N."
        )

    for n in range(1, num_time_steps + 1):
        t = n * dt
        T_new = T.copy()
        for i in range(1, N - 1):
            # Fluxes at i+1/2 and i-1/2
            F_plus = k_half[i + 1] * (T[i + 1] - T[i]) / dx
            F_minus = k_half[i] * (T[i] - T[i - 1]) / dx
            dT_dt = (F_plus - F_minus) / (dx * rho_cp[i])
            T_new[i] = T[i] + dT_dt * dt
        T_new[0], T_new[-1] = Thot, Tcold
        T = T_new
        T_history[:, n] = T.copy()

    t_array = np.linspace(0, t_end, num_time_steps + 1)
    return T_history, t_array, dt


def compute_radiation_intensity(T_real, epsilon_real, sigma_sb):
    """
    True radiation intensity: I = epsilon_real * sigma * T^4.
    真实辐射强度 I = ε_real * σ * T^4。
    """
    # epsilon_real (N,) broadcasts over time with T_real (N, n_steps)
    eps = np.asarray(epsilon_real).reshape(-1, 1) if T_real.ndim == 2 else epsilon_real
    I = eps * sigma_sb * (T_real ** 4)
    return I


def compute_camera_temperature(I, epsilon_set, sigma_sb):
    """
    Camera infers temperature using wrong emissivity: T_cam = (I/(epsilon_set*sigma))^0.25.
    相机用错误发射率反算温度。
    """
    denom = epsilon_set * sigma_sb
    denom = np.maximum(denom, 1e-20)  # Avoid division by zero
    T_cam = (I / denom) ** 0.25
    return T_cam


def emissivity_correction(T_cam, epsilon_set, epsilon_real):
    """
    Theoretical correction: T_corr = (epsilon_set/epsilon_real)^0.25 * T_cam.
    理论修正公式。
    """
    eps_real = np.asarray(epsilon_real).reshape(-1, 1) if T_cam.ndim == 2 else epsilon_real
    ratio = np.maximum(epsilon_set / np.maximum(eps_real, 1e-10), 1e-10)
    T_corr = (ratio ** 0.25) * T_cam
    return T_corr


def plot_temperature_profiles(x, T_history, t_array, plot_times,
                              T_cam_history, T_corr_history,
                              L1, L2):
    """
    Plot T_real, T_cam, T_corr at several time steps (same axes for comparison).
    在若干时刻绘制 T_real、T_cam、T_corr 对比。
    """
    L_total = L1 + L2
    n_plots = len(plot_times)
    fig, axes = plt.subplots(1, n_plots, figsize=(4 * n_plots, 4), sharey=True)
    if n_plots == 1:
        axes = [axes]

    for idx, t_plot in enumerate(plot_times):
        # Find time index closest to t_plot
        i_t = np.argmin(np.abs(t_array - t_plot))
        t_actual = t_array[i_t]

        ax = axes[idx]
        ax.plot(x, T_history[:, i_t], 'k-', lw=2, label='$T_{real}$')
        ax.plot(x, T_cam_history[:, i_t], 'b--', lw=1.5, label='$T_{cam}$')
        ax.plot(x, T_corr_history[:, i_t], 'r:', lw=1.5, label='$T_{corr}$')
        ax.axvline(L1, color='gray', ls='--', alpha=0.7)
        ax.set_xlabel('$x$ [m]')
        ax.set_title(f'$t$ = {t_actual:.1f} s')
        ax.legend(loc='best', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, L_total)

    axes[0].set_ylabel('Temperature [K]')
    fig.suptitle('Temperature profiles: $T_{real}$, $T_{cam}$, $T_{corr}$', fontsize=12)
    plt.tight_layout()
    return fig


def plot_IR_colormap(x, t_array, T_field, title, L1, L2):
    """
    Pseudo-IR colormap: imshow with cmap='inferno'.
    x horizontal, time vertical (t increasing downward).
    伪红外热图：x 为横轴，时间纵轴（向下增加）。
    """
    L_total = L1 + L2
    fig, ax = plt.subplots(figsize=(6, 4))
    # imshow: first dim is y (time), second is x (space)
    extent = [x[0], x[-1], t_array[-1], t_array[0]]  # t top=0, bottom=t_end
    im = ax.imshow(T_field, aspect='auto', extent=extent, cmap='inferno',
                   interpolation='bilinear')
    ax.axvline(L1, color='cyan', ls='--', alpha=0.8)
    ax.set_xlabel('$x$ [m]')
    ax.set_ylabel('Time [s]')
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label='Temperature [K]')
    plt.tight_layout()
    return fig


def main():
    """Run full pipeline: grid -> properties -> solve -> IR -> correct -> plot."""
    # Grid
    x, dx, L_total = setup_grid(L1, L2, N)
    k, rho_cp, alpha, epsilon_real = assign_material_properties(
        x, L1, L2, k_A, rho_A, cp_A, k_B, rho_B, cp_B,
        epsilon_real_A, epsilon_real_B
    )

    # Solve heat conduction 求解导热
    T_history, t_array, dt = solve_heat_conduction(
        x, dx, k, rho_cp, Thot, Tcold, T0, t_end, num_time_steps
    )
    # T_history shape: (N, num_time_steps+1); need (time, space) for imshow
    T_real_imshow = T_history.T   # (time, space)

    # IR camera simulation: I from T_real, then T_cam from I
    # epsilon_real is (N,) so we broadcast over time
    I_history = compute_radiation_intensity(T_history, epsilon_real, sigma_sb)
    T_cam_history = compute_camera_temperature(I_history, epsilon_set, sigma_sb)
    T_corr_history = emissivity_correction(T_cam_history, epsilon_set, epsilon_real)

    T_cam_imshow = T_cam_history.T
    T_corr_imshow = T_corr_history.T

    # Plots: temperature profiles at selected times
    fig_profiles = plot_temperature_profiles(
        x, T_history, t_array, plot_times,
        T_cam_history, T_corr_history, L1, L2
    )
    fig_profiles.savefig('temperature_profiles.png', dpi=150, bbox_inches='tight')
    plt.show()

    # Pseudo-IR colormaps
    fig_real = plot_IR_colormap(x, t_array, T_real_imshow,
                               'Pseudo-IR: $T_{real}(x,t)$', L1, L2)
    fig_real.savefig('colormap_T_real.png', dpi=150, bbox_inches='tight')

    fig_cam = plot_IR_colormap(x, t_array, T_cam_imshow,
                               'Pseudo-IR: $T_{cam}(x,t)$', L1, L2)
    fig_cam.savefig('colormap_T_cam.png', dpi=150, bbox_inches='tight')

    fig_corr = plot_IR_colormap(x, t_array, T_corr_imshow,
                                'Pseudo-IR: $T_{corr}(x,t)$', L1, L2)
    fig_corr.savefig('colormap_T_corr.png', dpi=150, bbox_inches='tight')

    plt.show()
    print("Done. Check temperature_profiles.png and colormap_*.png.")


if __name__ == "__main__":
    main()
