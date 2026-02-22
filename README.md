# 1D Transient Heat Conduction + IR Camera Simulation + Emissivity Correction  
# 一维瞬态导热与红外相机测量模拟及发射率修正

A teaching-oriented numerical simulation that couples **1D transient heat conduction** in a two-layer composite with **infrared camera measurement** (emissivity error) and **emissivity correction**.  
面向教学的数值模拟：将**一维瞬态导热**（双层复合材料）、**红外相机测量**（发射率误差）与**发射率修正**串联在一起。

---

## Overview | 概述

**English:**  
This project simulates heat flow in a 1D bar made of two materials in contact (e.g. metal + insulator). The “true” temperature field is computed with the heat equation, then a simplified IR camera model is applied: the camera sees radiation intensity and infers temperature using an assumed (wrong) emissivity. The script then applies the theoretical emissivity correction and compares true, measured, and corrected temperatures in both profile plots and pseudo-IR colormaps.

**中文：**  
本项目模拟由两种材料接触组成的一维杆中的热传导（如金属 + 绝缘体）。先用导热方程计算“真实”温度场，再施加简化的红外相机模型：相机根据辐射强度、用设定（错误）的发射率反算温度。脚本随后应用理论发射率修正，并在曲线图和伪红外热图中对比真实温度、相机测量温度与修正后温度。

---

## Physics | 物理模型

### Heat conduction | 导热

- **Governing equation 控制方程**
  <img width="964" height="338" alt="image" src="https://github.com/user-attachments/assets/e99daa13-8c47-442d-99f4-769cc6d178a6" />



- **Domain 计算域**
  - **Material A 材料 A:** $x \in [0, L_1]$ — thermal conductivity $k_A$, density $\rho_A$, specific heat $c_{p,A}$. 导热系数 $k_A$，密度 $\rho_A$，比热 $c_{p,A}$。
  - **Material B 材料 B:** $x \in [L_1, L_1+L_2]$ — $k_B$, $\rho_B$, $c_{p,B}$. $k_B$，$\rho_B$，$c_{p,B}$。
  - The interface at $x = L_1$ is handled by a spatially varying $k$ (harmonic mean at cell faces). 界面在 $x = L_1$，由空间变化的 $k$（界面处调和平均）自然满足。

- **Boundary conditions 边界条件**
  - $T(0, t) = T_{\mathrm{hot}}$ (hot end 热端)
  - $T(L_1 + L_2, t) = T_{\mathrm{cold}}$ (cold end 冷端)

- **Initial condition 初始条件**
  - $T(x, 0) = T_0$ (uniform 均匀).

### IR camera model | 红外相机模型

- **True radiation intensity 真实辐射强度** (Stefan–Boltzmann, gray body)  
  $
  I(x,t) = \varepsilon_{\mathrm{real}}(x) \sigma T_{\mathrm{real}}^4(x,t).
  $
  $\varepsilon_{\mathrm{real}}$ depends on material (A or B); $\sigma$ is the Stefan–Boltzmann constant. $\varepsilon_{\mathrm{real}}$ 随材料（A 或 B）变化；$\sigma$ 为斯特藩–玻尔兹曼常数。

- **Camera-inferred temperature 相机反算温度**  
  The camera assumes a single emissivity $\varepsilon_{\mathrm{set}}$ and computes. 相机假定单一发射率 $\varepsilon_{\mathrm{set}}$，并计算：
  $
  T_{\mathrm{cam}} = \left( \frac{I}{\varepsilon_{\mathrm{set}} \sigma} \right)^{1/4}.
  $
  If $\varepsilon_{\mathrm{set}} \neq \varepsilon_{\mathrm{real}}$, then $T_{\mathrm{cam}} \neq T_{\mathrm{real}}$. 若 $\varepsilon_{\mathrm{set}} \neq \varepsilon_{\mathrm{real}}$，则 $T_{\mathrm{cam}} \neq T_{\mathrm{real}}$。

### Emissivity correction | 发射率修正

- **Theoretical correction 理论修正**  
  From $I = \varepsilon_{\mathrm{real}} \sigma T_{\mathrm{real}}^4$ and $T_{\mathrm{cam}}^4 = I/(\varepsilon_{\mathrm{set}} \sigma)$ we get the correction. 由 $I = \varepsilon_{\mathrm{real}} \sigma T_{\mathrm{real}}^4$ 与 $T_{\mathrm{cam}}^4 = I/(\varepsilon_{\mathrm{set}} \sigma)$ 可得修正公式：
  $
  T_{\mathrm{corr}} = \left( \frac{\varepsilon_{\mathrm{set}}}{\varepsilon_{\mathrm{real}}} \right)^{1/4} T_{\mathrm{cam}}.
  $
  So $T_{\mathrm{corr}} = T_{\mathrm{real}}$ when the model (gray body, no reflection) holds. 在灰体、无反射等假设下，有 $T_{\mathrm{corr}} = T_{\mathrm{real}}$。

---

## Numerical method | 数值方法

- **Spatial discretization 空间离散:** 1D grid with $N$ points on $[0, L_1+L_2]$, spacing $\Delta x$. 在 $[0, L_1+L_2]$ 上均匀布 $N$ 个点，间距 $\Delta x$。

- **Time discretization 时间离散:** Explicit forward Euler; time step $\Delta t = t_{\mathrm{end}}/n$ with $n$ = `num_time_steps`. 显式前向欧拉；步长 $\Delta t = t_{\mathrm{end}}/n$，其中 $n$ 为 `num_time_steps`。

- **Finite-difference scheme 差分格式:**  
  At interior points, flux at “half” faces uses harmonic-mean conductivity. 内点处，界面热流用调和平均导热系数：
  $$
  k_{i+1/2} = \frac{2 k_i k_{i+1}}{k_i + k_{i+1}}, \quad
  F_{i+1/2} = -k_{i+1/2} \frac{T_{i+1} - T_i}{\Delta x}.
  $$  
  Then $\rho c_p \frac{\partial T}{\partial t}\big|_i \approx (F_{i+1/2} - F_{i-1/2})/\Delta x$, and $T_i^{n+1} = T_i^n + \Delta t (\partial T/\partial t)|_i$. 即 $\rho c_p \partial T/\partial t|_i \approx (F_{i+1/2} - F_{i-1/2})/\Delta x$，$T_i^{n+1} = T_i^n + \Delta t (\partial T/\partial t)|_i$。

- **Stability 稳定性:** The script checks $\Delta t \le 0.4 \Delta x^2 / \max(\alpha)$. If this fails, increase `num_time_steps` or decrease `N`. 脚本会检查 $\Delta t \le 0.4 \Delta x^2 / \max(\alpha)$；若不满足，请增大 `num_time_steps` 或减小 `N`。

---

## File structure | 文件结构

| File | Description |
|------|-------------|
| `thermo_ir_simulation.py` | Main script: grid, heat solve, IR simulation, correction, plotting. 主程序：网格、导热求解、红外模拟、修正与绘图。 |
| `README.md` | This file. 本说明。 |

**Outputs (after running 运行后生成):**

| File | Description |
|------|-------------|
| `temperature_profiles.png` | $T_{\mathrm{real}}$, $T_{\mathrm{cam}}$, $T_{\mathrm{corr}}$ vs $x$ at several times. 若干时刻的温度曲线对比。 |
| `colormap_T_real.png` | Pseudo-IR colormap of $T_{\mathrm{real}}(x,t)$. 真实温度的伪红外热图。 |
| `colormap_T_cam.png` | Pseudo-IR colormap of $T_{\mathrm{cam}}(x,t)$. 相机温度的伪红外热图。 |
| `colormap_T_corr.png` | Pseudo-IR colormap of $T_{\mathrm{corr}}(x,t)$. 修正后温度的伪红外热图。 |

---

## Requirements | 依赖

- **Python** 3.x  
- **NumPy**  
- **Matplotlib**  

Install (if needed):  
安装（如需要）：  
```bash
pip install numpy matplotlib
```

---

## Usage | 使用方法

Run the script from the project directory:  
在项目目录下运行：

```bash
python thermo_ir_simulation.py
```

The script will compute the solution, show figures, and save the four PNG files above.  
脚本将完成计算、弹出图形窗口并保存上述四个 PNG 文件。

---

## Main parameters | 主要参数

All parameters are defined at the top of `thermo_ir_simulation.py`. Key ones:  
所有参数在 `thermo_ir_simulation.py` 顶部定义。主要参数如下：

| Parameter | Meaning | Default |
|-----------|---------|---------|
| `L1`, `L2` | Lengths of material A and B [m] 材料 A、B 长度 | 0.05, 0.05 |
| `N` | Number of spatial grid points 空间网格点数 | 51 |
| `k_A`, `rho_A`, `cp_A` | Material A: conductivity, density, specific heat 材料 A：导热系数、密度、比热 | 50, 7800, 500 |
| `k_B`, `rho_B`, `cp_B` | Material B 材料 B | 0.2, 1000, 2000 |
| `Thot`, `Tcold`, `T0` | Hot/cold BC and initial temperature [K] 热端/冷端边界及初始温度 | 400, 300, 350 |
| `t_end` | Simulation end time [s] 模拟结束时间 | 100 |
| `num_time_steps` | Number of time steps 时间步数 | 1500 |
| `plot_times` | Times for profile plots [s] 绘制曲线的时刻 | [5, 20, 50, 100] |
| `epsilon_real_A`, `epsilon_real_B` | Real emissivity of A and B 材料 A、B 的真实发射率 | 0.85, 0.95 |
| `epsilon_set` | Emissivity assumed by camera 相机设定发射率 | 0.90 |

---

## Code structure (functions) | 代码结构（函数）

| Function | Role |
|----------|------|
| `setup_grid(L1, L2, N)` | Build 1D grid $x$, $\Delta x$, $L_{\mathrm{total}}$. 建立一维网格。 |
| `assign_material_properties(...)` | Assign $k$, $\rho c_p$, $\alpha$, $\varepsilon_{\mathrm{real}}$ on the grid. 在网格上分配物性。 |
| `solve_heat_conduction(...)` | Explicit FDM for transient conduction; returns $T(x,t)$ history. 显式 FDM 求解瞬态导热，返回温度历程。 |
| `compute_radiation_intensity(T_real, epsilon_real, sigma_sb)` | $I = \varepsilon_{\mathrm{real}}\, \sigma\, T_{\mathrm{real}}^4$. 计算辐射强度。 |
| `compute_camera_temperature(I, epsilon_set, sigma_sb)` | $T_{\mathrm{cam}} = (I/(\varepsilon_{\mathrm{set}}\sigma))^{1/4}$. 相机反算温度。 |
| `emissivity_correction(T_cam, epsilon_set, epsilon_real)` | $T_{\mathrm{corr}} = (\varepsilon_{\mathrm{set}}/\varepsilon_{\mathrm{real}})^{1/4}\, T_{\mathrm{cam}}$. 发射率修正。 |
| `plot_temperature_profiles(...)` | Plot $T_{\mathrm{real}}$, $T_{\mathrm{cam}}$, $T_{\mathrm{corr}}$ vs $x$ at selected times. 绘制选定时刻的温度曲线。 |
| `plot_IR_colormap(...)` | Pseudo-IR image (e.g. `inferno` colormap). 伪红外热图。 |

---

## Teaching notes | 教学说明

- **Heat equation:** Reinforces $\rho c_p \partial T/\partial t = \nabla\cdot(k\nabla T)$ and interface treatment with variable $k$.  
  **导热方程：** 巩固 $\rho c_p \partial T/\partial t = \nabla\cdot(k\nabla T)$ 及变 $k$ 的界面处理。

- **IR measurement:** Shows how a single “wrong” $\varepsilon_{\mathrm{set}}$ leads to systematic error when $\varepsilon_{\mathrm{real}}$ varies in space.  
  **红外测量：** 说明当 $\varepsilon_{\mathrm{real}}$ 随空间变化时，单一的“错误”$\varepsilon_{\mathrm{set}}$ 如何造成系统误差。

- **Correction:** Demonstrates that knowing $\varepsilon_{\mathrm{real}}(x)$ allows exact recovery of $T_{\mathrm{real}}$ from $T_{\mathrm{cam}}$ under the gray-body model.  
  **修正：** 在灰体模型下，已知 $\varepsilon_{\mathrm{real}}(x)$ 即可从 $T_{\mathrm{cam}}$ 精确恢复 $T_{\mathrm{real}}$。

- **Stability:** The explicit FDM stability condition is a good opportunity to discuss CFL-type constraints in parabolic PDEs.  
  **稳定性：** 显式 FDM 的稳定性条件是讨论抛物型 PDE 中 CFL 型约束的好例子。

---

## License | 许可

Use freely for teaching and learning. 可自由用于教学与学习。
