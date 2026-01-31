from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from cycler import cycler
from scipy.constants import gas_constant
from scipy.interpolate import interp1d

from postprocessing.read import read_physical_quantities, read_mesh, read_physical_quantity, read_solution_matrices


colors = list(plt.cm.tab10.colors)
linestyles = ['-', '--', '-.']

# Repeat patterns to match the number of colors (10)
ls10 = (linestyles * (len(colors) // len(linestyles) + 1))[:len(colors)]

custom_cycler = (
    cycler(color=colors) 
)

# ==================================================================================================
# UTILITY FUNCTIONS
# ==================================================================================================

def get_wave_properties(x_grid, Mach_profile, Temp_profile, target_M_values):
    """
    Interpolates to find the position x and Temperature T 
    where the flow reaches specific target Mach numbers.
    """
    # Ensure the Mach profile is sorted for interpolation (required by interp1d)
    # We sort by Mach number to handle the function x = f(M)
    sort_idx = np.argsort(Mach_profile)
    M_sorted = Mach_profile[sort_idx]
    x_sorted = x_grid[sort_idx]
    T_sorted = Temp_profile[sort_idx]

    # Create interpolators
    # bounds_error=False returns NaN if the Mach number isn't found in this step
    f_x = interp1d(M_sorted, x_sorted, kind='linear', bounds_error=False, fill_value=np.nan)
    f_T = interp1d(M_sorted, T_sorted, kind='linear', bounds_error=False, fill_value=np.nan)

    return f_x(target_M_values), f_T(target_M_values)

# ==================================================================================================
# SINGLE PROFILES
# ==================================================================================================

def draw_density_profile(x, folder, xlims = None, ylims = None, horizontal_line=None, save_path=None, show=False):
    """
    Draw the density profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    
    density, _, _ = read_physical_quantities(folder + "physical_quantities.txt")
    
    fig, ax = plt.subplots()
    ax.plot(x, density, label='Density')
    
    if horizontal_line is not None:
        ax.axhline(y=horizontal_line, color='k', linestyle='--', linewidth=1.0, label=f'Reference Line y={horizontal_line:.6f}')
    
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$\rho \ / \ \rho_w$')
    ax.set_title('Density Profile')
    ax.legend()
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_velocity_profile(x, folder, T_infty_w, xlims = None, ylims = None, horizontal_line=None, save_path=None, show=False):
    """
    Draw the velocity profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    
    _, velocity, _ = read_physical_quantities(folder + "physical_quantities.txt")
    
    fig, ax = plt.subplots()
    ax.plot(x, - velocity * (1 / sqrt((5.0 / 6.0) * T_infty_w)), label='Velocity')
    
    if horizontal_line is not None:
        ax.axhline(y=horizontal_line, color='k', linestyle='--', linewidth=1.0, label=f'Reference Line y={horizontal_line:.6f}')
    
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ - v_1 \ / \ a_\infty$')
    ax.set_title('Velocity Profile')
    ax.legend()
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

def draw_velocity_profile2(x, folder, T_infty_w, xlims = None, ylims = None, horizontal_line=None, save_path=None, show=False):
    """
    Draw the velocity profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    
    _, velocity, _ = read_physical_quantities(folder + "physical_quantities.txt")
    
    fig, ax = plt.subplots()
    ax.plot(x, velocity, label='Velocity')
    
    if horizontal_line is not None:
        ax.axhline(y=horizontal_line, color='k', linestyle='--', linewidth=1.0, label=f'Reference Line y={horizontal_line:.6f}')
    
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ v_1 \ / \ \sqrt{{2RT_w}}$')
    ax.set_title('Velocity Profile')
    ax.legend()
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_temperature_profile(x, folder, xlims = None, ylims = None, horizontal_line=None, save_path=None, show=False):
    """
    Draw the temperature profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    
    _, _, temperature = read_physical_quantities(folder + "physical_quantities.txt")
    
    fig, ax = plt.subplots()
    ax.plot(x, temperature, label='Temperature')
    
    if horizontal_line is not None:
        ax.axhline(y=horizontal_line, color='k', linestyle='--', linewidth=1.0, label=f'Reference Line y={horizontal_line:.6f}')
    
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ T \ / \ T_w$')
    ax.set_title('Temperature Profile')
    ax.legend()
    
    if xlims:   
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_pressure_profile(x, folder, xlims = None, ylims = None, horizontal_line=None, save_path=None, show=False):
    """
    Draw the pressure profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    
    density, _, temperature = read_physical_quantities(folder + "physical_quantities.txt")
    
    fig, ax = plt.subplots()
    ax.plot(x, density * temperature, label='Pressure')
    
    if horizontal_line is not None:
        ax.axhline(y=horizontal_line, color='k', linestyle='--', linewidth=1.0, label=f'Reference Line y={horizontal_line:.6f}')
    
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ p \ / \ p_w$')
    ax.set_title('Pressure Profile')
    ax.legend()
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# ==================================================================================================
# EVOLUTION PLOTS
# ==================================================================================================

def draw_density_evolution(x, timesteps, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the density evolution over time using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    fig, ax = plt.subplots()
    for i in timesteps:
        density, _, _ = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, density, label=rf'$\overline{{t}} = {i * dt:.2f}$')
        ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
        ax.set_ylabel(rf'$\rho \ / \ \rho_w$')
    ax.set_title('Density Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_temperature_evolution(x, timesteps, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the temperature evolution over time using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    fig, ax = plt.subplots()
    for i in timesteps:
        _, _, temp = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, temp, label=rf'$\overline{{t}} = {i * dt:.2f}$')
        ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
        ax.set_ylabel(rf'$ T \ / \ T_w$')
    ax.set_title('Temperature Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_velocity_evolution(x, timesteps, T_infty_w, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the velocity evolution over time using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    fig, ax = plt.subplots()
    for i in timesteps:
        _, vel, _ = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, - vel * (1 / sqrt((5.0 / 6.0) * T_infty_w)),
                label=rf'$\overline{{t}} = {i * dt:.2f}$')
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ - v_1 \ / \ a_\infty$')
    ax.set_title('Velocity Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)
        
def draw_velocity_evolution2(x, timesteps, T_infty_w, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the velocity evolution over time using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    fig, ax = plt.subplots()
    for i in timesteps:
        _, vel, _ = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, vel, label=rf'$\overline{{t}} = {i * dt:.2f}$')
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ v_1 \ / \ \sqrt{{2RT_w}}$')
    ax.set_title('Velocity Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_pressure_evolution(x, timesteps, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the pressure evolution over time using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    fig, ax = plt.subplots()
    for i in timesteps:
        density, _, temp = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, density * temp, label=rf'$\overline{{t}} = {i * dt:.2f}$')
    ax.set_xlabel(rf'$ X_1 \ / \ l_w$')
    ax.set_ylabel(rf'$ p \ / \ p_w$')
    ax.set_title('Pressure Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_solution_discontinuity_evolution(x, zeta, x_indices, timesteps, dt, zeta_lim, folder, save_path=None, show=False):
    """ 
    Produces N plots (one per timestep). Each plot shows the distribution 
    over zeta for all specified spatial indices (x_indices).
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        zeta (array-like) : The transformed velocity variable grid points.
        x_indices (list of int) : The spatial indices to plot (e.g., [10, 20, 30]).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "data").
        save_path (str, optional): Full path including filename to save the figure (e.g., "plots/solution_discontinuity.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    
    plt.rcParams['axes.prop_cycle'] = custom_cycler

    for i in timesteps:
        # Create a new figure for each timestep
        fig, (ax_g, ax_h) = plt.subplots(
            nrows=1, ncols=2, figsize=(12, 5), sharex=True
        )

        # 1. Read the data ONCE per timestep for efficiency
        g_data = read_solution_matrices(f"{folder}g_iter_{i}.txt")
        h_data = read_solution_matrices(f"{folder}h_iter_{i}.txt")

        # 2. Loop through each spatial index to plot multiple lines on these axes
        for idx in x_indices:
            label_text = rf'$x = {x[idx]:.2f}$'
            ax_g.plot(zeta, g_data[:, idx], label=label_text)
            ax_h.plot(zeta, h_data[:, idx], label=label_text)

        # ---- Formatting: g plot ----
        ax_g.set_xlabel(rf'$\zeta$')
        ax_g.set_ylabel(rf'$g$')
        ax_g.set_title(rf'Distributions of $g$ at $\overline{{t}} = {i * dt:.2f}$')
        ax_g.legend(loc="best", fontsize='small')
        if zeta_lim:
            ax_g.set_xlim(zeta_lim)
            ax_h.set_xlim(zeta_lim)

        # ---- Formatting: h plot ----
        ax_h.set_xlabel(rf'$\zeta$')
        ax_h.set_ylabel(rf'$h$')
        ax_h.set_title(rf'Distributions of $h$ at $\overline{{t}} = {i * dt:.2f}$')
        ax_h.legend(loc="best", fontsize='small')

        fig.tight_layout()

        # ---- Saving Logic ----
        if save_path:
            p = Path(save_path)
            # Append the iteration number to the filename so files don't overwrite
            unique_save_path = p.parent / f"{p.stem}_iter_{i}{p.suffix}"
            p.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(unique_save_path, dpi=300, bbox_inches='tight')

        if show:
            plt.show()
        else:
            plt.close(fig) # Critical to prevent memory bloat
    

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_MachNumber_evolution(x, timesteps, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draws the spatial evolution of the Mach number over specified time steps.

    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "data").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        save_path (str, optional): Full path including filename to save the figure (e.g., "plots/mach_evolution.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    fig, ax = plt.subplots()
    for i in timesteps:
        density, velocity, temp = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        c = np.sqrt((6.0 / 5.0) / temp)
        Mach = np.abs(velocity) * c
        ax.plot(x, Mach, label=rf'$\overline{{t}} = {i * dt:.2f}$')
    ax.set_xlabel(r'$ X_1 \ / \ l_w$')
    ax.set_ylabel('Mach Number')
    ax.set_title('Mach Number Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

def draw_Vp_star_evolution(x, timesteps, m, dt, folder, target_Machs=np.linspace(0.98, 1.02, 100), xlims=None, ylims=None, save_path=None, show=False):
    """
    Calculates and plots $V_p^*$ against Mach number.
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The specific time steps (iterations) you want to visualize (plot).
        m (int) : The step lag used for the finite difference derivative. 
                  Approximation: dx/dt ~ (x(t) - x(t-m)) / (m * dt)
        dt (float) : Physical time per iteration step.
        folder (str) : Path to the folder containing the output files (e.g., "data").
        target_Machs (array-like) : Mach numbers at which to evaluate Vp*.
        xlims (tuple, optional) : (min, max) limits for the x-axis.
        ylims (tuple, optional) : (min, max) limits for the y-axis.
        save_path (str, optional) : Full path including filename to save the figure (e.g., "plots/Vp_star_evolution.png").
        show (bool, default False) : If True, calls plt.show(). If False, closes the figure to free memory.
    """
    folder_path = Path(folder)
    plt.rcParams['axes.prop_cycle'] = custom_cycler
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for t_curr in timesteps:
        t_prev = t_curr - m
        
        # Sanity check: We can't look back before step 0
        if t_prev < 0:
            print(f"Skipping t={t_curr}: Lagged step {t_prev} is negative.")
            continue
            
        # 1. Define Paths
        path_curr = folder_path / f"phys_iter_{t_curr}.txt"
        path_prev = folder_path / f"phys_iter_{t_prev}.txt"
        
        # 2. Read Data
        _, v_curr, T_curr = read_physical_quantities(path_curr)
        _, v_prev, T_prev = read_physical_quantities(path_prev)
        
        if v_curr is None or v_prev is None:
            print(f"Skipping t={t_curr}: Missing file (checked {t_curr} and {t_prev}).")
            continue

        # 3. Calculate Mach Profiles
        c_curr = np.sqrt((6.0 / 5.0) / T_curr)
        c_prev = np.sqrt((6.0 / 5.0) / T_prev)
        
        M_curr = np.abs(v_curr) * c_curr
        M_prev = np.abs(v_prev) * c_prev
        
        # 4. Interpolate x(M) and T(M)
        x_locs_curr, T_locs_curr = get_wave_properties(x, M_curr, T_curr, target_Machs)
        x_locs_prev, _           = get_wave_properties(x, M_prev, T_prev, target_Machs)
        
        # 5. Compute Derivative and V_p*
        # Physical time interval = m * dt
        delta_t_phys = m * dt
        
        # Wave speed = (x_curr - x_prev) / delta_t
        wave_speed = (x_locs_curr - x_locs_prev) / delta_t_phys
        
        # Normalize: V_p* = sqrt(6/5) * (1/sqrt(T)) * wave_speed
        normalization = np.sqrt(6.0 / 5.0) * (1.0 / np.sqrt(T_locs_curr))
        V_p_star = normalization * wave_speed
        
        # 6. Plot
        valid = ~np.isnan(V_p_star)
        label_str = rf'$\overline{{t}} = {t_curr * dt:.0f}$'
        ax.plot(target_Machs[valid], V_p_star[valid], label=label_str, linewidth=1.5)

    # --- Reference Line: 1 - M ---
    subsonic = target_Machs <= 1.0
    ax.plot(target_Machs[subsonic], 1.0 - target_Machs[subsonic], 
            color='black', linestyle='--', linewidth=2, label=r'$1-M$')

    # --- Formatting ---
    ax.set_xlabel('Mach Number ($M$)', fontsize=12)
    ax.set_ylabel(r'$V_p^*(\overline{t}, M)$', fontsize=12)
    ax.set_title(f'Normalized Wave Velocity (lag $m={m}$)', fontsize=14)
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=10)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



if __name__ == "__main__":
    # Example usage
    folder = "output/first_test/"
    x = read_mesh(folder + "space_mesh.txt")

#    density, velocity, temperature = read_physical_quantities(
#        folder + "physical_quantities.txt")
#
#    draw_temperature_profile(
#        x, temperature, save_path=folder + "temperature_profile.png")
#
#    draw_density_profile(
#        x, density, save_path=folder + "density_profile.png")
#
#    draw_velocity_profile(
#        x, velocity, save_path=folder + "velocity_profile.png")

    timesteps = [1, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
    dt = 0.35
    draw_temperature_evolution(
        x, timesteps, dt, folder, save_path=folder + "temperature_evolution.png")
    draw_velocity_evolution(
        x, timesteps, dt, folder, save_path=folder + "velocity_evolution.png")
    print("Plots saved.")
