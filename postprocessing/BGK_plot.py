from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from cycler import cycler
from scipy.constants import gas_constant
from scipy.interpolate import interp1d

from postprocessing.BGK_read import read_physical_quantities, read_mesh, read_physical_quantity, read_solution_matrices


colors = list(plt.cm.tab10.colors)
linestyles = ['-', '--', '-.']

# Repeat patterns to match the number of colors (10)
ls10 = (linestyles * (len(colors) // len(linestyles) + 1))[:len(colors)]

custom_cycler = (
    cycler(color=colors) 
)


# ==================================================================================================
# SINGLE PROFILES
# ==================================================================================================

def draw_density_profile(x, folder, xlims = None, ylims = None, horizontal_line=None, save_path=None, show=False):
    """
    Draw the density profile using Matplotlib. If save_path is provided, save the figure as a PNG there.
    (supposes in the folder there is a file "physical_quantities.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like): The spatial grid points (X1/lw).
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        horizontal_line (float, optional): y-value at which to draw a horizontal reference line.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/density_profile.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    
    density, _, _ = read_physical_quantities(folder + "/physical_quantities.txt")
    
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
    Draw the velocity profile normalized by the sound speed at infinity using Matplotlib. If save_path is provided, 
    save the figure as a PNG there. (supposes in the folder there is a file "physical_quantities.txt" with density, 
    velocity, temperature columns)
    
    Args:
        x (array-like): The spatial grid points (X1/lw).
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        T_infty_w (float): The reference temperature for normalization.
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        horizontal_line (float, optional): y-value at which to draw a horizontal reference line.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/velocity_profile.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    
    _, velocity, _ = read_physical_quantities(folder + "/physical_quantities.txt")
    
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
    Draw the velocity profile normalized by the local sound speed at wall temperature using Matplotlib. If save_path is provided, 
    save the figure as a PNG there. (supposes in the folder there is a file "physical_quantities.txt" with density, 
    velocity, temperature columns)
    
    Args:
        x (array-like): The spatial grid points (X1/lw).
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        T_infty_w (float): The reference temperature for normalization.
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        horizontal_line (float, optional): y-value at which to draw a horizontal reference line.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/velocity_profile.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    
    _, velocity, _ = read_physical_quantities(folder + "/physical_quantities.txt")
    
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
    Draw the normalized temperature profile using Matplotlib. If save_path is provided, save the figure as a PNG there.
    (supposes in the folder there is a file "physical_quantities.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like): The spatial grid points (X1/lw).
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        horizontal_line (float, optional): y-value at which to draw a horizontal reference line.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/temperature_profile.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    
    _, _, temperature = read_physical_quantities(folder + "/physical_quantities.txt")
    
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
    Draw the normalized pressure profile using Matplotlib. If save_path is provided, save the figure as a PNG there.
    (supposes in the folder there is a file "physical_quantities.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like): The spatial grid points (X1/lw).
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        horizontal_line (float, optional): y-value at which to draw a horizontal reference line.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/pressure_profile.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
    """
    
    density, _, temperature = read_physical_quantities(folder + "/physical_quantities.txt")
    
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
    Draw the density evolution over time using Matplotlib. If save_path is provided, save the figure as a PNG there.
    (supposes in the folder there are files "phys_iter_#.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/density_evolution.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
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
    Draw the temperature evolution over time using Matplotlib. If save_path is provided, save the figure as a PNG there.
    (supposes in the folder there are files "phys_iter_#.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/temperature_evolution.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
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
    Draw the velocity normalized by the sound speed at infinity over time using Matplotlib. If save_path is provided, 
    save the figure as a PNG there. (supposes in the folder there are files "phys_iter_#.txt" with density, velocity, 
    temperature columns)
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        T_infty_w (float): The reference temperature for normalization.
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/velocity_evolution.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
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
    Draw the velocity evolution normalized by the local sound speed at the wall over time using Matplotlib. 
    If save_path is provided, save the figure as a PNG there. (supposes in the folder there are files 
    "phys_iter_#.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        T_infty_w (float): The reference temperature for normalization.
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/velocity_evolution.png").
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
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
    Draw the normalized pressure evolution over time using Matplotlib. If save_path is provided, save the figure 
    as a PNG there. (supposes in the folder there are files "phys_iter_#.txt" with density, velocity, temperature columns)
    
    Args:
        x (array-like) : The spatial grid points (X1/lw).
        timesteps (list of int) : The iteration numbers corresponding to the files to be read (e.g. [0, 500, 1000]).
        dt (float): The time step size (delta t) used to compute the physical time label.
        folder (str): Path to the folder containing the output files (e.g., "output/test1").
        xlims (tuple, optional): (min, max) limits for the x-axis.
        ylims (tuple, optional): (min, max) limits for the y-axis.
        save_path (str, optional): Full path including filename to save the figure (e.g., "output/test1/pressure_evolution.png").  
        show (bool, default False): If True, calls plt.show(). If False, closes the figure to free memory.
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
    Produces N plots (one per timestep). Each plot shows the distribution over zeta for all specified spatial 
    indices (x_indices). (Supposes in the folder there are files "g_iter_#.txt" and "h_iter_#.txt" with solution matrices)
    
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
    (supposes in the folder there are files "phys_iter_#.txt" with density, velocity, temperature columns)

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




# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------



if __name__ == "__main__":

    import postprocessing.BGK_read as rd
    import json
    from math import sqrt
    import numpy as np

    with open("data/Cond/type1.json", "r", encoding="utf-8") as f:
        problem_data = json.load(f)

    folder = problem_data["general"]["saving_folder_name"] + "/"

    # ============== EVOLUTION ==============

    # Type 1 
    timesteps = [0, 20, 80, 200, 400,700]
    dt = problem_data["simulation"]["time_step"]
    xlim = (0, 1.2)
    temp_y_lim = (0.98, 1.03)
    velocity_y_lim = (1.16, 1.25)
    pressure_y_lim = (16.2, 17.8)
    density_y_lim = (16.2, 17.7)
    
    x = rd.read_mesh(folder + "space_mesh.txt")

    draw_temperature_evolution(
        x, 
        timesteps, 
        dt, 
        folder, 
        xlims=xlim, 
        ylims=temp_y_lim, 
        save_path=folder + "temperature_evolution.png"
        )

    draw_velocity_evolution(
        x, 
        timesteps, 
        problem_data["physical"]["T_infty_w"], 
        dt, 
        folder, 
        xlims=xlim, 
        ylims=velocity_y_lim, 
        save_path=folder + "velocity_evolution.png"
        )

    draw_pressure_evolution(
        x, 
        timesteps, 
        dt, 
        folder, 
        xlims=xlim, 
        ylims=pressure_y_lim, 
        save_path=folder + "pressure_evolution.png"
        )

    draw_density_evolution(
        x, 
        timesteps, 
        dt, 
        folder, 
        xlims=xlim, 
        ylims=density_y_lim,
        save_path=folder + "density_evolution.png"
        )