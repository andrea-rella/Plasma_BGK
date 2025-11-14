from math import sqrt
import time
import matplotlib.pyplot as plt
from pathlib import Path
from cycler import cycler
from scipy.constants import gas_constant

from postprocessing.read import read_physical_quantities, read_space_mesh, read_physical_quantity

# ==================================================================================================
# SINGLE PROFILES
# ==================================================================================================

def draw_density_profile(x, density, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the density profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    fig, ax = plt.subplots()
    ax.plot(x, density, label='Density')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('Density')
    ax.set_title('Density Profile')
    ax.legend()
    ax.grid()
    
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

def draw_velocity_profile(x, velocity, T_infty_w, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the velocity profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    fig, ax = plt.subplots()
    ax.plot(x, - velocity * (1 / sqrt((5.0 / 6.0) * T_infty_w)), label='Velocity')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('- v1 / a_infty')
    ax.set_title('Velocity Profile')
    ax.legend()
    ax.grid()
    
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

def draw_temperature_profile(x, temperature, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the temperature profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    fig, ax = plt.subplots()
    ax.plot(x, temperature, label='Temperature')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('T / T_w')
    ax.set_title('Temperature Profile')
    ax.legend()
    ax.grid()
    
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

def draw_pressure_profile(x, density, temperature, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the pressure profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    fig, ax = plt.subplots()
    ax.plot(x, density * temperature, label='Pressure')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('p / p_w')
    ax.set_title('Pressure Profile')
    ax.legend()
    ax.grid()
    
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

def draw_temperature_evolution(x, timesteps, dt, folder, xlims = None, ylims = None, save_path=None, show=False):
    """
    Draw the temperature evolution over time using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    plt.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    for i in timesteps:
        _, _, temp = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, temp, label=f'Time {i * dt:.2f}')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('T / T_w')
    ax.set_title('Temperature Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.grid()
    
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
    plt.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    for i in timesteps:
        _, vel, _ = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, - vel * (1 / sqrt((5.0 / 6.0) * T_infty_w)),
                label=f'Time {i * dt:.2f}')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('- v1 / a_infty')
    ax.set_title('Velocity Profile')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.grid()
    
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
    plt.rcParams['axes.prop_cycle'] = cycler(color=plt.cm.tab20.colors)
    fig, ax = plt.subplots()
    for i in timesteps:
        density, _, temp = read_physical_quantities(
            folder + "/phys_iter_" + str(i) + ".txt")
        ax.plot(x, density * temp, label=f'Time {i * dt:.2f}')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('p / p_w')
    ax.set_title('Pressure Evolution')
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    ax.grid()
    
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




if __name__ == "__main__":
    # Example usage
    folder = "output/first_test/"
    x = read_space_mesh(folder + "space_mesh.txt")

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
