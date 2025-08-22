import matplotlib.pyplot as plt
from pathlib import Path

from read import read_physical_quantities, read_space_mesh


def draw_temperature_profile(x, temperature, save_path=None, show=False):
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

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)


def draw_density_profile(x, density, save_path=None, show=False):
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

    if save_path:
        p = Path(save_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)


def draw_velocity_profile(x, velocity, save_path=None, show=False):
    """
    Draw the velocity profile using Matplotlib.
    If save_path is provided, save the figure as a PNG there.
    """
    fig, ax = plt.subplots()
    ax.plot(x, velocity, label='Velocity')
    ax.set_xlabel('X_1 / l_w')
    ax.set_ylabel('Velocity')
    ax.set_title('Velocity Profile')
    ax.legend()
    ax.grid()

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

    density, velocity, temperature = read_physical_quantities(
        folder + "physical_quantities.txt")

    draw_temperature_profile(
        x, temperature, save_path=folder + "temperature_profile.png")

    draw_density_profile(
        x, density, save_path=folder + "density_profile.png")

    draw_velocity_profile(
        x, velocity, save_path=folder + "velocity_profile.png")
