import numpy as np
import os
import math
import matplotlib.pyplot as plt

# ==================================================================================================
# Default Aoki Wall Functions (x = 0)
# ==================================================================================================

def aoki_g_wall(zeta):
    """
    Computes g at the wall (x=0) for zeta > 0.
    Formula: g = (1/sqrt(pi)) * exp(-zeta^2)
    """
    return (1.0 / np.sqrt(np.pi)) * np.exp(-zeta**2)

# ==================================================================================================

def aoki_h_wall(zeta):
    """
    Computes h at the wall (x=0) for zeta > 0.
    Formula: h = (1/sqrt(pi)) * exp(-zeta^2)
    """
    return (1.0 / np.sqrt(np.pi)) * np.exp(-zeta**2)

# ==================================================================================================
# Default Aoki Downstream Functions (x = infinity)
# ==================================================================================================

def aoki_g_downstream(zeta, M_inf, p_ratio, T_ratio):
    """
    Computes g downstream (x -> inf) for zeta < 0.
    """
    # Shift velocity u_inf based on Mach number
    shift = np.sqrt((5.0 / 6.0) * T_ratio) * M_inf
    
    # Exponent: - (zeta - u_inf)^2 / (T_inf / T_w)
    exponent = - (zeta - shift)**2 / T_ratio
    
    # Pre-factor: (1/sqrt(pi)) * (p_inf/p_w) * (T_inf/T_w)^(-3/2)
    pre_factor = (1.0 / np.sqrt(np.pi)) * p_ratio * (T_ratio**(-1.5))
    
    return pre_factor * np.exp(exponent)

# ==================================================================================================

def aoki_h_downstream(zeta, M_inf, p_ratio, T_ratio):
    """
    Computes h downstream (x -> inf) for zeta < 0.
    """
    # Shift velocity u_inf based on Mach number
    shift = np.sqrt((5.0 / 6.0) * T_ratio) * M_inf
    
    # Exponent: - (zeta - u_inf)^2 / (T_inf / T_w)
    exponent = - (zeta - shift)**2 / T_ratio
    
    # Pre-factor: (1/sqrt(pi)) * (p_inf/p_w) * (T_inf/T_w)^(-5/2)
    # Note the power is -2.5 here vs -1.5 for g
    pre_factor = (1.0 / np.sqrt(np.pi)) * p_ratio * (T_ratio**(-2.5))
    
    return pre_factor * np.exp(exponent)

# ==================================================================================================
# Compute theoretical equilibrium parameters
# ==================================================================================================

def compute_ytrehus_parameters(M_inf):
    """ Computes the parameters needed for the ytrehus formulation

    Args:
        M_inf (float): v_infty / sqrt(5/3 R T_w)
    """
    
    S_inf = M_inf * np.sqrt(5.0 / 6.0)
    sqrt_pi = math.sqrt(math.pi) 
    
    F_minus = sqrt_pi * S_inf * (math.erf(S_inf)-1) + math.exp(-S_inf**2)
    G_minus = (2 * S_inf**2 + 1) * (1 - math.erf(S_inf)) - (2 / sqrt_pi) * S_inf * math.exp(-S_inf**2)
    
    return S_inf, F_minus, G_minus
    

def compute_equilibrium_phys(M_inf):
    """
    Compute theoretical equilibrium parameters for the Ytrehus model.

    Parameters
    ----------
    M_inf : float
        Downstream Mach number.

    Returns
    -------
    S_inf : float
        Non-dimensional stream parameter, ``S_inf = M_inf * sqrt(5/6)``.
    p_ratio : float
        Pressure ratio ``p_inf / p_w``.
    T_ratio : float
        Temperature ratio ``T_inf / T_w``.
    rho_ratio : float
        Density ratio ``rho_inf / rho_w``.
    beta : float
        Backscatter factor.
    phi1 : float
        Auxiliary function ``phi1``.
    phi2 : float
        Auxiliary function ``phi2``.
    P : float
        Collision/decay parameter.
    r : float
        Auxiliary function ``r``.
    """
    
    sqrt_pi = math.sqrt(math.pi)
    S_inf, F_minus, G_minus = compute_ytrehus_parameters(M_inf)
    
    # T_inf / T_w
    T_ratio = ((- sqrt_pi / 8) * S_inf + math.sqrt(1 + (math.pi / 64.0)* S_inf**2))**2
    
    # p_inf / p_w
    p_ratio_inv = 2 * math.exp(-S_inf**2) / (F_minus + math.sqrt(T_ratio) * G_minus)
    p_ratio = 1.0 / p_ratio_inv
    
    # rho_inf / rho_w
    rho_ratio = p_ratio / T_ratio
    
    # beta
    beta = (2 * (2 * S_inf**2 + 1) * math.sqrt(T_ratio) - 2 * sqrt_pi * S_inf) / \
            (F_minus + math.sqrt(T_ratio) * G_minus)
            
    # phi1
    phi1 = (1 / rho_ratio - 2 + beta * (1-math.erf(S_inf))) / (beta - 1)
    
    # phi2
    phi2 = (1 / p_ratio - 2 + beta * (1-math.erf(S_inf))) / (beta - 1)
    
    # P
    P = ((math.pi / 12) * rho_ratio * (beta - 1)**2 * phi1 * phi2) / \
        ((1 / p_ratio) * (1 - T_ratio))
    
    # r
    r = 1 - (2 / phi1) + (4 * S_inf**2 / phi2)
    
    return S_inf, p_ratio, T_ratio, rho_ratio, beta, phi1, phi2, P, r
    
    
def plot_ytrehus_3d(x_limit, show = False, save_path = None):
    
    from scipy.special import erf
    
    # 1. Coordinate Setup
    # S_inf is the independent variable
    s_inf = np.linspace(0, x_limit, 500)
    sqrt_pi = np.sqrt(np.pi)
    
    # 2. Temperature Ratio Calculation (Eq 1)
    # Solve for sqrt(T_inf/T_w) first
    term1 = -(sqrt_pi / 8) * s_inf
    term2 = np.sqrt(1 + (np.pi / 64) * (s_inf**2))
    sqrt_t_ratio = term1 + term2
    t_ratio = sqrt_t_ratio**2  # T_inf / T_w
    
    # 3. Pressure Ratio Calculation (Eq 2)
    # Using your specific definitions for F_minus and G_minus
    f_minus = sqrt_pi * s_inf * (erf(s_inf) - 1) + np.exp(-(s_inf**2))
    g_minus = (2 * s_inf**2 + 1) * (1 - erf(s_inf)) - (2 / sqrt_pi) * s_inf * np.exp(-(s_inf**2))
    
    # Calculate p_w / p_inf based on the formula provided
    pw_pinf = (2 * np.exp(-(s_inf**2))) / (f_minus + g_minus * sqrt_t_ratio)
    
    # Invert to get p_inf / p_w as requested
    p_ratio = 1 / pw_pinf
    
    # 4. 3D Plotting
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the curve
    ax.plot(s_inf, t_ratio, p_ratio, lw=2.5, color='blue', label='Steady State Ytrehus Relations')
    
    # Labeling
    ax.set_xlabel(rf'$S_\infty$')
    ax.set_ylabel(rf'$T_\infty / T_w$')
    ax.set_zlabel(rf'$p_\infty / p_w$')
    ax.set_title(rf'Ytrehus Relations for $S_\infty \in [0, {x_limit}]$')
    
    # Add a grid and improve viewing angle
    ax.view_init(elev=20, azim=45)
    plt.legend()
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path + "ytrehus_3d_plot.png")
    if show:
        plt.show()
    else:
        plt.close()

# ==================================================================================================
# Main Solver Function
# ==================================================================================================  


def compute_gh(M_inf, 
               x_mesh, 
               zeta_mesh, 
               output_path=None, 
               g_wall_func=None, h_wall_func=None,
               g_downstream_func=None, h_downstream_func=None
               ):
    """
    Computes g and h using the Ytrehus model with separate boundary functions.
    
    Args:
        M_inf (float): Flow parameter.
        x_mesh (numpy.ndarray): Space grid.
        zeta_mesh (numpy.ndarray): Velocity grid.
        output_path (str, optional): Path to save .txt files.
        g_wall_func (callable): Function for g at the wall (x=0). 
        h_wall_func (callable): Function for h at the wall (x=0).
        g_downstream_func (callable): Function for g downstream (x -> inf).
        h_downstream_func (callable): Function for h downstream (x -> inf).
        
    Returns:
        g (numpy.ndarray): Computed g values (shape: len(zeta_mesh) x len(x_mesh)).
        h (numpy.ndarray): Computed h values (shape: len(zeta_mesh) x len(x_mesh)).
    """
    
    # 1. Set Defaults
    if g_wall_func is None: g_wall_func = aoki_g_wall
    if h_wall_func is None: h_wall_func = aoki_h_wall
    if g_downstream_func is None: g_downstream_func = aoki_g_downstream
    if h_downstream_func is None: h_downstream_func = aoki_h_downstream

    g = np.zeros((len(zeta_mesh), len(x_mesh)))
    h = np.zeros((len(zeta_mesh), len(x_mesh)))

    S_inf, p_ratio, T_ratio, rho_ratio, beta, phi1, phi2, P, r = compute_equilibrium_phys(M_inf)

    if M_inf >= 0.999:
        a_inf_minus = 1 + (beta - 1) / (1 + P * (beta - 1) * x_mesh)
    else:
        E_x = np.exp(-P * (1 - r) * x_mesh) 
        a_inf_minus = ((beta - r) - r * (beta - 1) * E_x) / \
                      ((beta - r) - (beta - 1) * E_x)

    a_e_plus = (a_inf_minus - 1) / (beta - 1)
    a_inf_plus = (beta - a_inf_minus) / (beta - 1)

    g_wall_vals = g_wall_func(zeta_mesh)
    h_wall_vals = h_wall_func(zeta_mesh)
    g_inf_vals = g_downstream_func(zeta_mesh, M_inf, p_ratio, T_ratio)
    h_inf_vals = h_downstream_func(zeta_mesh, M_inf, p_ratio, T_ratio)
    
    g_e_plus_dist = g_wall_vals[:, np.newaxis]
    h_e_plus_dist = h_wall_vals[:, np.newaxis]
    g_inf_dist = g_inf_vals[:, np.newaxis]
    h_inf_dist = h_inf_vals[:, np.newaxis]

    A_e_plus = a_e_plus[np.newaxis, :]
    A_inf_plus = a_inf_plus[np.newaxis, :]
    A_inf_minus = a_inf_minus[np.newaxis, :]

    g_pos_full = A_e_plus * g_e_plus_dist + A_inf_plus * g_inf_dist
    h_pos_full = A_e_plus * h_e_plus_dist + A_inf_plus * h_inf_dist
    
    g_neg_full = A_inf_minus * g_inf_dist
    h_neg_full = A_inf_minus * h_inf_dist
    
    pos_mask = (zeta_mesh > 0)
    neg_mask = (zeta_mesh < 0)
    
    g[pos_mask, :] = g_pos_full[pos_mask, :]
    g[neg_mask, :] = g_neg_full[neg_mask, :]
    
    h[pos_mask, :] = h_pos_full[pos_mask, :]
    h[neg_mask, :] = h_neg_full[neg_mask, :]
    
    if output_path is not None:
        os.makedirs(output_path, exist_ok=True)
        file_g = os.path.join(output_path, "eval_g_ytrehus.txt")
        file_h = os.path.join(output_path, "eval_h_ytrehus.txt")
        np.savetxt(file_g, g, header="Rows: Velocity Mesh, Cols: Space Mesh")
        np.savetxt(file_h, h, header="Rows: Velocity Mesh, Cols: Space Mesh")
        print(f"Files saved to: {output_path}")
    return g, h


def compute_physical_quantities(M_inf, 
                                x_mesh, 
                                output_path=None):
    """
    Computes macroscopic physical quantities (rho/rho_w, T/T_w) using the Ytrehus model.
    
    Args:
        M_inf (float): Downstream Mach number.
        x_mesh (numpy.ndarray): Space grid (X1 / l_w).
        output_path (str, optional): Path to save the .txt file.
        
    Returns:
        rho_rel_w (numpy.ndarray): Density ratio rho / rho_w.
        T_rel_w (numpy.ndarray): Temperature ratio T / T_w.
    """

    # 1. Compute parameters
    S_inf, p_ratio, T_ratio, rho_ratio, beta, phi1, phi2, P, r = compute_equilibrium_phys(M_inf)
    
    # 2. Amplitude function a_inf_minus 
    if M_inf >= 0.999:
        # Sonic case: algebraic relaxation 
        a_inf_minus = 1 + (beta - 1) / (1 + P * (beta - 1) * x_mesh)
    else:
        # Subsonic case: exponential relaxation 
        E_x = np.exp(-P * (1 - r) * x_mesh) 
        a_inf_minus = ((beta - r) - r * (beta - 1) * E_x) / ((beta - r) - (beta - 1) * E_x)

    # 3. Macroscopic Quantities 
    # Density ratio rho / rho_inf
    rho_rel_inf = (phi1 * (a_inf_minus - 1) / 2) + 1
    
    # Temperature ratio T / T_inf 
    T_rel_inf = (1/3) * (1/rho_rel_inf) * (2*S_inf**2 + 3 + phi2*(a_inf_minus - 1) - 2*S_inf**2*(1/rho_rel_inf))
    
    # Final ratios relative to wall values 
    rho_rel_w = rho_rel_inf * rho_ratio
    T_rel_w = T_rel_inf * T_ratio
    
    # Placeholder for mean velocity
    mean_velocity = np.ones_like(x_mesh)

    # 4. Saving to file
    if output_path is not None:
        os.makedirs(output_path, exist_ok=True)
        file_path = os.path.join(output_path, "ytrehus_physical_quantities.txt")
        
        # Stack data: index, rho/rho_w, u (placeholder), T/T_w
        indices = np.arange(len(x_mesh))
        data = np.column_stack((indices, rho_rel_w, mean_velocity, T_rel_w))
        
        header = "Physical quantities at spatial grid points:\n"
        header += "index -- Density -- Mean Velocity -- Temperature\n"
        header += "--------------------------------------------------------"
        
        np.savetxt(file_path, data, fmt='%d %.8f %.8f %.8f', header=header)
        print(f"Data successfully saved to: {file_path}")

    return rho_rel_w, T_rel_w

# ==================================================================================================
# Compare plots function
# ==================================================================================================


def compare_ytrehus_numerical(x_mesh,
                              x_lim,
                              density_lim,
                              pressure_lim,
                              temperature_lim,
                              output_path = None):
    """
    Plots comparison between Ytrehus theoretical and numerical solutions.
    """
    
    density_num, _, temperature_num = read_physical_quantities(folder + "physical_quantities.txt")
    density_ytr, _, temperature_ytr = read_physical_quantities(folder + "ytrehus/ytrehus_physical_quantities.txt")

    # Compute series
    pressure_num = density_num * temperature_num
    pressure_ytr = density_ytr * temperature_ytr

    # Percentage errors (relative to numerical), with epsilon to avoid divide-by-zero
    eps = 1e-12
    perc_err_density = np.abs(density_ytr - density_num) / np.maximum(np.abs(density_num), eps) * 100.0
    perc_err_pressure = np.abs(pressure_ytr - pressure_num) / np.maximum(np.abs(pressure_num), eps) * 100.0
    perc_err_temperature = np.abs(temperature_ytr - temperature_num) / np.maximum(np.abs(temperature_num), eps) * 100.0

    max_perc_err_density = float(np.max(perc_err_density[:-10]))
    max_perc_err_pressure = float(np.max(perc_err_pressure[:-10]))
    max_perc_err_temperature = float(np.max(perc_err_temperature[:-10]))

    # Plot
    plt.figure(figsize=(15, 5))
    plt.subplot(1, 3, 1)
    plt.xlim(x_lim)
    plt.ylim(density_lim)
    plt.plot(x_mesh, density_num, label='Numerical', color='blue')
    plt.plot(x_mesh, density_ytr, label='Ytrehus Theoretical', color='red', linestyle='--')
    plt.title('Density Comparison')
    plt.xlabel(rf'$ X_1 \ / \ l_w$')
    plt.ylabel(rf'$ \rho \ / \ \rho_w$')
    plt.legend()
    plt.grid()
    
    plt.subplot(1, 3, 2)
    plt.xlim(x_lim)
    plt.ylim(pressure_lim)
    plt.plot(x_mesh, pressure_num, label='Numerical', color='blue')
    plt.plot(x_mesh, pressure_ytr, label='Ytrehus Theoretical', color='red', linestyle='--')
    plt.title('Pressure Comparison')
    plt.xlabel(rf'$ X_1 \ / \ l_w$')
    plt.ylabel(rf'$ p \ / \ p_w$')
    plt.legend()
    plt.grid()

    plt.subplot(1, 3, 3)
    plt.xlim(x_lim)
    plt.ylim(temperature_lim)
    plt.plot(x_mesh, temperature_num, label='Numerical', color='blue')
    plt.plot(x_mesh, temperature_ytr, label='Ytrehus Theoretical', color='red', linestyle='--')
    plt.title('Temperature Comparison')
    plt.xlabel(rf'$ X_1 \ / \ l_w$')
    plt.ylabel(rf'$ T \ / \ T_w$')
    plt.legend()
    plt.grid()

    if output_path is not None:
        plt.savefig(os.path.join(output_path, "ytrehus/ytrehus_comparison.png"))
        print(f"Comparison plot saved to: {output_path}ytrehus/ytrehus_comparison.png")
        
    print("Maximum percentage errors (relative to numerical):")
    print(f"- Density: {max_perc_err_density:.6f}%")
    print(f"- Pressure: {max_perc_err_pressure:.6f}%")
    print(f"- Temperature: {max_perc_err_temperature:.6f}%")

    return {
        "max_percent_error": {
            "density": max_perc_err_density,
            "pressure": max_perc_err_pressure,
            "temperature": max_perc_err_temperature
        }
    }



# ==================================================================================================
# Example Usage
# ==================================================================================================

if __name__ == "__main__":
    
    from postprocessing.read import read_mesh, read_physical_quantities
    import json 
    
    with open("data/Evap/type1.json", "r", encoding="utf-8") as f:
        problem_data = json.load(f)
        
    folder = problem_data["general"]["saving_folder_name"] + "/"
    x = read_mesh(folder + "space_mesh.txt")
    zeta = read_mesh(folder + "velocity_mesh.txt")
    
    plot_ytrehus_3d(1, show=False, save_path="./")
    
    #M_inf = math.sqrt(6.0 / 5.0) * 0.1
    M_inf = 0.8
    #M_inf = problem_data["physical"]["M_infty"]
    
    S_inf, p_ratio, T_ratio, rho_ratio, beta, phi1, phi2, P, r = compute_equilibrium_phys(M_inf)
    print("\n")
    print("--------------------------------------------")
    print("Ytrehus Theoretical Equilibrium Parameters:")
    print("--------------------------------------------")
    print("M_inf:", M_inf)
    print(f"S_inf: {S_inf}")
    print(f"z_e: {1 / p_ratio}")
    print(f"rho_ratio: {rho_ratio}")
    print(f"T_ratio: {T_ratio}")
    print(f"beta: {beta} \n")
    print(f"p_ratio: {p_ratio}")
    print("--------------------------------------------")
    print(f"r: {r}")
    print(f"P: {P}")
    print(f"phi1: {phi1}")
    print(f"phi2: {phi2}")  
    print("--------------------------------------------")
    
    #evap2
    #x_lim_knu = (0, 100)
    #pressure_lim_knu = (0.15, 0.7)
    #temperature_lim_knu = (0.65, 1.1)
    
    #evap5
    #x_lim_knu = (0, 150)
    #pressure_lim_knu = (0.8, 0.88)#(0.83, 0.88)
    #temperature_lim_knu = (0.958, 0.97)#(0.963, 0.97)
    
    #evap6
    #x_lim_knu = (0,240) #(0, 16)
    #density_lim_knu = (0.3, 0.9) #(0.92, 1.02)
    #pressure_lim_knu = (0.3, 0.65) #(0.40, 0.85) #(0.40, 0.62)
    #temperature_lim_knu = (0.75, 0.9) #(0.8, 1.3) #(0.8, 0.88)
    
    #evap7
    x_lim_knu = (0, 130) 
    density_lim_knu = (0.3, 0.9)
    pressure_lim_knu = (0.1, 0.7)
    temperature_lim_knu = (0.65, 0.95)
    
    #evap8
    #x_lim_knu = (0, 100)
    #pressure_lim_knu = (0.6, 1.1)
    #temperature_lim_knu = (0.8, 1.1)
    
    #compute_physical_quantities(M_inf, x, output_path=folder + "ytrehus/")
    #compute_gh(M_inf, x, zeta, output_path=folder + "ytrehus/")
    #compare_ytrehus_numerical(x, x_lim_knu, density_lim_knu, pressure_lim_knu, temperature_lim_knu, output_path=folder)