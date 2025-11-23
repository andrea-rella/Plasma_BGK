import postprocessing.plot as pl
import postprocessing.read as rd
import json
from math import sqrt

with open("data/fig9f.json", "r", encoding="utf-8") as f:
    problem_data = json.load(f)

folder = "output/" + problem_data["general"]["saving_folder_name"] + "/"

# ============== EVOLUTION ==============

# fig4a 
#timesteps = [0,1,20,40,60,80,100, 200, 400, 540,700]
# fig4b
# timesteps = [0,1,40,100, 250, 750, 1500, 2500, 3500, 4500, 5000]
# fig4c
#timesteps = [0, 20000, 40000, 60000, 80000, 100000, 120000, 140000, 160000, 180000, 200000, 220000]
# fig4c_test
#timesteps = [50000]
# fig4c_test2
#timesteps = [2000, 8000, 14000, 20000, 25000]
# fig9a
#timesteps = [0,1, 20, 100, 250, 750]

# fig9c
#timesteps = [0,1, 60, 100, 500, 1000, 2000]

# fig9d
# timesteps = [0,1, 60, 100, 500, 1000, 3000, 6500]

# fig9e
# timesteps = [0, 100, 1000, 3000, 5000, 8000, 10000]

# fig9f
timesteps = [0, 100, 1000, 5000]

dt = 0.05

# fig4a
#xlim = (0, 1.2)
#temp_y_lim = (0.98, 1.03)
#velocity_y_lim = (1.16, 1.25)
#pressure_y_lim = (16.2, 17.8)

# fig4b
#xlim = (0, 6)
#temp_y_lim = (0.987, 1.01)
#velocity_y_lim = (1.025, 1.052)
#pressure_y_lim = (16.5, 17.1)

# fig4c
# xlim = (0, 180)
# temp_y_lim = (0.968, 1.005)
# velocity_y_lim = (0.96, 1.02)
# pressure_y_lim = (14.8, 16.2)
# mach_y_lim = (0.96, 1.04)
# target_Machs=np.linspace(0.98, 1.02, 100)
# Vp_y_lim = (0.0, 0.02)

# fig9a
# xlim = (0, 1.5)
# temp_y_lim = (0.45, 0.65)
# velocity_y_lim = (0.4, 0.8)
# pressure_y_lim = (2.4, 3.4)
#
# t_horizontal = 0.512201
# v_horizontal = float(- (-0.331606) * (1 / sqrt((5.0 / 6.0) * problem_data["physical"]["T_infty_w"])))
# p_horizontal = float(6.21985 * t_horizontal)

# fig9c
# xlim = (0, 6)
# temp_y_lim = (1.65, 2.2)
# velocity_y_lim = (0.46, 0.57)
# pressure_y_lim = (3.0, 3.4)
# 
# t_horizontal = 2.0194
# v_horizontal = float(- (-0.691114) * (1 / sqrt((5.0 / 6.0) * problem_data["physical"]["T_infty_w"])))
# p_horizontal = float(1.62351 * t_horizontal)

# fig9d
# xlim = (0, 18)
# temp_y_lim = (2.4, 4.2)
# velocity_y_lim = (0.25, 0.55)
# pressure_y_lim = (3.0, 3.4)
# 
# t_horizontal = 4.11163
# v_horizontal = float(- (-0.926062) * (1 / sqrt((5.0 / 6.0) * problem_data["physical"]["T_infty_w"])))
# p_horizontal = float(0.781162 * t_horizontal)

# fig9e
# xlim = (0, 30)
# temp_y_lim = (0.98, 1.04)
# velocity_y_lim = (0.092, 0.104)
# pressure_y_lim = (1.173, 1.22)
# 
# t_horizontal = 1.00301
# v_horizontal = float(- (-0.0872349) * (1 / sqrt((5.0 / 6.0) * problem_data["physical"]["T_infty_w"])))
# p_horizontal = float(1.20522 * t_horizontal)

# fig9f
xlim = (0, 0.9)
temp_y_lim = (0.94, 1.01)
velocity_y_lim = (0.88, 1.04)
pressure_y_lim = (8, 11)

t_horizontal = 0.99389
v_horizontal = float(- (-0.829854) * (1 / sqrt((5.0 / 6.0) * problem_data["physical"]["T_infty_w"])))
p_horizontal = float(9.91186 * t_horizontal)


x = rd.read_space_mesh(folder + "space_mesh.txt")

pl.draw_temperature_evolution(
    x, 
    timesteps, 
    dt, 
    folder, 
    xlims=xlim, 
    ylims=temp_y_lim, 
    save_path=folder + "temperature_evolution.png"
    )

pl.draw_velocity_evolution(
    x, 
    timesteps, 
    problem_data["physical"]["T_infty_w"], 
    dt, 
    folder, 
    xlims=xlim, 
    ylims=velocity_y_lim, 
    save_path=folder + "velocity_evolution.png"
    )

pl.draw_pressure_evolution(
    x, 
    timesteps, 
    dt, 
    folder, 
    xlims=xlim, 
    ylims=pressure_y_lim, 
    save_path=folder + "pressure_evolution.png"
    )
    
# pl.draw_MachNumber_evolution(
#     x, 
#     timesteps, 
#     dt, 
#     folder, 
#     xlims=xlim, 
#     ylims=mach_y_lim, 
#     save_path=folder + "mach_evolution.png"
#     )

# pl.draw_Vp_star_evolution(
#     x = x, 
#     timesteps = timesteps,
#     m = 1000,
#     dt = dt, 
#     folder = folder, 
#     target_Machs=target_Machs, 
#     xlims=(0.98, 1.02), 
#     ylims=Vp_y_lim, 
#     save_path=folder + "Vp_star_evolution.png"
#     )

# ============== FINAL STATE ==============

pl.draw_temperature_profile(x, 
                            folder, 
                            xlims=xlim, 
                            ylims=temp_y_lim,
                            horizontal_line=t_horizontal,
                            save_path=folder + "final_temperature.png")

pl.draw_velocity_profile(x,  
                         folder, 
                         problem_data["physical"]["T_infty_w"],
                         xlims=xlim, 
                         ylims=velocity_y_lim,
                         horizontal_line=v_horizontal,
                         save_path=folder + "final_velocity.png")

pl.draw_pressure_profile(x, 
                        folder, 
                        xlims=xlim, 
                        ylims=pressure_y_lim,
                        horizontal_line=p_horizontal,
                        save_path=folder + "final_pressure.png")

print("Plots saved.")