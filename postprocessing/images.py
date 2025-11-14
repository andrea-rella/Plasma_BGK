from operator import neg
import postprocessing.plot as pl
import postprocessing.read as rd
import json


with open("data/fig4.json", "r", encoding="utf-8") as f:
    problem_data = json.load(f)

folder = "output/first_test/"

# ============== EVOLUTION ==============

# timesteps = [1, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]
timesteps = [0,1,20,40,60,80,100, 200, 400, 540,700]
# timesteps = [1, 100, 200]
dt = 0.05

x = rd.read_space_mesh(folder + "space_mesh.txt")

xlim = (0, 1.2)
temp_y_lim = (0.98, 1.03)
velocity_y_lim = (1.16, 1.25)
pressure_y_lim = (16.2, 17.8)

pl.draw_temperature_evolution(
    x, timesteps, dt, folder, xlims=xlim, ylims=temp_y_lim, save_path=folder + "temperature_evolution.png")

pl.draw_velocity_evolution(
    x, timesteps, problem_data["physical"]["T_infty_w"], dt, folder, xlims=xlim, ylims=velocity_y_lim, save_path=folder + "velocity_evolution.png")

pl.draw_pressure_evolution(
    x, timesteps, dt, folder, xlims=xlim, ylims=pressure_y_lim, save_path=folder + "pressure_evolution.png")

print("Plots saved.")
