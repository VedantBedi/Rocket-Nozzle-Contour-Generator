import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import interp1d

# Constants
gamma = 1.4   # Specific heat 
R = 287.0     # Gas constant (J/(kg·K))

def area_ratio(M):
    #Isentropic area ratio A/A* for given Mach M
    term = (2.0/(gamma+1)) * (1 + (gamma-1)/2 * M**2)
    exp = (gamma+1)/(2*(gamma-1))
    return (1.0/M) * (term ** exp)

def mach_from_area_ratio(AR_value, subsonic=False):
    #Mach number from area ratio
    if abs(AR_value - 1.0) < 1e-6:
        return 0.99 if subsonic else 1.01
    M_low, M_high = (0.01, 0.99) if subsonic else (1.01, 10.0)
    for _ in range(50):
        M_mid = (M_low + M_high) / 2.0
        AR_mid = area_ratio(M_mid)
        if abs(AR_mid - AR_value) < 1e-6:
            return M_mid
        if subsonic:
            if AR_mid > AR_value:
                M_low = M_mid
            else:
                M_high = M_mid
        else:
            if AR_mid > AR_value:
                M_high = M_mid
            else:
                M_low = M_mid
    return (M_low + M_high) / 2.0

def create_nozzle_profile(Me, R_throat, L_con=None, L_div=None, R_chamber=None):
    #nozzle profile (x, radius)
    L_con = L_con or 2.5 * R_throat
    L_div = L_div or 6.0 * R_throat
    R_chamber = R_chamber or 3.0 * R_throat
    AR_exit = area_ratio(Me)
    R_exit = R_throat * np.sqrt(AR_exit)
    #Convergent section
    a = 2.0 * (R_chamber - R_throat) / (L_con**3)
    b = 3.0 * (R_chamber - R_throat) / (L_con**2)
    x_con = np.linspace(-L_con, 0, 100)
    r_con = a * x_con**3 + b * x_con**2 + R_throat
    #Divergent section
    x_div = np.linspace(0, L_div, 200)
    r_div = R_throat + (R_exit - R_throat) * np.sin((np.pi/2) * (x_div / L_div))**2
    x_full = np.concatenate((x_con, x_div[1:]))
    r_full = np.concatenate((r_con, r_div[1:]))
    return x_full, r_full

def calculate_mach_distribution(x, r, R_throat):
    #Mach number
    M_vals = np.zeros_like(x)
    for i in range(len(x)):
        AR_local = (r[i] / R_throat)**2
        if x[i] < 0:
            M_vals[i] = mach_from_area_ratio(AR_local, subsonic=True)
        elif abs(x[i]) < 1e-6:
            M_vals[i] = 1.0
        else:
            M_vals[i] = mach_from_area_ratio(AR_local, subsonic=False)
    return M_vals

def plot_nozzle(x, r, M, Me):
    #Plot the nozzle wall Mach number
    x_min, x_max = np.min(x), np.max(x)
    r_max = np.max(r)
    res_x, res_y = 300, 150
    x_grid = np.linspace(x_min, x_max, res_x)
    y_grid = np.linspace(-r_max, r_max, res_y)
    X, Y = np.meshgrid(x_grid, y_grid)
    wall_interp = interp1d(x, r, kind='cubic', bounds_error=False, fill_value="extrapolate")
    mach_interp = interp1d(x, M, kind='cubic', bounds_error=False, fill_value=np.nan)
    Mach_grid = np.full(X.shape, np.nan)
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            xv = X[i, j]
            yv = Y[i, j]
            if abs(yv) <= wall_interp(xv):
                Mach_grid[i, j] = mach_interp(xv)
    cmap = LinearSegmentedColormap.from_list("MachMap", ["green", "yellow", "orange", "red"])
    plt.figure(figsize=(14, 6))
    cp = plt.contourf(X, Y, Mach_grid, levels=1000, cmap=cmap, extend='both')
    plt.plot(x, r, "k-", lw=2, label="Nozzle Contour")
    plt.plot(x, -r, "k-", lw=2)
    plt.fill_between(x, r, -r, color="lightblue", alpha=0.3)
    plt.xlabel("Length (mm)")
    plt.ylabel("Radius (mm)")
    plt.title("Nozzle".format(Me))
    plt.axis("equal")
    plt.colorbar(cp, label="Mach Number")
    plt.legend()
    plt.tight_layout()
    plt.show()

def export_contour_to_xyz(x, r, filename="nozzle_contour_points.csv"):
    #Export nozzle contour points (x, y, 0) in mm.
    points = np.column_stack((x, r, np.zeros_like(x)))
    np.savetxt(filename, points, delimiter=",", header="x,y,z", comments='', fmt="%.6f")
    print(f"Exported contour points to {filename}")

def main():
    #inputs, nozzle profile, plot, CSV.
    try:
        P0 = float(input("Chamber Pressure P0 (Pa 1e6-5e6): "))
        T0 = float(input("Chamber Temperature T0 (K 2000-3500): "))
        mdot = float(input("Mass Flow Rate mdot (kg/s 5-250): "))
        Me = float(input("Exit Mach Number Me (2.5-5.0): "))
    except Exception:
        print("Invalid input; using defaults.")
        P0, T0, mdot, Me = 2e6, 3000, 100, 3.5

    mdot_star = P0 * np.sqrt(gamma/(R*T0)) * ((gamma+1)/2) ** (-(gamma+1)/(2*(gamma-1)))
    A_throat = mdot / mdot_star
    R_throat = np.sqrt(A_throat / np.pi)
    x_nozzle, r_nozzle = create_nozzle_profile(Me, R_throat)
    M_distribution = calculate_mach_distribution(x_nozzle, r_nozzle, R_throat)
    
    # Convert profile from meters to millimeters
    x_mm = x_nozzle * 1000
    r_mm = r_nozzle * 1000

    plot_nozzle(x_mm, r_mm, M_distribution, Me)
    export_contour_to_xyz(x_mm, r_mm)

if __name__ == "__main__":
    main()
