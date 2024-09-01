from Equations import (
    Structure,
    gibbs_minimizer,
)  # Import custom classes and functions for thermodynamic calculations
import Calculate_PCI  # Module for calculating Pressure-Composition Isotherms
import numpy as np  # For numerical operations
import matplotlib.pyplot as plt  # For plotting


# Initialize structure with given solids and their properties --------------------------------
solids = [[["LA"], ["NI"]], [["CE"], ["NI"]]]  # Define solids involved in the structure
color = ["blue", "orange"]
for n, solid in enumerate(solids):
    multiplicity = [1, 5, 1, 6]  # Define the multiplicity of sites in the structure
    site_fractions = [[1], [1]]  # Initial site fractions for metal atom sites
    Initial_Structure = Structure(
        solid, multiplicity, site_fractions
    )  # Create a structure instance
    T = 273  # Set temperature (in Kelvin)

    # Perform Gibbs energy minimization to find equilibrium compositions --------------------------
    Calculate_PCI.gm, Calculate_PCI.x_H, fractions = gibbs_minimizer(
        T, Initial_Structure
        )  # Minimize Gibbs energy for the structure at temperature T
    HV_fractions = [[i[0] for i in fractions], [i[1] for i in fractions]]

    # Plot Gibbs Energy vs. Hydrogen Mole Fraction -------------------------------------------------
    plt.figure(1)
    plt.plot(
        Calculate_PCI.x_H,
        Calculate_PCI.gm,
        label=f"{solid[0][0]}Ni$_5$-H",
        color=color[n],
    )
    plt.ylabel("Gibbs Energy, $\t{\Delta G}$ / J/mol")
    plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
    plt.legend()
    plt.tight_layout()  # Adjust layout
    # Plot site_fractions vs. Hydrogen Mole Fraction -------------------------------------------------
    plt.figure(2)
    for i, site_fraction in enumerate(HV_fractions):
        if i == 0:
            plt.plot(
                Calculate_PCI.x_H,
                site_fraction,
                "--",
                label=f"{solid[0][0]}Ni$_5$-H: 2b site",
                color=color[n],
            )
        else:
            plt.plot(
                Calculate_PCI.x_H,
                site_fraction,
                label=f"{solid[0][0]}Ni$_5$-H: 6c site",
                color=color[n],
            )
    plt.ylabel("site occupancy of interstitial sites, $\t{y}$ / ")
    plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
    plt.tight_layout()  # Adjust layout
    plt.legend()

    # Calculate Plateau Pressures and related properties for Pressure-Composition Isotherms (PCIs) -
    BOUND_MAX = Calculate_PCI.x_H[-1]  # Get the maximum hydrogen mole fraction

    CTs, plateau_slopes, plateau_intercepts, plateau_pressures = Calculate_PCI.plateau(
        T
    )  # Calculate plateau pressures and related properties at temperature T

    # Plotting Pressure vs. Hydrogen Mole Fraction in a semilog plot -------------------------------
    pl_root = np.log10(
        plateau_pressures[0]
    )  # Calculate the logarithmic base for plateau pressures
    p_plot = np.logspace(
        pl_root - 2, pl_root + 2, 20
    )  # Generate pressures for plotting
    # Adjust pressure points for detailed plotting around plateau pressures
    for k in plateau_pressures:
        p_plot = np.append(p_plot, k + k / 200000000)
        p_plot = np.append(p_plot, k - k / 200000000)
    p_plot = np.sort(p_plot)  # Sort the pressures for plotting
    x_plot = []  # Initialize list for hydrogen mole fractions
    # Calculate hydrogen mole fraction for each pressure point
    for p in p_plot:
        result, _, _ = Calculate_PCI.calculate_x_from_p(p, plateau_pressures, CTs, T)
        x_plot.append(Calculate_PCI.x_H[result])
    plt.figure(3)
    plt.semilogy(
        x_plot, p_plot, label=f"{solid[0][0]}Ni$_5$-H", color=color[n]
    )  # Create a semilogarithmic plot
    plt.ylabel("Hydrogen Pressure, $\t{P}$ / Pa")
    plt.xlabel("Hydrogen mole fraction, $\t{x_{H}}$ / ")
    plt.tight_layout()  # Adjust
    plt.legend()
    plt.grid()
plt.figure(1)
plt.savefig('plot_gibbs.png')  
plt.figure(2)
plt.savefig('plot_sf.png')
plt.figure(3)
plt.savefig('plot_PCI.png')
