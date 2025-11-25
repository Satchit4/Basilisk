# Basilisk Bubble Simulations

This repository contains Basilisk C setups for simulating bubble dynamics in different configurations.

## Project Overview

The project is divided into two main simulation types:

1.  **2D Bubble Simulation**: Located in `bubble 2D/`. This setup simulates a bubble in a 2D planar domain, likely focusing on interaction with boundaries or flow conditions.
2.  **Axisymmetric Bubble Simulation**: Located in `axi/axi/`. This setup uses axisymmetric coordinates to simulate a bubble, which is useful for representing 3D bubbles with rotational symmetry (e.g., rising in a tube).

## Directory Structure

-   **`bubble 2D/`**: Contains the source code for the 2D simulation.
    -   `bubble.c`: Main simulation file.
    -   `iforce_palas.h`: Header for interfacial forces.
    -   `tension_palas.h`: Header for surface tension.
    -   `output_vtu_foreach.h`: Helper for VTU output.
-   **`axi/`**: Contains the axisymmetric simulation files.
    -   `axi/bubble.c`: Main simulation file for the axisymmetric case.
    -   `axi/output_vtu_foreach.h`: Helper for VTU output.

## Prerequisites

-   **Basilisk C**: You need to have Basilisk installed and the `qcc` compiler available in your path.
-   **MPI** (Optional): For running simulations in parallel.
-   **Paraview** (Optional): For visualizing the `.vtu` output files.

## Usage

### Compiling and Running

Navigate to the directory of the simulation you want to run.

**Example for 2D Bubble:**

```bash
cd "bubble 2D"
# Compile
qcc -O2 -Wall bubble.c -o bubble -lm
# Run (requires an argument for MAXLEVEL, e.g., 10)
./bubble 10
```

**Example for Axisymmetric Bubble:**

```bash
cd axi/axi
# Compile
qcc -O2 -Wall bubble.c -o bubble -lm
# Run (requires an argument for MAXLEVEL, e.g., 10)
./bubble 10
```

### Parallel Execution

If compiled with MPI support (using `CC99='mpicc -std=c99' qcc ... -D_MPI=1`):

```bash
mpirun -np 4 ./bubble 10
```

### Visualization

The simulations generate `.vtu` or `.pvtu` files (e.g., `bessel-*.vtu`, `sol_*.vtu`). Open these files in Paraview to visualize the fields (Volume fraction `f`, Velocity `u`, Pressure `p`).

## Simulation Details

-   **Physical Parameters**: The simulations use CGS units (or consistent non-dimensional units). Check the `bubble.c` files for specific values of density (`rho1`, `rho2`), viscosity (`mu1`, `mu2`), and surface tension (`f.sigma`).
-   **Geometry**:
    -   The 2D case involves specific boundary conditions for inflow/outflow.
    -   The axisymmetric case simulates a bubble in a channel/tube configuration.
