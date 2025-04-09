# -----------------------------
# Mesh Definition: Defines the simulation domain.
# -----------------------------
[Mesh]
    type = GeneratedMesh           # Use a generated 1D mesh.
    dim = 1                        # 1D simulation.
    nx = 10                        # Divide the 1 mm thick sample into 10 elements.
    xmin = 0.0                   # Start position (m).
    xmax = 1.0e-3                  # End position (m) => 1 mm thickness.
    elem_type = EDGE2              # 2-node line elements for 1D.
[]
  
# -----------------------------
# Global Parameters: Sets overall simulation parameters.
# -----------------------------
[GlobalParams]
    temperature = 295              # Temperature in Kelvin (room temperature: 22°C).
    material = W                   # Use tungsten material properties.
[]
  
# -----------------------------
# Variables: Declare the diffusing species.
# -----------------------------
[Variables]
    [tritium]
        # This variable will store the tritium concentration profile.
    []
[]
  
# -----------------------------
# Initial Conditions (ICs): Specify the starting state.
# -----------------------------
[ICs]
    [tritium_ic]
        type = ConstantIC            # Initialize with a constant value.
        variable = tritium           # Applies to the tritium variable.
        value = 0.0                  # Start with zero tritium inside the sample.
    []
[]
  
# -----------------------------
# Kernels: Define the governing equations (here, diffusion).
# -----------------------------
[Kernels]
    [diff]
        type = LMDiffusion              # The diffusion kernel for tritium.
        primal_variable = tritium       # Variable from which the Laplacian is computer
        variable = tritium              # Refers to the quantity being solved for, tritium concentration in this case
        diffusivity = 5.6e-7
        block = 0
    []
[]
  
# -----------------------------
# Boundary Conditions (BCs): Define how the boundaries behave.
# -----------------------------
[BCs]
    [left]
        type = DirichletBC           # Fixed-value boundary condition.
        variable = tritium           # Applies to the tritium variable.
        boundary = left              # On the left side of the sample.
        value = 0.0                  # Tritium concentration fixed at zero.
    []
    [right]                      # This boundary simulates the incoming plasma flux.
        type = DiffusionFluxBC            # Imposes a flux boundary condition.
        variable = tritium           # Applies to the tritium variable.
        boundary = right             # On the right side of the sample.
        value = 4.8e21               # Incoming flux of tritium (atoms/m^2/s).
    []
[]
  
# -----------------------------
# Materials: Define material properties for tungsten.
# -----------------------------
# [Materials]
#    
# []
  
# -----------------------------
# Executioner: Set up the transient simulation parameters.
# -----------------------------
[Executioner]
    type = Transient               # Transient (time-dependent) simulation.
    scheme = bdf2                  # 2nd-order backward differentiation formula.
    dt = 100.0                     # Time step size in seconds.
    end_time = 2500                # Total simulation time (in seconds).
    solve_type = 'NEWTON'          # Use a Newton solver.
[]
  
# -----------------------------
# Outputs: Configure simulation outputs.
# -----------------------------
[Outputs]
    exodus = true                # Generate an Exodus file for visualization.
    csv = true                   # Output data in CSV format.
    print_linear_residuals = false  # Option to print linear residuals (disabled here).
[]
  