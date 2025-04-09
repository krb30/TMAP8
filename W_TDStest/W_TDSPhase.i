# -----------------------------
# Mesh Definition
# -----------------------------
[Mesh]
    type = GeneratedMesh           # Simple 1D mesh generation
    dim = 1                        # 1D geometry (x-direction only)
    nx = 10                        # Number of elements
    xmin = 0.0
    xmax = 1.0e-3                  # 1 mm sample thickness
    elem_type = EDGE2              # 2-node 1D elements
[]
  
# -----------------------------
# Global Parameters
# -----------------------------
[GlobalParams]
    temperature = 295              # Initial temperature (used as fallback)
    material = W                   # Tungsten material properties
[]

# -----------------------------
# Auxiliary Variable for Temperature Ramp
# -----------------------------
[AuxVariables]
    [temp]
        order = FIRST                # First-order interpolation
        family = LAGRANGE            # Lagrange polynomial basis
    []
[]
  
# -----------------------------
# Apply a Time-Dependent Temperature Function
# -----------------------------
[AuxKernels]
    [tds_temp]
        type = FunctionAux           # Assigns the temp variable a value using a function
        variable = temp
        function = ramp_temp         # References the ramp function defined below
    []
[]
  
# -----------------------------
# Define the Temperature Ramp Function
# -----------------------------
[Functions]
    [ramp_temp]
        type = PiecewiseLinear
        x = '0   1000'               # Time in seconds
        y = '295 1000'               # Temperature in Kelvin
                                    # Linear ramp from 22°C to 727°C over 1000 seconds
    []
[]
  
# -----------------------------
# Main Diffusing Species
# -----------------------------
[Variables]
    [tritium]
    []
[]
  
# -----------------------------
# Initial Conditions (ICs)
# -----------------------------
[ICs]
    [tritium_ic]
        type = FileIC                # Load initial tritium profile from exposure simulation
        variable = tritium
        file = W_ExposurePhase_out.e # Exodus output from W_ExposurePhase.i run
        timeset = 0                  # Use the final time step from that simulation
    []
[]
  
# -----------------------------
# Physics Kernels
# -----------------------------
[Kernels]
    [diff]
        type = Diffusion         # Fickian diffusion kernel
        variable = tritium
        temperature = temp           # Use the ramped temperature
    []
[]
  
# -----------------------------
# Boundary Conditions (Allow Desorption)
# -----------------------------
[BCs]
    [left]
        type = DirichletBC
        variable = tritium
        boundary = left
        value = 0.0                  # Zero tritium concentration at left boundary
                                    # Simulates tritium freely escaping (desorbing)
    []
    
    [right]
        type = DirichletBC
        variable = tritium
        boundary = right
        value = 0.0                  # Same for right boundary
    []
[]
  
# -----------------------------
# Material Model (W Properties with Temperature Dependence)
# -----------------------------
[Materials]
    [w_props]
        type = TMAPMaterial
        block = 0
        temperature = temp           # Passes the ramped temperature to material model
    []
[]
  
# -----------------------------
# Time Stepping and Solver Settings
# -----------------------------
[Executioner]
    type = Transient
    scheme = bdf2                  # 2nd-order backward differentiation
    dt = 1.0                       # Time step size in seconds
    end_time = 1000.0              # Total simulation time (match the heating ramp)
    solve_type = 'NEWTON'
[]
  
# -----------------------------
# Outputs (Visualization and Data)
# -----------------------------
[Outputs]
    exodus = true                  # Create .e file for post-processing in Paraview
    csv = true                     # Save CSV output for easier plotting
    print_linear_residuals = false
[]
  