# Attempting a simplified version of the divertor_monoblock example

### Nomenclatures
### C_mobile_j      mobile H concentration in "j" material, where j = CuCrZr, Cu, W
### C_trapped_j     trapped H concentration in "j" material, where j = CuCrZr, Cu, W
### C_total_j       total H concentration in "j" material, where j = CuCrZr, Cu, W
###
### S_empty_j       empty site concentration in "j" material, where j = CuCrZr, Cu, W
### S_trapped_j     trapped site concentration in "j" material, where j = CuCrZr, Cu, W
### S_total_j       total site H concentration in "j" material, where j = CuCrZr, Cu, W
###
### F_permeation    permeation flux
### F_recombination recombination flux
###
### Sc_             Scaled
### Int_            Integrated
### ScInt_          Scaled and integrated

tungsten_atomic_density = ${units 6.338e28 m^-3}

# Can use a simple rectangle mesh, simulating that you're looking at the side of the button
# and the plasma is impacting normal to the surface from the right and the button is stood
# up on the edge
[Mesh]
    type = GeneratedMesh           # Use a generated 1D mesh.
    dim = 2                        # 2D simulation.
    nx = 10                        # Divide the 1 mm thick sample into 10 elements.
    ny = 10
    xmin = 0.0                     # Start position (m).
    xmax = 1.0e-3                  # End position (m) => 1 mm thickness.
    ymin = 0.0
    ymax = 6.0e-3
    elem_type = QUAD4              # 2D quadrilateral elements
[]

# Skipping the ReferenceResidualProblem included in monoblock example. Relying on default convergence tests.

[Variables]
    [temperature]
        order = FIRST
        family = LAGRANGE
        initial_condition = ${units 295 K}  #22C room temp
    []
    [C_mobile_W]
        order  = FIRST
        family = LAGRANGE
        # initial_condition = ${units ### m^-3}
        block =  0
    []
    [C_trapped_W]
        order = FIRST
        family = LAGRANGE
        # initial_condition = ${units ### m^-3}
        block = 0
    []
[]


# Used to track quantities that aren't solved for by the PDEs, 
# but are needed to compute material properties, or are desirable to obtain as outputs.
[AuxVariables]      
    [flux_y]
        order = FIRST
        family = MONOMIAL
    []
    [Sc_C_mobile_W]     # Sc = scaled
        block = 0
    []
    [Sc_C_trapped_W]
        block = 0
    []
    [C_total_W]
        block = 0
    []
    [Sc_C_total_W]
        block = 0
    []
    [S_empty_W]
        block = 0
    []
    [Sc_S_empty_W]
        block = 0
    []
    [S_trapped_W]
        block = 0
    []
    [Sc_S_trapped_W]
        block = 0
    []
    [S_total_W]
        block = 0
    []
    [Sc_S_total_W]
        block = 0
    []
[]

[Kernels]
    [diff_W]
        type = ADMatDiffusion
        variable = C_mobile_W
        diffusivity = diffusivity_W
        block = 0
        # extra_vector_tags = ref
    []
    [time_diff_W]
        type = ADTimeDerivative
        variable = C_mobile_W
        block = 0
        # extra_vector_tags = ref
    []
    [coupled_time_W]
        type = ScaledCoupledTimeDerivative
        variable = C_mobile_W
        v = C_trapped_W
        factor = 1e0
        block = 0
        # extra_vector_tags = ref
    []
    [heat_conduction_W]
        type = HeatConduction
        variable = temperature
        diffusion_coefficient = thermal_conductivity_W
        block = 0
        # extra_vector_tags = ref
    []
    [time_heat_conduction_W]
        type = SpecificHeatConductionTimeDerivative
        variable = temperature
        specific_heat = specific_heat_W
        density = density_W
        block = 0
        # extra_vector_tags = ref
    []
[]

# Expressions that define the AuxVariables
[AuxKernels]
    [Scaled_mobile_W]
        variable = Sc_C_mobile_W
        type = NormalizationAux
        normal_factor = ${tungsten_atomic_density}
        source_variable = C_mobile_W
    []
    [Scaled_trapped_W]
        variable = Sc_C_trapped_W
        type = NormalizationAux
        normal_factor = ${tungsten_atomic_density}
        source_variable = C_trapped_W
    []
    [total_W]
        variable = C_total_W
        type = ParsedAux
        expression = 'C_mobile_W + C_trapped_W'
        coupled_variables = 'C_mobile_W C_trapped_W'
    []
    [Scaled_total_W]
        variable = Sc_C_total_W
        type = NormalizationAux
        normal_factor = ${tungsten_atomic_density}
        source_variable = C_total_W
    []
    [empty_sites_W]
        variable = S_empty_W
        type = EmptySitesAux
        N = ${units 1.0e0 m^-3}       # = ${tungsten_atomic_density} #/m^3 (W lattice density)
        # Ct0 = ${units 1.0e-4 m^-3}   # E.A. Hodille et al 2021 Nucl. Fusion 61 126003, trap 1
        Ct0 = ${units 1.0e-4 m^-3}    # E.A. Hodille et al 2021 Nucl. Fusion 61 1260033, trap 2
        trap_per_free = 1.0e0         # 1.0e1
        trapped_concentration_variables = C_trapped_W
    []
    [scaled_empty_W]
        variable = Sc_S_empty_W
        type = NormalizationAux
        normal_factor = ${tungsten_atomic_density}
        source_variable = S_empty_W
    []
    [trapped_sites_W]
        variable = S_trapped_W
        type = NormalizationAux
        normal_factor = 1e0
        source_variable = C_trapped_W
    []
    [scaled_trapped_W]
        variable = Sc_S_trapped_W
        type = NormalizationAux
        normal_factor = ${tungsten_atomic_density}
        source_variable = S_trapped_W
    []
    [total_sites_W]
        variable = S_total_W
        type = ParsedAux
        expression = 'S_trapped_W + S_empty_W'
        coupled_variables = 'S_trapped_W S_empty_W'
    []
    [scaled_total_W]
        variable = Sc_S_total_W
        type = NormalizationAux
        normal_factor = ${tungsten_atomic_density}
        source_variable = S_total_W
    []
    [flux_y_W]
        type = DiffusionFluxAux
        diffusivity = diffusivity_W
        variable = flux_y
        diffusion_variable = C_mobile_W
        component = x
        block = 0
    []
[]

# I think I can skip InterfaceKernels because there are no boundaries in my cases
# ADPenaltyInterfaceDiffusion is used to conserve the particle flux at the interface b/w
# 2 different solubilities.

# NodalKernels are used for non-diffusive variables (trapped sites)
[NodalKernels]
    [time_W]
        type = TimeDerivativeNodalKernel
        variable = C_trapped_W
    []
    [trapping_W]
        type = TrappingNodalKernel
        variable = C_trapped_W
        temperature = temperature
        alpha_t = 2.75e11      # 1e15
        N = 1.0e0  # = (1e0) x (${tungsten_atomic_density} #/m^3)
        # Ct0 = 1.0e-4                # E.A. Hodille et al 2021 Nucl. Fusion 61 126003, trap 1
        Ct0 = 1.0e-4                # E.A. Hodille et al 2021 Nucl. Fusion 61 1260033, trap 2
        trap_per_free = 1.0e0       # 1.0e1
        mobile_concentration = 'C_mobile_W'
        # extra_vector_tags = ref
    []
    [release_W]
        type = ReleasingNodalKernel
        alpha_r = 8.4e12    # 1.0e13
        temperature = temperature
        # detrapping_energy = 9863.9    # = 0.85 eV    E.A. Hodille et al 2021 Nucl. Fusion 61 126003, trap 1
        detrapping_energy = 11604.6   # = 1.00 eV    E.A. Hodille et al 2021 Nucl. Fusion 61 126003, trap 2
        variable = C_trapped_W
    []
[]

[BCs]
    [left]
        type = DirichletBC           # Fixed-value boundary condition.
        variable = C_mobile_W        # Applies to the tritium variable.
        boundary = left              # On the left side of the sample.
        value = 0.0                  # Tritium concentration fixed at zero.
    []
    [right]                      # This boundary simulates the incoming plasma flux.
        type = DiffusionFluxBC            # Imposes a flux boundary condition.
        variable = C_mobile_W           # Applies to the tritium variable.
        boundary = right             # On the right side of the sample.
        value = 4.8e21               # Incoming flux of tritium (atoms/m^2/s).
    []
    [vertical_flux]         # Nothing leaving the mesh from the top or bottom
        type = DiffusionFluxBC
        variable = C_mobile_W
        boundary = ' top bottom'
        value = 0.0
    []
[]

# I think I can skip the functions b/c my exposures are simpler and I'm using diff BCs

[Materials]
    [diffusivity_W]
        type = ADParsedMaterial
        property_name = diffusivity_W
        coupled_variables = 'temperature'
        block = 0
        #expression = '2.4e-7*exp(-4525.8/temperature)'    # H diffusivity in W
        expression = 5.6e-7
        outputs = all
    []
    [solubility_W]
        type = ADParsedMaterial
        property_name = solubility_W
        coupled_variables = 'temperature'
        block = 0
        # expression = '2.95e-5 *exp(-12069.0/temperature)'              # H solubility in W = (1.87e24)/(${tungsten_atomic_density}) [#/m^3]
        expression = '2.95e-5 *exp(-12069.0/temperature) + 4.95e-8 * exp(-6614.6/temperature)'    # H solubility in W = (1.87e24)/(${tungsten_atomic_density}) [#/m^3]
        outputs = all
    []
    [converter_to_regular_W]
        type = MaterialADConverter
        ad_props_in = 'diffusivity_W'
        reg_props_out = 'diffusivity_W_nonAD'
        block = 0
    []
    [heat_transfer_W]
        type = GenericConstantMaterial
        prop_names = 'density_W'
        prop_values = '19300'                # [g/m^3]
        block = 0
    []
    [specific_heat_W]
        type = ParsedMaterial
        property_name = specific_heat_W
        coupled_variables = 'temperature'
        block = 0
        expression = '1.16e2 + 7.11e-2 * temperature - 6.58e-5 * temperature^2 + 3.24e-8 * temperature^3 -5.45e-12 * temperature^4'    # ~ 132[J/kg-K]
        outputs = all
    []
    [thermal_conductivity_W]
        type = ParsedMaterial
        property_name = thermal_conductivity_W
        coupled_variables = 'temperature'
        block = 0
        # expression = '-7.8e-9 * temperature^3 + 5.0e-5 * temperature^2 - 1.1e-1 * temperature + 1.8e2'    # ~ 173.0 [ W/m-K]   from R. Delaporte-Mathurin et al 2021 Nucl. Fusion 61 036038,
        expression = '2.41e2 - 2.90e-1 * temperature + 2.54e-4 * temperature^2 - 1.03e-7 * temperature^3 + 1.52e-11 * temperature^4'    # ~ 173.0 [ W/m-K]
        outputs = all
    []
[]

[Postprocessors]
    [F_recombination]
        type = SideDiffusiveFluxAverage
        boundary = 'right'
        diffusivity = 5.01e-24   # (3.01604928)/(6.02e23)/[gram(T)/m^2]
        # diffusivity = 5.508e-19   # (1.0e3)*(1.0e3)/(6.02e23)/(3.01604928) [gram(T)/m^2]
        variable = Sc_C_total_W
    []
    [Int_C_mobile_W]
        type = ElementIntegralVariablePostprocessor
        variable = C_mobile_W
        block = 0
    []
    [ScInt_C_mobile_W]
        type = ScalePostprocessor
        value =  Int_C_mobile_W
        scaling_factor = 3.491e10   # (1.0e3)*(1.0e3)*(${tungsten_atomic_density})/(6.02e23)/(3.01604928) [gram(T)/m^2]
    []
    [Int_C_trapped_W]
        type = ElementIntegralVariablePostprocessor
        variable = C_trapped_W
        block = 0
    []
    [ScInt_C_trapped_W]
        type = ScalePostprocessor
        value = Int_C_trapped_W
        scaling_factor = 3.491e10   # (1.0e3)*(1.0e3)*(${tungsten_atomic_density})/(6.02e23)/(3.01604928) [gram(T)/m^2]
    []
    [Int_C_total_W]
        type = ElementIntegralVariablePostprocessor
        variable = C_total_W
        block = 0
    []
    [ScInt_C_total_W]
        type = ScalePostprocessor
        value = Int_C_total_W
        scaling_factor = 3.491e10   # (1.0e3)*(1.0e3)*(${tungsten_atomic_density})/(6.02e23)/(3.01604928) [gram(T)/m^2]
    []
    [dt]
        type = TimestepSize
    []
    #[temperature_top]
    #    type = PointValue
    #    variable = temperature
    #    point = '0 14.0e-3 0'
    #[]
    # limit timestep
    #[timestep_max_pp] # s
    #type = FunctionValuePostprocessor
    #function = timestep_function
    #[]
[]

# Skipping VectorPostprocessors for now
# ^Same with Preconditioning

[Executioner]
#    type = Transient
#    scheme = bdf2       # 2nd order backward differentiation formula time integration scheme
#    solve_type = NEWTON
#    petsc_options_iname = '-pc_type'
#    petsc_options_value = 'lu'
#    nl_rel_tol  = 1e-6 # nonlinear relative tolerance
#    nl_abs_tol  = 1e-7 # nonlinear absolute tolerance
#    end_time = 2500     # sec
#    automatic_scaling = true
#    line_search = 'none'
#    dtmin = 1e-4    # Min timestep size in an adaptive run
#    nl_max_its = 18 # max nonlinear iterations
#    [TimeStepper]
#        type = IterationAdaptiveDT
#        dt = 20
#        optimal_iterations = 15
#        iteration_window = 1
#        growth_factor = 1.2
#        cutback_factor = 0.8
#        timestep_limiting_postprocessor = timestep_max_pp
#    []
    type = Transient               # Transient (time-dependent) simulation.
    scheme = bdf2                  # 2nd-order backward differentiation formula.
    dt = 100.0                     # Time step size in seconds.
    end_time = 2500                # Total simulation time (in seconds).
    solve_type = 'NEWTON'          # Use a Newton solver.
[]

[Outputs]
    [exodus]
        type = Exodus
        sync_only = false
        # output at key moment in the first two cycles, and then at the end of the simulation
        sync_times = '110.0 480.0 590.0 1600.0 1710.0 2080.0 2190.0 3400.0 8.0e4'
    []
    csv = true
    perf_graph = true
[]