# General parameters
kB = '${units 1.380649e-23 J/K}' # Boltzmann constant (from PhysicalConstants.h - https://physics.nist.gov/cgi-bin/cuu/Value?r)

# Model parameters
TDS_initial_time = '${units 5.1e3 s}'
TDS_critial_time_1 = '${units 5400 s}'
TDS_critial_time_2 = '${units 5404 s}'
simulation_time = '${units 7.3e3 s}'   # 5100 s implantation + ~1946 s ramp to 1273 K + margin
outputs_initial_time = '${units 0 s}'
step_interval_max = 50 # (-)
step_interval_mid = 15 # (-)
step_interval_min = 6 # (-)
#bound_value_max = '${units 2e4 at/mum^3}'
#bound_value_min = 0.0 #'${units -1e-10 at/mum^3}'

# Diffusion parameters
flux_high = '${units 1.8e17 at/m^2/s -> at/mum^2/s}'    # incident flux (2e21) * implantation ratio (9e-5)
flux_low = '${units 0      at/mum^2/s}'
diffusivity_coefficient = '${units 5.05e-7 m^2/s -> mum^2/s}'   # 4.1e-7 * sqrt(3/2)
E_D = '${units 0.39 eV -> J}'
initial_concentration = '${units 1e-10 at/m^3 -> at/mum^3}'
width_source = '${units 3e-9 m -> mum}'
depth_source = '${units 4.6e-9 m -> mum}'

# Traps parameters
N = '${units 6.25e28 at/m^3 -> at/mum^3}'
#initial_concentration_trap_2 = 4.4e-10 # (-)
#initial_concentration_trap_3 = 1.4e-10 # (-)
trapping_energy = '${fparse ${units 0.39 eV -> J} / ${kB}}'
#detrapping_energy_1 = '${fparse ${units 1.2 eV -> J} / ${kB}}'
#detrapping_energy_2 = '${fparse ${units 1.6 eV -> J} / ${kB}}'
detrapping_energy_3 = '${fparse ${units 1.73 eV -> J} / ${kB}}'
#trapping_site_fraction_1 = 0.002156 # (-)
#trapping_site_fraction_2 = 0.00175 # (-)
#trapping_site_fraction_3 = 0.0020 # (-)
trapping_rate_prefactor = '${units 1.119e13 1/s}'    # 9.1316e12 * sqrt(3/2)
release_rate_profactor  = '${units 1.030e13 1/s}'    # 8.4e12 * sqrt(3/2)
#trap_per_free_1 = 1e6 # (-)
#trap_per_free_2 = 1e4 # (-)
trap_per_free_3 = 1e4 # (-)
#width_trap1 = '${units 10e-9 m -> mum}'
a3 = 4.6e-5
b3 = 1.5e-1
x0 = ${units 60e-9 m -> mum}

# thermal parameters
temperature_low = '${units 600 K}'  # temp during ion implantation
temperature_TDS_start = '${units 300 K}'
temperature_cool_duration = '${units 60 s}' 
temperature_high = '${units 1273 K}'
temperature_rate = '${units 0.5 K/s}'

Kr0 = '${units 3.2e-15 m^4/s -> mum^4/s}'
Ec  = '${fparse ${units 1.16 eV -> J}}'

[Mesh]
  active = 'cartesian_mesh'
  [cartesian_mesh]
    nx_scale = 2
    type = CartesianMeshGenerator
    dim = 1
    dx = '${fparse 10 * ${units 1.5e-9 m -> mum}}
          ${units 1e-9 m -> mum}    ${units 1e-8 m -> mum}    ${units 1e-7 m -> mum}
          ${units 4e-6 m -> mum}    ${units 3.9956e-3 m -> mum}'
    subdomain_id = '0 1 1 1 1 1'
    ix = '${fparse 10 * ${nx_scale}}
          ${fparse 1 * ${nx_scale}}    ${fparse 1 * ${nx_scale}}   ${fparse 1 * ${nx_scale}}
          ${fparse 6 * ${nx_scale}}    50'
  []

  [cartesian_mesh_coarse]
    nx_scale = 2
    type = CartesianMeshGenerator
    dim = 1
    dx = '${fparse 10 * ${units 1.5e-9 m -> mum}}
          ${units 1e-9 m -> mum}    ${units 1e-8 m -> mum}    ${units 1e-7 m -> mum}
          ${units 4e-6 m -> mum}    ${units 3.9956e-3 m -> mum}'
    subdomain_id = '0 1 1 1 1 1'
    ix = '${fparse 10 * ${nx_scale}}
          ${fparse 1 * ${nx_scale}}    ${fparse 1 * ${nx_scale}}   ${fparse 1 * ${nx_scale}}
          ${fparse 6 * ${nx_scale}}    50'
  []
[]

[Problem]
  type = ReferenceResidualProblem
  extra_tag_vectors = 'ref'
  reference_vector = 'ref'
[]

#[Bounds]
#  [concentration_lower_bound]
#    type = ConstantBounds
#    variable = bounds_dummy
#    bounded_variable = concentration
#    bound_type = lower
#    bound_value = ${bound_value_min}
#  []
#  [trapped_1_lower_bound]
#    type = ConstantBounds
#    variable = bounds_dummy
#    bounded_variable = trapped_1
#    bound_type = lower
#    bound_value = ${bound_value_min}
#  []
#  [trapped_2_upper_bound]
#    type = ConstantBounds
#    variable = bounds_dummy
#    bounded_variable = trapped_2
#    bound_type = upper
#    bound_value = ${bound_value_max}
#  []
#  [trapped_2_lower_bound]
#    type = ConstantBounds
#    variable = bounds_dummy
#    bounded_variable = trapped_2
#    bound_type = lower
#    bound_value = ${bound_value_min}
#  []
#  [trapped_3_upper_bound]
#    type = ConstantBounds
#    variable = bounds_dummy
#    bounded_variable = trapped_3
#    bound_type = upper
#    bound_value = ${bound_value_max}
#  []
#  [trapped_3_lower_bound]
#    type = ConstantBounds
#    variable = bounds_dummy
#    bounded_variable = trapped_3
#    bound_type = lower
#    bound_value = ${bound_value_min}
#  []
#[]

[Variables]
  [concentration]
    order = FIRST
    family = LAGRANGE
    initial_condition = ${initial_concentration}
  []
#  [trapped_1]
#    order = FIRST
#    family = LAGRANGE
#    block = 0
#    outputs = none
#  []
#  [trapped_2]
#    order = FIRST
#    family = LAGRANGE
#    initial_condition = '${fparse initial_concentration_trap_2 * trapping_site_fraction_2 * N}'
#    outputs = none
#  []
  [trapped_3]
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0 #'${fparse initial_concentration_trap_3 * trapping_site_fraction_3 * N}'
    outputs = none
  []
[]

[AuxVariables]
  [temperature]
  []
  [bounds_dummy]
    order = FIRST
    family = LAGRANGE
  []
  [trap_dist_aux]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [temperature_Aux]
    type = FunctionAux
    variable = temperature
    function = Temperature_function
    execute_on = 'initial timestep_end linear'
  []
  [trap_dist_aux_kernel]
    type = FunctionAux
    variable = trap_dist_aux
    function = trap_3_distribution_function
    execute_on = 'initial'
  []
[]

[Kernels]
  [diffusion_implantation]
    type = ADMatDiffusion
    variable = concentration
    diffusivity = Diffusivity
    extra_vector_tags = ref
  []
  [time_diffusion_implantation]
    type = ADTimeDerivative
    variable = concentration
    extra_vector_tags = ref
  []
  [source]
    type = ADBodyForce
    variable = concentration
    function = concentration_source_norm_function
  []

  # trapping kernel
#  [coupled_time_trap_1]
#    type = ADCoefCoupledTimeDerivative
#    variable = concentration
#    v = trapped_1
#    coef = ${trap_per_free_1}
#    block = 0
#    extra_vector_tags = ref
#  []
#  [coupled_time_trap_2_implantation]
#    type = ADCoefCoupledTimeDerivative
#    variable = concentration
#    v = trapped_2
#    coef = ${trap_per_free_2}
#    extra_vector_tags = ref
#  []
  [coupled_time_trap_3_implantation]
    type = ADCoefCoupledTimeDerivative
    variable = concentration
    v = trapped_3
    coef = ${trap_per_free_3}
    extra_vector_tags = ref
  []
[]

[NodalKernels]
  # First traps
#  [time_1]
#    type = TimeDerivativeNodalKernel
#    variable = trapped_1
#  []
#  [trapping_1]
#    type = TrappingNodalKernel
#    variable = trapped_1
#    mobile_concentration = concentration
#    alpha_t = '${trapping_rate_prefactor}'
#    trapping_energy = '${trapping_energy}'
#    N = '${N}'
#    Ct0 = 'trap_1_distribution_function'
#    temperature = 'temperature'
#    trap_per_free = ${trap_per_free_1}
#    extra_vector_tags = ref
#  []
#  [release_1]
#    type = ReleasingNodalKernel
#    variable = trapped_1
#    alpha_r = '${release_rate_profactor}'
#    detrapping_energy = '${detrapping_energy_1}'
#    temperature = 'temperature'
#  []

  # Second traps
#  [time_2]
#    type = TimeDerivativeNodalKernel
#    variable = trapped_2
#  []
#  [trapping_2_implantation]
#    type = TrappingNodalKernel
#    variable = trapped_2
#    mobile_concentration = concentration
#    alpha_t = '${trapping_rate_prefactor}'
#    trapping_energy = '${trapping_energy}'
#    N = '${N}'
#    Ct0 = '${trapping_site_fraction_2}'
#    temperature = 'temperature'
#    trap_per_free = ${trap_per_free_2}
#    extra_vector_tags = ref
#  []
#  [release_2_implantation]
#    type = ReleasingNodalKernel
#    variable = trapped_2
#    alpha_r = '${release_rate_profactor}'
#    detrapping_energy = '${detrapping_energy_2}'
#    temperature = 'temperature'
#  []

  # Third traps
  [time_3]
    type = TimeDerivativeNodalKernel
    variable = trapped_3
  []
  [trapping_3_implantation]
    type = TrappingNodalKernel
    variable = trapped_3
    mobile_concentration = concentration
    alpha_t = '${trapping_rate_prefactor}'
    trapping_energy = '${trapping_energy}'
    N = '${N}'
    Ct0 = 'trap_3_distribution_function'
    temperature = 'temperature'
    trap_per_free = ${trap_per_free_3}
    extra_vector_tags = ref
  []
  [release_3_implantation]
    type = ReleasingNodalKernel
    variable = trapped_3
    alpha_r = '${release_rate_profactor}'
    detrapping_energy = '${detrapping_energy_3}'
    temperature = 'temperature'
  []
[]

[BCs]
  [left_recombination]
    type = BinaryRecombinationBC
    variable = concentration
    v = concentration
    boundary = left
    Kr = recombination_rate
    #Kd = 0   # no dissociation from gas phase on implanted surface
  []
  [right]
    type = ADDirichletBC
    variable = concentration
    boundary = right
    value = 0
  []
[]

[Materials]
  [Diffusivity_implantation]
    type = ADDerivativeParsedMaterial
    property_name = 'Diffusivity'
    functor_names = 'Temperature_function'
    functor_symbols = 'Temperature_function'
    expression = '${diffusivity_coefficient} * exp(- ${E_D} / ${kB} / Temperature_function)'
    block = 0
    output_properties = 'Diffusivity'
  []
  [Diffusivity_Tungsten]
    type = ADDerivativeParsedMaterial
    property_name = 'Diffusivity'
    functor_names = 'Temperature_function'
    functor_symbols = 'Temperature_function'
    expression = '${diffusivity_coefficient} * exp(- ${E_D} / ${kB} / Temperature_function)'
    block = 1
    output_properties = 'Diffusivity'
  []
  [converter_to_regular]
    type = MaterialADConverter
    ad_props_in = 'Diffusivity'
    reg_props_out = 'Diffusivity_nonAD'
    outputs = none
  []
  [recombination_rate]
    type = ADDerivativeParsedMaterial
    property_name = 'recombination_rate'
    functor_names = 'Temperature_function'
    functor_symbols = 'T'
    expression = '${Kr0} * exp(-${Ec} / T)'
    block = 0
  []
[]

[Functions]
  [Temperature_function]
    type = ParsedFunction
    expression = 'if(t < ${TDS_initial_time} - ${temperature_cool_duration},
                      ${temperature_low},
                    if(t < ${TDS_initial_time},
                      ${temperature_low} + (${temperature_TDS_start} - ${temperature_low}) *
                        (t - (${TDS_initial_time} - ${temperature_cool_duration})) / ${temperature_cool_duration},
                    if(t < ${TDS_initial_time} + (${temperature_high} - ${temperature_TDS_start}) / ${temperature_rate},
                      ${temperature_TDS_start} + ${temperature_rate} * (t - ${TDS_initial_time}),
                      ${temperature_high})))'
  []

  [surface_flux_function]
    type = ParsedFunction
    expression = 'if(t < ${TDS_initial_time} - 10.0,
                    ${flux_high},
                    if(t < ${TDS_initial_time},
                    ${flux_high} * (${TDS_initial_time} - t) / 10.0,
                    ${flux_low}))'
  []

  [source_distribution_function]
    type = ParsedFunction
    expression = '1 / ( ${width_source} * sqrt(2 * pi) ) * exp(-0.5 * ((x - ${depth_source}) / ${width_source} ) ^ 2)'
  []

  [trap_3_distribution_function]
    type = ParsedFunction
    expression = '${a3} + ${b3} * exp(-x / ${x0})'
  []

#  [trap_1_distribution_function]
#    type = ParsedFunction
#    expression = ' ${trapping_site_fraction_1} / ( ${width_trap1} * sqrt(2 * pi) ) * exp(-0.5 * ((x - ${depth_source}) / ${width_trap1}) ^ 2)'
#  []

  [concentration_source_norm_function]
    type = ParsedFunction
    symbol_names = 'source_distribution_function surface_flux_function'
    symbol_values = 'source_distribution_function surface_flux_function'
    expression = 'source_distribution_function * surface_flux_function'
  []

  [max_dt_size_function]
    type = ParsedFunction
    expression = 'if(t < ${TDS_initial_time} - 20,  ${step_interval_mid},
                    if(t < ${TDS_initial_time} + 20,   ${step_interval_min},
                                                    ${step_interval_min}))'
  []

  [max_dt_size_function_coarse]
    type = ParsedFunction
    expression = 'if(t<${TDS_initial_time}  , ${step_interval_mid},
                  if(t<${TDS_critial_time_1}, ${step_interval_max},
                  if(t<${TDS_critial_time_2}, ${step_interval_min}, ${step_interval_max})))'
  []
[]

[Postprocessors]
  [flux_surface_left]
    type = SideDiffusiveFluxIntegral
    variable = concentration
    diffusivity = 'Diffusivity_nonAD'
    boundary = 'left'
    outputs = none
  []
  [scaled_flux_surface_left]
    type = ScalePostprocessor
    scaling_factor = '${units 1 m^2 -> mum^2}'
    value = flux_surface_left
    execute_on = 'initial nonlinear linear timestep_end'
    outputs = 'console csv exodus'
  []
  [temperature_sample]
    type = PointValue
    variable = temperature
    point = '0 0 0'
    execute_on = 'initial timestep_end'
    outputs = 'console csv'
  []
  [flux_surface_right]
    type = SideDiffusiveFluxIntegral
    variable = concentration
    diffusivity = 'Diffusivity_nonAD'
    boundary = 'right'
    outputs = none
  []
  [scaled_flux_surface_right]
    type = ScalePostprocessor
    scaling_factor = '${units 1 m^2 -> mum^2}'
    value = flux_surface_right
    execute_on = 'initial nonlinear linear timestep_end'
    outputs = none
  []
  [max_time_step_size]
    type = FunctionValuePostprocessor
    function = max_dt_size_function
    execute_on = 'initial nonlinear linear timestep_end'
    outputs = none
  []
  [total_mobile]
    type = ElementIntegralVariablePostprocessor
    variable = concentration
    execute_on = 'initial timestep_end'
    outputs = 'csv'
  []
  [total_trapped_3]
    type = ElementIntegralVariablePostprocessor
    variable = trapped_3
    execute_on = 'initial timestep_end'
    outputs = 'csv'
  []
  [total_trapped_3_physical]
    type = ScalePostprocessor
    scaling_factor = '${N}'
    value = total_trapped_3
    execute_on = 'initial timestep_end'
    outputs = none
  []
  [total_mobile_per_m2]
    type = ScalePostprocessor
    scaling_factor = '${units 1 m -> mum}'   # converts at/mum to at/m^2
    value = total_mobile
    execute_on = 'initial timestep_end'
    outputs = 'csv console'
  []
  [total_trapped_3_per_m2]
    type = ScalePostprocessor
    scaling_factor = '${units 1 m -> mum}'
    value = total_trapped_3_physical   # <-- changed from total_trapped_3
    execute_on = 'initial timestep_end'
    outputs = 'csv console'
  []
  [total_retained_per_m2]
    type = LinearCombinationPostprocessor
    pp_names  = 'total_mobile_per_m2 total_trapped_3_per_m2'
    pp_coefs  = '1.0                 1.0'
    execute_on = 'initial timestep_end'
    outputs = 'csv console'
  []
  [concentration_at_left]
    type = PointValue
    variable = concentration
    point = '0 0 0'
    execute_on = 'initial timestep_end'
    outputs = none
  []
  [temperature_at_left]
    type = PointValue
    variable = temperature
    point = '0 0 0'
    execute_on = 'initial timestep_end'
    outputs = none
  []
  [recombination_flux_left]
    type = ParsedPostprocessor
    pp_names = 'concentration_at_left temperature_at_left'
    expression = 'concentration_at_left * concentration_at_left * ${Kr0} * exp(-${Ec} / temperature_at_left)'
    execute_on = 'initial timestep_end'
    outputs = 'csv console'
  []
  [scaled_recombination_flux]
    type = ScalePostprocessor
    scaling_factor = '${units 1 m^2 -> mum^2}'
    value = recombination_flux_left
    execute_on = 'initial timestep_end'
    outputs = 'csv console'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  petsc_options_iname = '-pc_type'
  petsc_options_value  = 'lu'
  solve_type = NEWTON

  end_time = ${simulation_time}
  automatic_scaling = true
  nl_rel_tol = 1e-4
  nl_max_its = 100
  nl_abs_tol = 1e-6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1.0
    iteration_window = 5
    optimal_iterations = 26
    growth_factor = 1.1
    cutback_factor = 0.9
    cutback_factor_at_failure = 0.9
    timestep_limiting_postprocessor = max_time_step_size
  []
[]

[Outputs]
  file_base = 'pisces_val_600K_trap3_out'
  [csv]
    type = CSV
    start_time = ${outputs_initial_time}
  []
  [exodus]
    type = Exodus
    start_time = ${outputs_initial_time}
    output_material_properties = true
    time_step_interval = 20
  []
[]
