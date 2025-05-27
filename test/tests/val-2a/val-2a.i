nx_scale = 5
high_dt_max = 300
low_dt_max = 4
simulation_time = '${units 2e4 s}'
diffusivity_D = '${units 3e-10 m^2/s -> mum^2/s}'
recombination_parameter_enclos2 = '${units 2e-31 m^4/at/s -> mum^4/at/s}'
flux_high = '${units 4.9e19 at/m^2/s -> at/mum^2/s}'
flux_low =  '${units 0      at/mum^2/s}'
recombination_coefficient_parameter_enclos1_TMAP4 = '${units 1e-27 m^4/at/s -> mum^4/at/s}' # Specify no/perfect recom at downstream side
width = '${units 2.4e-9 m -> mum}'
depth = '${units 14e-9 m -> mum}'
time_1 = '${units 5820 s}'
time_2 = '${units 9056 s}'
time_3 = '${units 12062 s}'
time_4 = '${units 14572 s}'
time_5 = '${units 17678 s}'

[Variables]
  [concentration]  # (atoms/mum^3/s)
    order = FIRST
    family = LAGRANGE
  []
[]

[Mesh]    # Creates a very fine mesh near the surface (0.8 nm elements) then it gets coarser into the bulk
  [cartesian]
    type = CartesianMeshGenerator
    dim = 1
    dx = '${fparse 5 * ${units 4e-9 m -> mum}}  ${units 1e-8 m -> mum}  ${units 1e-7 m -> mum}
          ${units 1e-6 m -> mum}                ${units 1e-5 m -> mum}  ${fparse 10 * ${units 4.88e-5 m -> mum}}'
    ix = '${fparse 5 * ${nx_scale}}             ${nx_scale}             ${nx_scale}
          ${nx_scale}                           ${nx_scale}             ${fparse 10 * ${nx_scale}}'
  []
[]

[Kernels]
  [diffusion]   # Uses Fick's Law: ∂C/∂t = ∇·(D∇C)
    type = ADMatDiffusion
    variable = concentration
    diffusivity = ${diffusivity_D}
  []
  [time_diffusion]  # Adds the time derivative
    type = ADTimeDerivative
    variable = concentration
  []
  [source]  # Adds volumetric source term based on space and time based on concentration_source_norm_func function below
    type = ADBodyForce
    variable = concentration
    function = concentration_source_norm_func
  []
[]

[AuxVariables]
  [concentration_source]
  []
  [recombination_TMAP4]
  []
[]

[AuxKernels]
  [concentration_source_aux]  
    type = FunctionAux
    variable = concentration_source
    function = concentration_source_norm_func
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [recombination_aux_TMAP4]   # outputs recombination coefficient value over time
    type = FunctionAux
    variable = recombination_TMAP4
    function = '${recombination_coefficient_parameter_enclos1_TMAP4}'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]

[BCs]
  [left]
    type = ADMatNeumannBC
    variable = concentration
    boundary = left
    value = 1
    boundary_material = flux_on_left  # Calculated in the [Materials] block
  []
  [right]
    type = ADMatNeumannBC
    variable = concentration
    boundary = right
    value = 1
    boundary_material = flux_on_right # Calculated in the [Materials] block
  []
[]

[Materials]
  [flux_on_left]
    type = ADDerivativeParsedMaterial
    coupled_variables = 'concentration'
    property_name = 'flux_on_left'
    functor_names = 'Kr_left_func'
    functor_symbols = 'Kr_left_func'
    expression = '- 2 * Kr_left_func * concentration ^ 2'
  []
  [flux_on_right]
    type = ADDerivativeParsedMaterial
    coupled_variables = 'concentration'
    property_name = 'flux_on_right'
    expression = '- 2 * ${recombination_parameter_enclos2} * concentration ^ 2'
  []
[]

[Functions]
  [Kr_left_func] # Recombination coefficient on left boundary w/ units [microns^4/at/s]
    type = ParsedFunction
    expression = '${recombination_coefficient_parameter_enclos1_TMAP4} * (1 - 0.9999 * exp(-6e-5 * t))'
  []

  [surface_flux_func] # Describes the time varying particle flux at the surface w/ units [atoms/mum^2/s]  
    type = ParsedFunction
    expression = 'if(t < ${time_1}, ${flux_high},
                  if(t < ${time_2}, ${flux_low},
                  if(t < ${time_3},  ${flux_high},
                  if(t < ${time_4},  ${flux_low},
                  if(t < ${time_5},  ${flux_high}, ${flux_low}))))) * 0.75'   # Multiplied by 0.75 b/c TRIM calculation showed only 75% of incident flux remained in sample
  []

  [source_distribution] # Likely trapping site density
    type = ParsedFunction
    expression = '1.5 / (${width} * sqrt(2 * pi)) * exp(-0.5 * ((x - ${depth}) / ${width})^2)'
  []

  [concentration_source_norm_func] # Multiplies spatial source dist w/ time-varying flux, units [atoms/microns^2/s]
    type = ParsedFunction
    symbol_names = 'source_distribution surface_flux_func'
    symbol_values = 'source_distribution surface_flux_func'
    expression = 'source_distribution * surface_flux_func'
  []

  [max_dt_size_func] # Dynamically changes the max time step size to ensure numerical stability during rapid changes in BCs, units [s]
    type = ParsedFunction
    expression = 'if(t<${time_1}-100,  ${high_dt_max},
                  if(t<${time_1}+100,  ${low_dt_max},
                  if(t<${time_2}-100,  ${high_dt_max},
                  if(t<${time_2}+100,  ${low_dt_max},
                  if(t<${time_3}-100,  ${high_dt_max},
                  if(t<${time_3}+100,  ${low_dt_max},
                  if(t<${time_4}-100,  ${high_dt_max},
                  if(t<${time_4}+100,  ${low_dt_max},
                  if(t<${time_5}-100,  ${high_dt_max},
                  if(t<${time_5}+100,  ${low_dt_max}, ${high_dt_max}))))))))))'
  []
[]

[Postprocessors]
  [dcdx_left]
    type = ADSideAverageMaterialProperty
    boundary = left
    property = flux_on_left
    outputs = none
  []
  [scaled_recombination_flux_left]
    type = ScalePostprocessor
    scaling_factor = '${fparse -1 * ${units 1 m^2 -> mum^2}}'
    value = dcdx_left
    execute_on = 'initial nonlinear linear timestep_end'
    outputs = 'console csv exodus'
  []
  [dcdx_right]
    type = ADSideAverageMaterialProperty
    boundary = right
    property = flux_on_right
    outputs = none
  []
  [scaled_recombination_flux_right]
    type = ScalePostprocessor
    scaling_factor = '${fparse -1 * ${units 1 m^2 -> mum^2}}'
    value = dcdx_right
    execute_on = 'initial nonlinear linear timestep_end'
    outputs = 'console csv exodus'
  []
  [max_time_step_size]
    type = FunctionValuePostprocessor
    function = max_dt_size_func
    execute_on = 'initial nonlinear linear timestep_end'
    outputs = none
  []
[]

[Preconditioning]
  [SMP]  # Summetric Multi-Processing - used to optimize performance of solver
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2   # backward differentiation formula, order 2
  solve_type = NEWTON 
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  end_time = ${simulation_time}
  automatic_scaling = true
  nl_rel_tol = 5e-7
  [TimeStepper]
    type = IterationAdaptiveDT  # Chooses next time step size based on # of Newton iterations needed for last step
    dt = 1  # initial guess for time step
    optimal_iterations = 6  # if Newton ~6 iterations the dt is "just right"
    growth_factor = 1.1   # if Newton converges faster, inc dt by 10%
    cutback_factor_at_failure = 0.9   # if Newton fails to converge, dec dt by 10% and retry
    timestep_limiting_postprocessor = max_time_step_size
  []
[]

[Outputs]
  file_base = 'val-2a_out'
  csv = true
  [exodus]
    type = Exodus
    output_material_properties = true
    time_step_interval = 2
  []
[]
