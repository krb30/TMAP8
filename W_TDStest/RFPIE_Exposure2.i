on_time = 4 # 4s in the middle of a 6s shot
off_time = 17   
number_of_shots = 10    # 250 shots eventually
cycle_time = '${fparse on_time + off_time}'
simulation_time = '${fparse number_of_shots * cycle_time}'
  
nx_scale = 5
diffusivity_D = '${units 3e-10 m^2/s -> mum^2/s}'
recombination_parameter_enclos2 = '${units 2e-31 m^4/at/s -> mum^4/at/s}'
recombination_coefficient_parameter_enclos1_TMAP4 = '${units 1e-27 m^4/at/s -> mum^4/at/s}' # Specify no/perfect recom at downstream side
width = '${units 2.4e-9 m -> mum}'  # omega - characteristic width of the normal distribution
depth = '${units 14e-9 m -> mum}'   # mu - depth of the normal distribution from the upstream side

flux_base = '${units 1.3e21 at/m^2/s -> at/mum^2/s}'
flux_high = '${fparse 3 * flux_base}'   # 3*flux for He plasma

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
  [filled_trap1]
  []
  [filled_trap2]
  []
  [filled_trap3]
  []
  [temperature]
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
  [temperature_aux]
    type = ConstantAux
    variable = temperature
    value = 300   # K - Samples stay at room temperature    ########### This could be what I'd change if I want to look at temperature dependence since plasma inc temp ~50 C
  []
[]

[NodalKernels]
  [trap1_trapping]
    type = TrappingNodalKernel
    variable = concentration
    filled_trap = filled_trap1
    Ct0 = 1.0
    N = 6.3e28  # number density of W
    alpha_t = 1e6
    mobile_concentration = concentration
    temperature = 'temperature'
  []

  [trap2_trapping]
    type = TrappingNodalKernel
    variable = concentration
    filled_trap = filled_trap2
    Ct0 = 0.1
    N = 6.3e28
    alpha_t = 5e5
    mobile_concentration = concentration
    temperature = 'temperature'
  []

  [trap3_trapping]
    type = TrappingNodalKernel
    variable = concentration
    filled_trap = filled_trap3
    Ct0 = 0.01
    N = 6.3e28
    alpha_t = 1e5
    mobile_concentration = concentration
    temperature = 'temperature'
  []
[]


[BCs]
  [left]
    type = ADMatNeumannBC
    variable = concentration
    boundary = left
    value = 0
    boundary_material = flux_on_left  # Calculated in the [Materials] block
  []
  [right]
    type = ADMatNeumannBC
    variable = concentration
    boundary = right
    value = 0
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
    on_time = 4
    cycle_time = 21
    flux_high = ${flux_high}
    expression = 'if(((t / ${cycle_time}) - floor(t / ${cycle_time})) < ${on_time} / ${cycle_time}, ${flux_high}, 0)'
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
    outputs = 'console csv_out exodus'
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
    outputs = 'console csv_out exodus'
  []
  [total_H_mobile]   # This is the sum of mobile H in the lattice. Once a H atom is trapped, it no longer counts in concentration
    type = ElementIntegralVariablePostprocessor
    variable = concentration
    execute_on = 'timestep_end'
  []
  [filled_trap1_pp]
    type = ElementIntegralVariablePostprocessor
    variable = filled_trap1
    execute_on = 'timestep_end'
    outputs = 'console csv_out'
  []

  [filled_trap2_pp]  
    type = ElementIntegralVariablePostprocessor
    variable = filled_trap2
    execute_on = 'timestep_end'
    outputs = 'console csv_out'
  []

  [filled_trap3_pp]
    type = ElementIntegralVariablePostprocessor
    variable = filled_trap3
    execute_on = 'timestep_end'
    outputs = 'console csv_out'
  []

  [total_H_implanted]
    type = ElementIntegralVariablePostprocessor
    variable = concentration
    execute_on = 'timestep_end'
    # Include other trapped variables in "additional_variables"
    additional_variables = 'trap1 trap2 trap3'
    scale_factor = 1.0
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
[]

[Outputs]
  file_base = 'RFPIE_Exposure2_out'

  [exodus]
    type = Exodus
    output_material_properties = true
    time_step_interval = 2
  []

  [csv_out]
    type = CSV
    execute_postprocessors_on = 'timestep_end'
    time_step_interval = 1
  []

  [chk]
    type = Checkpoint
    num_files = 1
    execute_on = 'timestep_end'
  []
[]
