## Inputs from RFPIE_Exposure2.i
total_H = 1613172964.681    # total_H_mobile from RFPIE_Exposure2.i
initial_concentration_trap_1 = 4.4e-10 # (-)
initial_concentration_trap_2 = 4.4e-10 # (-)
initial_concentration_trap_3 = 1.4e-10 # (-)
########################################################################

# General parameters
kB = '${units 1.380649e-23 J/K}'  # Boltzmann constant (from PhysicalConstants.h - https://physics.nist.gov/cgi-bin/cuu/Value?r)

# Model parameters
simulation_time = '${units 1e4 s}'  # 10,000 s
outputs_initial_time = '${units 0 s}'
bound_value_min = '${units -1e-10 at/mum^3}'

nx_scale = 5

# Diffusion parameters
diffusivity_coefficient = '${units 4.1e-7 m^2/s -> mum^2/s}'    #max diffusivity coeff when x < 15e-9 m
E_D = '${units 0.39 eV -> J}'           # activity energy for diffusion
mesh_volume = '${units 1.4107e-8 m^3 -> mum^3}'
initial_concentration = '${fparse ${total_H} / ${mesh_volume}}'

# Traps parameters
N = '${units 6.25e28 at/m^3 -> at/mum^3}'
trapping_energy = '${fparse ${units 0.39 eV -> J} / ${kB}}'
detrapping_energy_1 = '${fparse ${units 1.2 eV -> J} / ${kB}}'  # Should release ~300-500 C
detrapping_energy_2 = '${fparse ${units 1.6 eV -> J} / ${kB}}'  # Releases at higher temperatures
detrapping_energy_3 = '${fparse ${units 3.1 eV -> J} / ${kB}}'  # Very deep traps, release ~700-800 C
trapping_site_fraction_1 = 0.002156 # (-)
trapping_site_fraction_2 = 0.00175 # (-)
trapping_site_fraction_3 = 0.0020 # (-)
trapping_rate_prefactor = '${units 9.1316e12 1/s}'
release_rate_profactor = '${units 8.4e12 1/s}'
trap_per_free_1 = 1e6 # (-)
trap_per_free_2 = 1e4 # (-)
trap_per_free_3 = 1e4 # (-)
#width_trap1 = '${units 10e-9 m -> mum}'

# Thermal parameters
#temperature_low = '${units 300 K}'  # room temp
#temperature_high = '${units 1173 K}'    # 900 C
#temperature_rate = '${units ${fparse 15 / 60} K/s}'     # 15 C/min

# Mesh copy-pasted from RFPIE_Exposure.i
[Mesh]
    [cartesian]
      type = CartesianMeshGenerator
      dim = 1
      dx = '${fparse 5 * ${units 4e-9 m -> mum}}  ${units 1e-8 m -> mum}  ${units 1e-7 m -> mum}
           ${units 1e-6 m -> mum}                ${units 1e-5 m -> mum}  ${fparse 10 * ${units 4.88e-5 m -> mum}}'
      ix = '${fparse 5 * ${nx_scale}}             ${nx_scale}             ${nx_scale}
           ${nx_scale}                           ${nx_scale}             ${fparse 10 * ${nx_scale}}'
    []
[]

[Problem]
    type = ReferenceResidualProblem
    extra_tag_vectors = 'ref'
    reference_vector = 'ref'
[]
  
[Bounds]
    [concentratio_lower]
      type = ConstantBounds
      bounded_variable = concentration
      bound_type = lower
      variable = bounds_dummy
      bound_value = 0.0
    []
    [trapped_1_lower]
      type = ConstantBounds
      bounded_variable = trapped_1
      bound_type = lower
      variable = bounds_dummy
      bound_value = ${bound_value_min}
    []
    [trapped_2_lower]
      type = ConstantBounds
      bounded_variable = trapped_2
      bound_type = lower
      variable = bounds_dummy
      bound_value = ${bound_value_min}
    []
    [trapped_3_lower]
      type = ConstantBounds
      bounded_variable = trapped_3
      bound_type = lower
      variable = bounds_dummy
      bound_value = ${bound_value_min}
    []
    [trapped_2_upper]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = trapped_2
      bound_type = upper
      bound_value = '${fparse ${N} * ${trapping_site_fraction_2}}'
    []
    [trapped_3_upper]
      type = ConstantBounds
      variable = bounds_dummy
      bounded_variable = trapped_3
      bound_type = upper
      bound_value = '${fparse ${N} * ${trapping_site_fraction_3}}'
    []
[]
  
[Variables]
    [concentration]
      order = FIRST
      family = LAGRANGE
      initial_condition = ${initial_concentration}
    []
    [trapped_1]
      order = FIRST
      family = LAGRANGE
      initial_condition = '${fparse initial_concentration_trap_1 * trapping_site_fraction_1 * N}'
      outputs = none
    []
    [trapped_2]
      order = FIRST
      family = LAGRANGE
      initial_condition = '${fparse initial_concentration_trap_2 * trapping_site_fraction_2 * N}'
      outputs = none
    []
    [trapped_3]
      order = FIRST
      family = LAGRANGE
      initial_condition = '${fparse initial_concentration_trap_3 * trapping_site_fraction_3 * N}'
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
[]
  
[AuxKernels]
    [temperature_aux]
      type = FunctionAux
      variable = temperature
      function = Temperature_function
      execute_on = 'initial timestep_begin'
    []
[]
  
[Kernels]
    [diffusion]
      type = ADMatDiffusion
      variable = concentration
      diffusivity = Diffusivity
      extra_vector_tags = ref
    []
    [time_derivative]
      type = ADTimeDerivative
      variable = concentration
      extra_vector_tags = ref
    []
    [coupled_trap_1]
      type = ADCoefCoupledTimeDerivative
      variable = concentration
      v = trapped_1
      coef = ${trap_per_free_1}
      extra_vector_tags = ref
    []
    [coupled_trap_2]
      type = ADCoefCoupledTimeDerivative
      variable = concentration
      v = trapped_2
      coef = ${trap_per_free_2}
      extra_vector_tags = ref
    []
    [coupled_trap_3]
      type = ADCoefCoupledTimeDerivative
      variable = concentration
      v = trapped_3
      coef = ${trap_per_free_3}
      extra_vector_tags = ref
    []
[]
  
[NodalKernels]
    [time_1]
      type = TimeDerivativeNodalKernel
      variable = trapped_1
    []
#    [trapping_1]
#      type = TrappingNodalKernel
#      variable = trapped_1
#      mobile_concentration = concentration
#      alpha_t = '${trapping_rate_prefactor}'
#      trapping_energy = '${trapping_energy}'
#      N = '${N}'
#      Ct0 = 'trap_1_distribution_function'
#      temperature = 'temperature'
#      trap_per_free = ${trap_per_free_1}
#      extra_vector_tags = ref
#    []
    [release_1]
      type = ReleasingNodalKernel
      variable = trapped_1
      alpha_r = '${release_rate_profactor}'
      detrapping_energy = '${detrapping_energy_1}'
      temperature = 'temperature'
    []
    [time_2]
      type = TimeDerivativeNodalKernel
      variable = trapped_2
    []
#    [trapping_2]
#      type = TrappingNodalKernel
#      variable = trapped_2
#      mobile_concentration = concentration
#      alpha_t = '${trapping_rate_prefactor}'
#      trapping_energy = '${trapping_energy}'
#      N = '${N}'
#      Ct0 = '${trapping_site_fraction_2}'
#      temperature = 'temperature'
#      trap_per_free = ${trap_per_free_2}
#      extra_vector_tags = ref
#    []
    [release_2]
      type = ReleasingNodalKernel
      variable = trapped_2
      alpha_r = '${release_rate_profactor}'
      detrapping_energy = '${detrapping_energy_2}'
      temperature = 'temperature'
    []
    [time_3]
      type = TimeDerivativeNodalKernel
      variable = trapped_3
    []
#    [trapping_3]
#      type = TrappingNodalKernel
#      variable = trapped_3
#      mobile_concentration = concentration
#      alpha_t = '${trapping_rate_prefactor}'
#      trapping_energy = '${trapping_energy}'
#      N = '${N}'
#      Ct0 = '${trapping_site_fraction_3}'
#      temperature = 'temperature'
#      trap_per_free = ${trap_per_free_3}
#      extra_vector_tags = ref
#    []
    [release_3]
      type = ReleasingNodalKernel
      variable = trapped_3
      alpha_r = '${release_rate_profactor}'
      detrapping_energy = '${detrapping_energy_3}'
      temperature = 'temperature'
    []
[]
  
[BCs]
    [left]
      type = ADDirichletBC  # Acts as a sink for outgassing
      variable = concentration
      boundary = left
      value = 0
    []
#    [right]    # commented because only 1 face outgases in Sandia set up
#      type = ADDirichletBC
#      variable = concentration
#      boundary = right
#      value = 0
#    []
[]
  
[Materials]
    [Diffusivity]
      type = ADDerivativeParsedMaterial
      property_name = 'Diffusivity'
      functor_names = 'Temperature_function'
      functor_symbols = 'Temperature_function'
      expression = '${diffusivity_coefficient} * exp(- ${E_D} / ${kB} / Temperature_function)'
      output_properties = 'Diffusivity'
    []
    [converter]
      type = MaterialADConverter
      ad_props_in = 'Diffusivity'
      reg_props_out = 'Diffusivity_nonAD'
      outputs = none
    []
[]
  
[Functions]
    [Temperature_function]
        type = PiecewiseLinear
        data_file = 'TDSTempRamp.csv'
        format = 'columns'
        x_title = 'time [s]'
        y_title = 'T [K]'
    []
  
    [trap_1_distribution_function]
      type = ParsedFunction
      expression = '${trapping_site_fraction_1}'  # Sets traps as a uniform fraction of N everywhere
    []

    [dt_cap_function]
      type = ParsedFunction
      expression = '20.0'  # in seconds, or whatever max timestep you want
    []

[]
  

[Postprocessors]
    [flux_surface_left]     # Gives desorbing flux (atoms/s) at surface
        type = SideDiffusiveFluxIntegral
        variable = concentration
        diffusivity = 'Diffusivity_nonAD'
        boundary = 'left'
        outputs = none
    []
    [flux_surface_left_area_norm]   # Flux normalized to atoms/m^2/s
        type = ScalePostprocessor
        scaling_factor = ${fparse 1.0/2.827e-5}   # 1/area
        value = flux_surface_left
        execute_on = 'initial nonlinear linear timestep_end'
        outputs = 'console csv_out exodus'
    []
    [scaled_flux_surface_left]
        type = ScalePostprocessor
        scaling_factor = '${units 1 m^2 -> mum^2}'
        value = flux_surface_left
        execute_on = 'initial nonlinear linear timestep_end'
        outputs = 'console csv_out exodus'
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
    [flux_surface_right_area_norm]
        type = ScalePostprocessor
        scaling_factor = ${fparse 1.0/2.827e-5}   # 1/area
        value = flux_surface_right
        execute_on = 'initial nonlinear linear timestep_end'
        outputs = 'console csv_out exodus'
    []
    [avg_temperature]     # Gives temp vs time
        type = ElementAverageValue
        variable = temperature
        outputs = 'console csv_out'
        execute_on = 'initial timestep_end'
    []
    [total_H]     # Multiple by mesh volume to get total H amount
        type = ElementIntegralVariablePostprocessor
        variable = concentration
        outputs = 'console csv_out'
        execute_on = 'initial timestep_end'
    []
    [d2_flux_surface_left]
        type = ScalePostprocessor
        scaling_factor = 0.5
        value = flux_surface_left_area_norm
        execute_on = 'initial nonlinear linear timestep_end'
        outputs = 'console csv_out exodus'
    []

    [dt_cap]
        type = FunctionValuePostprocessor
        function = dt_cap_function
        execute_on = 'initial nonlinear linear timestep_end'
    []

    # Surface flux at left surface (node 0)
    [surface_flux_left]
        type = SideDiffusiveFluxIntegral
        variable = concentration
        #node = 0
        diffusivity = 'Diffusivity_nonAD'
        boundary = 'left'
        outputs = 'console csv_out'
    []

    # Trapped populations at left surface (node 0)
    [trap1_left]
        type = ElementAverageValue
        variable = trapped_1
        #node = 0
        outputs = 'console csv_out'
    []

    [trap2_left]
        type = ElementAverageValue
        variable = trapped_2
        #node = 0
        outputs = 'console csv_out'
    []

    [trap3_left]
        type = ElementAverageValue
        variable = trapped_3
        #node = 0
        outputs = 'console csv_out'
    []

    # Bulk concentration (volume average)
    [bulk_concentration]
        type = ElementAverageValue
        variable = concentration
        outputs = 'console csv_out'
    []

    # Bulk trap populations (volume average)
    [bulk_trap1]
        type = ElementAverageValue
        variable = trapped_1
        outputs = 'console csv_out'
    []

    [bulk_trap2]
        type = ElementAverageValue
        variable = trapped_2
        outputs = 'console csv_out'
    []

    [bulk_trap3]
        type = ElementAverageValue
        variable = trapped_3
        outputs = 'console csv_out'
    []

    # Temperature at left surface (node 0)
    [temperature_left]
        type = ElementAverageValue
        variable = temperature
        #node = 0
        outputs = 'console csv_out'
    []

    # Bulk average temperature
    [avg_temperature_bulk]
        type = ElementAverageValue
        variable = temperature
        outputs = 'console csv_out'
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
    solve_type = NEWTON

    #petsc_options_iname = '-pc_type -snes_type'
    #petsc_options_value = 'lu vinewtonrsls'

    # Use VI-Newton (keeps bounds) but with backtracking
    petsc_options_iname  = '-snes_type -snes_linesearch_type -ksp_type -pc_type'
    petsc_options_value  = 'vinewtonrsls bt               gmres     lu'

    end_time = ${simulation_time}
    line_search = 'bt'  #backtracking
    automatic_scaling = true
    nl_rel_tol = 1e-6
    nl_abs_tol = 1e-12
    nl_max_its = 50
    [TimeStepper]
        type = IterationAdaptiveDT
        dt = 1e-10
        iteration_window = 5
        optimal_iterations = 12
        growth_factor = 1.25
        cutback_factor = 0.5
        cutback_factor_at_failure = 0.2
        timestep_limiting_postprocessor = dt_cap
    []
[]

[Outputs]
    file_base = 'TDS_out'
    #output_variables = 'temperature'
    [csv_out]
        type = CSV
        execute_postprocessors_on = 'nonlinear timestep_end'
        time_step_interval = 1
    []
    [exodus]
        type = Exodus
        start_time = ${outputs_initial_time}
        output_material_properties = true
        time_step_interval = 1
    []
[]