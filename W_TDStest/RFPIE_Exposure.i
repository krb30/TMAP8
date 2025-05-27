on_time = 4 # 4s in the middle of a 6s shot
off_time = 17   
number_of_shots = 10    # 250 shots eventually
cycle_time = '${fparse on_time + off_time}'
simulation_time = '${fparse number_of_shots * cycle_time}'
  
nx_scale = 5
diffusivity_D = '${units 3e-10 m^2/s -> mum^2/s}'
recombination_parameter_enclos2 = '${units 2e-31 m^4/at/s -> mum^4/at/s}'
recombination_coefficient_parameter_enclos1_TMAP4 = '${units 1e-27 m^4/at/s -> mum^4/at/s}' # Specify no/perfect recom at downstream side
width = '${units 2.4e-9 m -> mum}'
depth = '${units 14e-9 m -> mum}'

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
    [loop_source]  
        type     = ADBodyForce
        variable = concentration
        function = loop_flux
    []
[]

[AuxVariables]
    [concentration_source]
    []
    [recombination_TMAP4]
    []
    [loop_flux_field]
        family = MONOMIAL                    ## Constant within each element
        order  = CONSTANT
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
    [flux_kernel]
        type = FunctionAux
        variable = loop_flux_field
        function = loop_flux
    []
[]

[BCs]
    [left]
      type              = ADMatNeumannBC
      variable          = concentration
      boundary          = left
      value             = 1
      boundary_material = flux_on_left
    []
    [right]
      type              = ADMatNeumannBC
      variable          = concentration
      boundary          = right
      value             = 1
      boundary_material = flux_on_right
    []
  []
  

[Materials]
    [flux_on_left]
        type = ADDerivativeParsedMaterial
        coupled_variables = 'concentration'
        property_name = 'flux_on_left'
        #functor_names = 'Kr_left_func'
        #functor_symbols = 'Kr_left_func'
        #expression = '- 2 * Kr_left_func * concentration ^ 2'
        expression = '- 2 * ${recombination_parameter_enclos2} * concentration ^ 2'     # Assuming constant recombination across surface, left and right
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

    [loop_flux]
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
        symbol_names = 'source_distribution loop_flux'
        symbol_values = 'source_distribution loop_flux'
        expression = 'source_distribution * loop_flux'
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
    [loop_flux_pp]
        type = FunctionValuePostprocessor
        function = loop_flux
        outputs = 'csv'
        execute_on = 'initial timestep_end'
    []
    [flux_pp]
        type = FunctionValuePostprocessor
        function = loop_flux
        execute_on = 'timestep_end'
        outputs = 'csv'
    [] 
    [concentration_pp]
        type = FunctionValuePostprocessor
        function = concentration_source_norm_func
        execute_on = 'initial nonlinear timestep_end'
        outputs = 'csv'
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
    scheme = bdf2
    dt = 0.1
    start_time = 0.0
    end_time = ${simulation_time}
[]  

[Outputs]
    csv = false
    [csv_flux_only]
        type = CSV
        file_base = 'rfpie_flux'
        time_column = true
        postprocessors = 'flux_pp'
        execute_on = 'initial timestep_end'
    []

    [csv_full_data]
        type = CSV
        file_base = 'rfpie_full'
        time_column = true
        postprocessors = 'flux_pp scaled_recombination_flux_left scaled_recombination_flux_right concentration_pp'
        execute_on = 'initial timestep_end'
    []

    [exodus]
        type = Exodus
        output_material_properties = true
        time_step_interval = 2
    []
[]

