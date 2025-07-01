"""
	function build_mc_opf(
		pm::AbstractUnbalancedPowerModel
	)

Constructor for Optimal Power Flow
"""
function build_mc_opf(pm::PMD.AbstractUnbalancedPowerModel)
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_power(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_power(pm)
    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_storage_power(pm)

    PMD.constraint_mc_model_voltage(pm)

    for i in PMD.ids(pm, :ref_buses)
        # PMD.constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_reference(pm, i)
    end

    # generators should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_power_balance(pm, i)
        constraint_mc_bus_voltage_balance(pm, i)

        if i ∉ PMD.ids(pm, :ref_buses)
            constraint_mc_bus_voltage_angle_difference(pm, i)
        end
    end

    for i in PMD.ids(pm, :storage)
        PMD.constraint_storage_state(pm, i)
        PMD.constraint_storage_complementarity_nl(pm, i)
        PMD.constraint_mc_storage_losses(pm, i)
        PMD.constraint_mc_storage_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_ohms_yt_from(pm, i)
        PMD.constraint_mc_ohms_yt_to(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
        PMD.constraint_mc_ampacity_from(pm, i)
        PMD.constraint_mc_ampacity_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
        PMD.constraint_mc_switch_ampacity(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

    PMD.objective_mc_min_fuel_cost(pm)
end


"""
	function build_mc_opf(
		pm::AbstractUnbalancedIVRModel
	)

constructor for OPF in current-voltage variable space
"""
function build_mc_opf(pm::PMD.AbstractUnbalancedIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_switch_current(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_load_current(pm)

    # Constraints
    for i in PMD.ids(pm, :ref_buses)
        # constraint_mc_voltage_reference(pm, i)
        # PMD.constraint_mc_theta_ref(pm, i)
        constraint_mc_voltage_reference(pm, i)
    end

    # gens should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
    end

    # loads should be constrained before KCL, or Pd/Qd undefined
    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
        constraint_mc_bus_voltage_balance(pm, i)  # This is added to the IVR

        if i ∉ PMD.ids(pm, :ref_buses)
            constraint_mc_bus_voltage_angle_difference(pm, i)
        end
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)

        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_voltage_angle_difference(pm, i)

        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
        constraint_mc_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_state(pm, i)
        PMD.constraint_mc_switch_current_limit(pm, i)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_power(pm, i)
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end


"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralIVRModel
	)

constructor for OPF in current-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::PMD.AbstractExplicitNeutralIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.ariable_mc_load_current(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_current(pm)
    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_current(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_current(pm)

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
            # constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        PMD.constraint_mc_generator_current(pm, id)
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
        PMD.constraint_mc_load_current(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_current(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_current_from(pm, i)
        PMD.constraint_mc_current_to(pm, i)
        PMD.constraint_mc_bus_voltage_drop(pm, i)

        PMD.constraint_mc_branch_current_limit(pm, i)
        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_current(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
        constraint_mc_bus_voltage_balance(pm, i)  # This is added to the en_IVR
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end



"""
	function build_mc_opf(
		pm::AbstractExplicitNeutralACRModel
	)

constructor for OPF in power-voltage variable space with explicit neutrals
"""
function build_mc_opf(pm::PMD.AbstractExplicitNeutralACRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_power(pm)
    PMD.variable_mc_load_power(pm)
    PMD.variable_mc_generator_power(pm)
    PMD.variable_mc_transformer_power(pm)
    PMD.variable_mc_switch_power(pm)

    # Constraints
    for i in PMD.ids(pm, :bus)

        if i in PMD.ids(pm, :ref_buses)
            PMD.constraint_mc_voltage_reference(pm, i)
        end

        PMD.constraint_mc_voltage_absolute(pm, i)
        PMD.constraint_mc_voltage_pairwise(pm, i)
    end

    # components should be constrained before KCL, or the bus current variables might be undefined

    for id in PMD.ids(pm, :gen)
        PMD.constraint_mc_generator_power(pm, id)
        # PMD.constraint_mc_generator_current(pm, id)
    end

    for id in PMD.ids(pm, :load)
        PMD.constraint_mc_load_power(pm, id)
    end

    for i in PMD.ids(pm, :transformer)
        PMD.constraint_mc_transformer_voltage(pm, i)
        PMD.constraint_mc_transformer_power(pm, i)

        PMD.constraint_mc_transformer_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :branch)
        PMD.constraint_mc_ohms_yt_from(pm, i)
        PMD.constraint_mc_ohms_yt_to(pm, i)

        PMD.constraint_mc_branch_current_limit(pm, i)
        PMD.constraint_mc_thermal_limit_from(pm, i)
        PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    for i in PMD.ids(pm, :switch)
        PMD.constraint_mc_switch_power(pm, i)
        PMD.constraint_mc_switch_state(pm, i)

        PMD.constraint_mc_switch_current_limit(pm, i)
        PMD.constraint_mc_switch_thermal_limit(pm, i)
    end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_power_balance(pm, i)
        constraint_mc_bus_voltage_balance(pm, i)  # This is added to the ACR
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end
