
function run_mc_opf(math, formulation, post, optimizer)
    pm = PMD.instantiate_mc_model(
        math,
        formulation,
        post;
        # ref_extensions=ref_extensions,
        multinetwork=false,
        # kwargs...
    )
    result = PMD.solve_mc_opf(math, formulation, optimizer)
    return result
end


function build_mc_opf_alt_4w(pm::PMD.AbstractExplicitNeutralIVRModel)
    # Variables
    PMD.variable_mc_bus_voltage(pm)
    PMD.variable_mc_branch_current(pm)
    PMD.variable_mc_load_current(pm)
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
        # PMD.constraint_mc_thermal_limit_from(pm, i)
        # PMD.constraint_mc_thermal_limit_to(pm, i)
    end

    # for i in PMD.ids(pm, :switch)
    #     PMD.constraint_mc_switch_current(pm, i)
    #     PMD.constraint_mc_switch_state(pm, i)

    #     PMD.constraint_mc_switch_current_limit(pm, i)
    #     PMD.constraint_mc_switch_thermal_limit(pm, i)
    # end

    for i in PMD.ids(pm, :bus)
        PMD.constraint_mc_current_balance(pm, i)
    end

    # Objective
    PMD.objective_mc_min_fuel_cost(pm)
end


