function find_voltage_source_branch_bus(math)
    for (b, branch) in math["branch"]
        if branch["source_id"] == "voltage_source.source"
            return b, branch["f_bus"]
        end
    end
    return error()
end


