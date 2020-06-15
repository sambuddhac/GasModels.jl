import XMLDict
import ZipFile

function _read_gaslib_file(gaslib_zip, extension::String)
    index = findfirst(x -> occursin(extension, x.name), gaslib_zip.files)
    subdata_xml = ZipFile.read(gaslib_zip.files[index], String)
    return XMLDict.xml_dict(subdata_xml)
end

function _correct_ids(data::Dict{String,<:Any})
    new_data = deepcopy(data)
    junction_names = keys(data["junction"])
    junction_mapping = Dict(k=>i for (i, k) in enumerate(junction_names))

    for (junction_name, junction) in data["junction"]
        i = junction_mapping[junction_name]
        new_data["junction"][string(i)] = junction
        new_data["junction"][string(i)]["id"] = i
        new_data["junction"][string(i)]["index"] = i
        delete!(new_data["junction"], junction_name)
    end

    for node_type in ["delivery", "receipt"]
        node_names = keys(data[node_type])

        for (i, node_name) in enumerate(node_names)
            new_data[node_type][string(i)] = data[node_type][node_name]
            new_data[node_type][string(i)]["id"] = i
            new_data[node_type][string(i)]["index"] = i
            new_data[node_type][string(i)]["junction_id"] = junction_mapping[node_name]
            delete!(new_data[node_type], node_name)
        end
    end

    for edge_type in ["compressor", "pipe", "resistor", "short_pipe", "valve"]
        edge_id = 1

        for (a, edge) in data[edge_type]
            fr_junction, to_junction = edge["fr_junction"], edge["to_junction"]
            edge["fr_junction"] = junction_mapping[fr_junction]
            edge["to_junction"] = junction_mapping[to_junction]
            edge["id"] = edge["index"] = edge_id
            new_data[edge_type][string(edge_id)] = edge
            delete!(new_data[edge_type], a)
            edge_id += 1
        end
    end

    return new_data
end

function _add_auxiliary_junctions!(junctions, compressors)
    new_junctions = Dict{String,Any}()
    new_compressors = Dict{String,Any}()

    for (a, compressor) in compressors
        junction_aux_name = a * "_aux_junction"
        fr_junction = junctions[compressor["fr_junction"]]
        to_junction = junctions[compressor["to_junction"]]
        junction_aux = deepcopy(fr_junction)

        compressor_reverse_name = a * "_reverse"
        compressor_reverse = deepcopy(compressor)
        compressor_reverse["fr_junction"] = junction_aux_name
        compressor["to_junction"] = junction_aux_name

        push!(new_junctions, junction_aux_name=>junction_aux)
        push!(new_compressors, compressor_reverse_name=>compressor_reverse)
    end

    junctions = _IM.update_data!(junctions, new_junctions)
    compressors = _IM.update_data!(compressors, new_compressors)
end

function _get_max_drive_power(drive)
    coeff_dict = filter(x -> occursin("power_fun_coeff", string(x.first)), drive)
    coeffs = [parse(Float64, coeff_dict["power_fun_coeff_" *
        string(i)][:value]) for i in 1:length(coeff_dict)]

    if length(coeffs) == 3
        return 1.0e3 * (coeffs[3] - (coeffs[2]^2 * inv(4.0 * coeffs[1])))
    elseif length(coeffs) == 9
        return 1.0e3 * coeffs[1]
    else
        Memento.error(_LOGGER, "Cannot compute maximum compressor power.")
    end
end

function _get_compressor_entry(compressor, stations)
    station = stations[findfirst(x -> compressor[:id] == x[:id], stations)]
    fr_junction, to_junction = compressor[:from], compressor[:to]
    inlet_p_min = parse(Float64, compressor["pressureInMin"][:value]) * 1.0e5
    inlet_p_max = parse(Float64, compressor["pressureOutMax"][:value]) * 1.0e5
    outlet_p_min = parse(Float64, compressor["pressureInMin"][:value]) * 1.0e5
    outlet_p_max = parse(Float64, compressor["pressureOutMax"][:value]) * 1.0e5
    flow_min = parse(Float64, compressor["flowMin"][:value]) * inv(3.6)
    flow_max = parse(Float64, compressor["flowMax"][:value]) * inv(3.6)

    if "diameterIn" in keys(compressor) && "diameterOut" in keys(compressor)
        diameter_in = parse(Float64, compressor["diameterIn"][:value]) * 1.0e-3
        diameter_out = parse(Float64, compressor["diameterOut"][:value]) * 1.0e-3
        diameter = 0.5 * (diameter_in + diameter_out) # TODO: Can we assume this?
    else
        diameter = 0.0 # TODO: What should be done, here?
    end

    c_ratio_min, c_ratio_max = 1.0, 5.0 # TODO: Can we derive better bounds?
    operating_cost = 10.0 # TODO: Is this derivable?
    power_max = 0.0 # Assume a pump's maximum power is zero to begin.

    for (comp_name, comp) in station["drives"]
        if comp isa Array
            for drive in comp
                power_max += abs(_get_max_drive_power(drive) * 1.0e3)
            end
        else
            power_max += abs(_get_max_drive_power(comp) * 1.0e3)
        end
    end

    return Dict{String,Any}("is_per_unit"=>false, "fr_junction"=>fr_junction,
        "to_junction"=>to_junction, "inlet_p_min"=>inlet_p_min,
        "inlet_p_max"=>inlet_p_max, "outlet_p_min"=>outlet_p_min,
        "outlet_p_max"=>outlet_p_max, "flow_min"=>flow_min,
        "flow_max"=>flow_max, "diameter"=>diameter,
        "is_per_unit"=>0, "directionality"=>2,
        "status"=>1, "is_si_units"=>1, "is_english_units"=>0,
        "c_ratio_min"=>c_ratio_min, "c_ratio_max"=>c_ratio_max,
        "power_max"=>power_max, "operating_cost"=>operating_cost)
end

function _get_delivery_entry(delivery)
    if delivery["flow"] isa Array
        min_id = findfirst(x -> x[:bound] == "lower", delivery["flow"])
        max_id = findfirst(x -> x[:bound] == "upper", delivery["flow"])
        withdrawal_min = parse(Float64, delivery["flow"][min_id][:value]) * inv(3.6)
        withdrawal_max = parse(Float64, delivery["flow"][max_id][:value]) * inv(3.6)
    else
        withdrawal_min = 0.0
        withdrawal_max = parse(Float64, delivery["flow"][:value]) * inv(3.6)
    end

    return Dict{String,Any}("withdrawal_min"=>withdrawal_min,
        "withdrawal_max"=>withdrawal_max, "withdrawal_nominal"=>withdrawal_max,
        "is_dispatchable"=>0, "is_per_unit"=>0, "status"=>1, "is_si_units"=>1,
        "is_english_units"=>0, "junction_id"=>delivery[:id])
end

function _get_junction_entry(junction)
    lat = parse(Float64, junction[:geoWGS84Lat])
    lon = parse(Float64, junction[:geoWGS84Long])
    p_min = parse(Float64, junction["pressureMin"][:value]) * 1.0e5
    p_max = parse(Float64, junction["pressureMax"][:value]) * 1.0e5

    return Dict{String,Any}("lat"=>lat, "lon"=>lon, "p_min"=>p_min, "p_max"=>p_max,
        "is_dispatchable"=>0, "is_per_unit"=>0, "status"=>1, "junction_type"=>0,
        "is_si_units"=>1, "is_english_units"=>0, "edi_id"=>junction[:id],
        "id"=>junction[:id], "index"=>junction[:id], "pipeline_id"=>"")
end

function _get_pipe_entry(pipe)
    fr_junction, to_junction = pipe[:from], pipe[:to]
    p_max = parse(Float64, pipe["pressureMax"][:value]) * 1.0e5
    diameter = parse(Float64, pipe["diameter"][:value]) * inv(1000.0)
    roughness = parse(Float64, pipe["roughness"][:value]) * inv(1000.0)
    length = parse(Float64, pipe["length"][:value]) * 1000.0
    friction_factor = (2.0 * log(3.7 * diameter / roughness))^(-2)

    return Dict{String,Any}("fr_junction"=>fr_junction, "to_junction"=>to_junction,
        "diameter"=>diameter, "length"=>length, "is_per_unit"=>0, "status"=>1,
        "is_si_units"=>1, "is_english_units"=>0, "p_min"=>0.0,
        "p_max"=>p_max, "friction_factor"=>friction_factor)
end

function _get_short_pipe_entry(short_pipe)
    fr_junction, to_junction = short_pipe[:from], short_pipe[:to]

    return Dict{String,Any}("fr_junction"=>fr_junction, "to_junction"=>to_junction,
        "is_bidirectional"=>1, "is_per_unit"=>0, "status"=>1, "is_si_units"=>1,
        "is_english_units"=>0)
end

function _get_receipt_entry(receipt)
    if receipt["flow"] isa Array
        min_id = findfirst(x -> x[:bound] == "lower", receipt["flow"])
        max_id = findfirst(x -> x[:bound] == "upper", receipt["flow"])
        injection_min = parse(Float64, receipt["flow"][min_id][:value]) * inv(3.6)
        injection_max = parse(Float64, receipt["flow"][max_id][:value]) * inv(3.6)
    else
        injection_min = 0.0
        injection_max = parse(Float64, receipt["flow"][:value]) * inv(3.6)
    end

    return Dict{String,Any}("injection_min"=>injection_min, "injection_max"=>injection_max,
        "injection_nominal"=>injection_max, "is_dispatchable"=>0, "is_per_unit"=>0,
        "status"=>1, "is_si_units"=>1, "is_english_units"=>0, "junction_id"=>receipt[:id])
end

function _get_resistor_entry(resistor)
    fr_junction, to_junction = resistor[:from], resistor[:to]
    diameter = parse(Float64, resistor["diameter"][:value]) * inv(1000.0)
    drag = parse(Float64, resistor["dragFactor"][:value])

    return Dict{String,Any}("fr_junction"=>fr_junction, "to_junction"=>to_junction,
        "drag"=>drag, "diameter"=>diameter, "is_per_unit"=>0, "status"=>1,
        "is_si_units"=>1, "is_english_units"=>0, "is_bidirectional"=>1)
end

function _get_valve_entry(valve)
    fr_junction, to_junction = valve[:from], valve[:to]

    return Dict{String,Any}("fr_junction"=>fr_junction, "to_junction"=>to_junction,
        "is_per_unit"=>0, "status"=>1, "is_si_units"=>1, "is_english_units"=>0,
        "is_bidirectional"=>1)
end

function _get_gaslib_compressors(topology, compressor_stations)
    compressors = topology["network"]["connections"]["compressorStation"]
    compressors = Dict{String,Any}(x[:id] => x for x in compressors)
    stations = compressor_stations["compressorStations"]["compressorStation"]
    return Dict{String,Any}(i => _get_compressor_entry(x, stations) for (i, x) in compressors)
end

function _get_gaslib_deliveries(topology, nomination)
    ids = [x[:id] for x in topology["network"]["nodes"]["sink"]]
    scenario = nomination["boundaryValue"]["scenario"]
    deliveries = Dict{String,Any}(x[:id] => x for x in scenario["node"] if x[:id] in ids)
    return Dict{String,Any}(i => _get_delivery_entry(x) for (i, x) in deliveries)
end

function _get_gaslib_junctions(topology)
    nodes = topology["network"]["nodes"]["innode"]
    nodes = vcat(nodes, topology["network"]["nodes"]["sink"])
    nodes = vcat(nodes, topology["network"]["nodes"]["source"])
    junctions = Dict{String,Any}(x[:id] => x for x in nodes)
    return Dict{String,Any}(i => _get_junction_entry(x) for (i, x) in junctions)
end

function _get_gaslib_short_pipes(topology)
    if "shortPipe" in keys(topology["network"]["connections"])
        short_pipes = topology["network"]["connections"]["shortPipe"]
        short_pipes = Dict{String,Any}(x[:id] => x for x in short_pipes)
        return Dict{String,Any}(i => _get_short_pipe_entry(x) for (i, x) in short_pipes)
    else
        return Dict{String,Any}()
    end
end

function _get_gaslib_pipes(topology)
    pipes = topology["network"]["connections"]["pipe"]
    pipes = Dict{String,Any}(x[:id] => x for x in pipes)
    return Dict{String,Any}(i => _get_pipe_entry(x) for (i, x) in pipes)
end

function _get_gaslib_receipts(topology, nomination)
    ids = [x[:id] for x in topology["network"]["nodes"]["source"]]
    scenario = nomination["boundaryValue"]["scenario"]
    receipts = Dict{String,Any}(x[:id] => x for x in scenario["node"] if x[:id] in ids)
    return Dict{String,Any}(i => _get_receipt_entry(x) for (i, x) in receipts)
end

function _get_gaslib_resistors(topology)
    if "resistor" in keys(topology["network"]["connections"])
        resistors = topology["network"]["connections"]["resistor"]
        resistors = Dict{String,Any}(x[:id] => x for x in resistors)
        return Dict{String,Any}(i => _get_resistor_entry(x) for (i, x) in resistors)
    else
        return Dict{String,Any}()
    end
end

function _get_gaslib_valves(topology)
    if "controlValve" in keys(topology["network"]["connections"])
        valves = topology["network"]["connections"]["controlValve"]
        valves = Dict{String,Any}(x[:id] => x for x in valves)
        return Dict{String,Any}(i => _get_valve_entry(x) for (i, x) in valves)
    else
        return Dict{String,Any}()
    end
end

function parse_gaslib(zip_path::Union{IO, String})
    # Read in the relevant data files as dictionaries.
    gaslib_zip = ZipFile.Reader(zip_path)
    topology = _read_gaslib_file(gaslib_zip, ".net")
    nominations = _read_gaslib_file(gaslib_zip, ".scn")
    compressor_stations = _read_gaslib_file(gaslib_zip, ".cs")

    # Create a dictionary for receipt components.
    compressors = _get_gaslib_compressors(topology, compressor_stations)
    deliveries = _get_gaslib_deliveries(topology, nominations)
    junctions = _get_gaslib_junctions(topology)
    receipts = _get_gaslib_receipts(topology, nominations)
    pipes = _get_gaslib_pipes(topology)
    resistors = _get_gaslib_resistors(topology)
    short_pipes = _get_gaslib_short_pipes(topology)
    valves = _get_gaslib_valves(topology)

    # Add auxiliary nodes for bidirectional compressors.
    _add_auxiliary_junctions!(junctions, compressors)

    # TODO: What is a general "sound_speed" that should be used?
    data = Dict{String,Any}("compressor"=>compressors, "delivery"=>deliveries,
        "junction"=>junctions, "receipt"=>receipts, "pipe"=>pipes,
        "resistor"=>resistors, "short_pipe"=>short_pipes, "valve"=>valves,
        "is_si_units"=>true, "per_unit"=>false, "sound_speed"=>312.8060)

    # Assign nodal IDs in place of string IDs.
    data = _correct_ids(data)

    # Return the dictionary.
    return data
end
