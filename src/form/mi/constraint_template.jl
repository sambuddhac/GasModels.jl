######################################################################################
# Constraints associated witn cutting planes on the direction variables
######################################################################################

" Template: Constraints for ensuring that at least one direction is set to take flow away from a junction (typically used on source nodes) "
function constraint_source_flow(gm::GenericGasModel{T}, n::Int, i) where T <: AbstractMIForms
    f_pipes          = ref(gm,n,:f_pipes,i)
    t_pipes          = ref(gm,n,:t_pipes,i)
    f_compressors    = ref(gm,n,:f_compressors,i)
    t_compressors    = ref(gm,n,:t_compressors,i)
    f_resistors      = ref(gm,n,:f_resistors,i)
    t_resistors      = ref(gm,n,:t_resistors,i)
    f_short_pipes    = ref(gm,n,:f_short_pipes,i)
    t_short_pipes    = ref(gm,n,:t_short_pipes,i)
    f_valves         = ref(gm,n,:f_valves,i)
    t_valves         = ref(gm,n,:t_valves,i)
    f_control_valves = ref(gm,n,:f_control_valves,i)
    t_control_valves = ref(gm,n,:t_control_valves,i)

    constraint_source_flow(gm, n, i, f_pipes, t_pipes, f_compressors, t_compressors, f_resistors, t_resistors, f_short_pipes, t_short_pipes, f_valves, t_valves, f_control_valves, t_control_valves)
end
constraint_source_flow(gm::GenericGasModel, i::Int) = constraint_source_flow(gm, gm.cnw, i)

" Template: Constraints for ensuring that at least one direction is set to take flow away from a junction (typically used on source nodes) "
function constraint_source_flow_ne(gm::GenericGasModel{T}, n::Int, i) where T <: AbstractMIForms
    f_pipes          = ref(gm,n,:f_pipes,i)
    t_pipes          = ref(gm,n,:t_pipes,i)
    f_compressors    = ref(gm,n,:f_compressors,i)
    t_compressors    = ref(gm,n,:t_compressors,i)
    f_resistors      = ref(gm,n,:f_resistors,i)
    t_resistors      = ref(gm,n,:t_resistors,i)
    f_short_pipes    = ref(gm,n,:f_short_pipes,i)
    t_short_pipes    = ref(gm,n,:t_short_pipes,i)
    f_valves         = ref(gm,n,:f_valves,i)
    t_valves         = ref(gm,n,:t_valves,i)
    f_control_valves = ref(gm,n,:f_control_valves,i)
    t_control_valves = ref(gm,n,:t_control_valves,i)
    f_ne_pipes       = ref(gm,n,:f_ne_pipes,i)
    t_ne_pipes       = ref(gm,n,:t_ne_pipes,i)
    f_ne_compressors = ref(gm,n,:f_ne_compressors,i)
    t_ne_compressors = ref(gm,n,:t_ne_compressors,i)

    constraint_source_flow_ne(gm, n, i, f_pipes, t_pipes, f_compressors, t_compressors, f_resistors, t_resistors, f_short_pipes, t_short_pipes, f_valves, t_valves, f_control_valves, t_control_valves, f_ne_pipes, t_ne_pipes, f_ne_compressors, t_ne_compressors)
end
constraint_source_flow_ne(gm::GenericGasModel, i::Int) = constraint_source_flow_ne(gm, gm.cnw, i)

" Template: Constraints for ensuring that at least one direction is set to take flow away from a junction (typically used on sink nodes) "
function constraint_sink_flow(gm::GenericGasModel{T}, n::Int, i) where T <: AbstractMIForms
    f_pipes          = ref(gm,n,:f_pipes,i)
    t_pipes          = ref(gm,n,:t_pipes,i)
    f_compressors    = ref(gm,n,:f_compressors,i)
    t_compressors    = ref(gm,n,:t_compressors,i)
    f_resistors      = ref(gm,n,:f_resistors,i)
    t_resistors      = ref(gm,n,:t_resistors,i)
    f_short_pipes    = ref(gm,n,:f_short_pipes,i)
    t_short_pipes    = ref(gm,n,:t_short_pipes,i)
    f_valves         = ref(gm,n,:f_valves,i)
    t_valves         = ref(gm,n,:t_valves,i)
    f_control_valves = ref(gm,n,:f_control_valves,i)
    t_control_valves = ref(gm,n,:t_control_valves,i)

    constraint_sink_flow(gm, n, i, f_pipes, t_pipes, f_compressors, t_compressors, f_resistors, t_resistors, f_short_pipes, t_short_pipes, f_valves, t_valves, f_control_valves, t_control_valves)
end
constraint_sink_flow(gm::GenericGasModel, i::Int) = constraint_sink_flow(gm, gm.cnw, i)

" Template: Constraints for ensuring that at least one direction is set to take flow away from a junction (typically used on sink nodes) "
function constraint_sink_flow_ne(gm::GenericGasModel{T}, n::Int, i) where T <: AbstractMIForms
    f_pipes          = ref(gm,n,:f_pipes,i)
    t_pipes          = ref(gm,n,:t_pipes,i)
    f_compressors    = ref(gm,n,:f_compressors,i)
    t_compressors    = ref(gm,n,:t_compressors,i)
    f_resistors      = ref(gm,n,:f_resistors,i)
    t_resistors      = ref(gm,n,:t_resistors,i)
    f_short_pipes    = ref(gm,n,:f_short_pipes,i)
    t_short_pipes    = ref(gm,n,:t_short_pipes,i)
    f_valves         = ref(gm,n,:f_valves,i)
    t_valves         = ref(gm,n,:t_valves,i)
    f_control_valves = ref(gm,n,:f_control_valves,i)
    t_control_valves = ref(gm,n,:t_control_valves,i)
    f_ne_pipes       = ref(gm,n,:f_ne_pipes,i)
    t_ne_pipes       = ref(gm,n,:t_ne_pipes,i)
    f_ne_compressors = ref(gm,n,:f_ne_compressors,i)
    t_ne_compressors = ref(gm,n,:t_ne_compressors,i)

    constraint_sink_flow_ne(gm, n, i, f_pipes, t_pipes, f_compressors, t_compressors, f_resistors, t_resistors, f_short_pipes, t_short_pipes, f_valves, t_valves, f_control_valves, t_control_valves, f_ne_pipes, t_ne_pipes, f_ne_compressors, t_ne_compressors)
end
constraint_sink_flow_ne(gm::GenericGasModel, i::Int) = constraint_sink_flow_ne(gm, gm.cnw, i)

" Template: Constraints to ensure that flow is the same direction through a node with degree 2 and no production or consumption "
function constraint_conserve_flow(gm::GenericGasModel{T}, n::Int, idx) where T <: AbstractMIForms
    f_pipes          = Dict(i => ref(gm,n,:pipe,i)["t_junction"] for i in ref(gm,n,:f_pipes,idx))
    t_pipes          = Dict(i => ref(gm,n,:pipe,i)["f_junction"] for i in ref(gm,n,:t_pipes,idx))
    f_compressors    = Dict(i => ref(gm,n,:compressor,i)["t_junction"] for i in ref(gm,n,:f_compressors,idx))
    t_compressors    = Dict(i => ref(gm,n,:compressor,i)["f_junction"] for i in ref(gm,n,:t_compressors,idx))
    f_resistors      = Dict(i => ref(gm,n,:resistor,i)["t_junction"] for i in ref(gm,n,:f_resistors,idx))
    t_resistors      = Dict(i => ref(gm,n,:resistor,i)["f_junction"] for i in ref(gm,n,:t_resistors,idx))
    f_short_pipes    = Dict(i => ref(gm,n,:short_pipe,i)["t_junction"] for i in ref(gm,n,:f_short_pipes,idx))
    t_short_pipes    = Dict(i => ref(gm,n,:short_pipe,i)["f_junction"] for i in ref(gm,n,:t_short_pipes,idx))
    f_valves         = Dict(i => ref(gm,n,:valve,i)["t_junction"] for i in ref(gm,n,:f_valves,idx))
    t_valves         = Dict(i => ref(gm,n,:valve,i)["f_junction"] for i in ref(gm,n,:t_valves,idx))
    f_control_valves = Dict(i => ref(gm,n,:control_valve,i)["t_junction"] for i in ref(gm,n,:f_control_valves,idx))
    t_control_valves = Dict(i => ref(gm,n,:control_valve,i)["f_junction"] for i in ref(gm,n,:t_control_valves,idx))

    constraint_conserve_flow(gm, n, idx, f_pipes, t_pipes, f_compressors, t_compressors, f_resistors, t_resistors, f_short_pipes, t_short_pipes, f_valves, t_valves, f_control_valves, t_control_valves)
end
constraint_conserve_flow(gm::GenericGasModel, i::Int) = constraint_conserve_flow(gm, gm.cnw, i)

" Template: Constraints to ensure that flow is the same direction through a node with degree 2 and no production or consumption "
function constraint_conserve_flow_ne(gm::GenericGasModel{T}, n::Int, idx) where T <: AbstractMIForms
    f_pipes          = Dict(i => ref(gm,n,:pipe,i)["t_junction"] for i in ref(gm,n,:f_pipes,idx))
    t_pipes          = Dict(i => ref(gm,n,:pipe,i)["f_junction"] for i in ref(gm,n,:t_pipes,idx))
    f_compressors    = Dict(i => ref(gm,n,:compressor,i)["t_junction"] for i in ref(gm,n,:f_compressors,idx))
    t_compressors    = Dict(i => ref(gm,n,:compressor,i)["f_junction"] for i in ref(gm,n,:t_compressors,idx))
    f_resistors      = Dict(i => ref(gm,n,:resistor,i)["t_junction"] for i in ref(gm,n,:f_resistors,idx))
    t_resistors      = Dict(i => ref(gm,n,:resistor,i)["f_junction"] for i in ref(gm,n,:t_resistors,idx))
    f_short_pipes    = Dict(i => ref(gm,n,:short_pipe,i)["t_junction"] for i in ref(gm,n,:f_short_pipes,idx))
    t_short_pipes    = Dict(i => ref(gm,n,:short_pipe,i)["f_junction"] for i in ref(gm,n,:t_short_pipes,idx))
    f_valves         = Dict(i => ref(gm,n,:valve,i)["t_junction"] for i in ref(gm,n,:f_valves,idx))
    t_valves         = Dict(i => ref(gm,n,:valve,i)["f_junction"] for i in ref(gm,n,:t_valves,idx))
    f_control_valves = Dict(i => ref(gm,n,:control_valve,i)["t_junction"] for i in ref(gm,n,:f_control_valves,idx))
    t_control_valves = Dict(i => ref(gm,n,:control_valve,i)["f_junction"] for i in ref(gm,n,:t_control_valves,idx))
    f_ne_pipes       = Dict(i => ref(gm,n,:ne_pipe,i)["t_junction"] for i in ref(gm,n,:f_ne_pipes,idx))
    t_ne_pipes       = Dict(i => ref(gm,n,:ne_pipe,i)["f_junction"] for i in ref(gm,n,:t_ne_pipes,idx))
    f_ne_compressors = Dict(i => ref(gm,n,:ne_compressor,i)["t_junction"] for i in ref(gm,n,:f_ne_compressors,idx))
    t_ne_compressors = Dict(i => ref(gm,n,:ne_compressor,i)["f_junction"] for i in ref(gm,n,:t_ne_compressors,idx))

    constraint_conserve_flow_ne(gm, n, idx, f_pipes, t_pipes, f_compressors, t_compressors, f_resistors,
                                    t_resistors, f_short_pipes, t_short_pipes, f_valves, t_valves, f_control_valves,
                                    t_control_valves, f_ne_pipes, t_ne_pipes, f_ne_compressors, t_ne_compressors)

end
constraint_conserve_flow_ne(gm::GenericGasModel, i::Int) = constraint_conserve_flow_ne(gm, gm.cnw, i)

" Template: Constraints which ensure that parallel lines have flow in the same direction "
function constraint_parallel_flow(gm::GenericGasModel{T}, n::Int, idx) where T <: AbstractMIForms
    connection = ref(gm,n,:connection,idx)
    connections = ref(gm,n,:connection)

    i = min(connection["f_junction"], connection["t_junction"])
    j = max(connection["f_junction"], connection["t_junction"])

    parallel_connections = ref(gm,n,:parallel_connections,(i,j))

    f_connections = filter(i -> connections[i]["f_junction"] == connection["f_junction"], parallel_connections)
    t_connections = filter(i -> connections[i]["f_junction"] != connection["f_junction"], parallel_connections)

    if length(parallel_connections) <= 1
        return nothing
    end

    constraint_parallel_flow(gm, n, idx, i, j, f_connections, t_connections)
end
constraint_parallel_flow(gm::GenericGasModel, i::Int) = constraint_parallel_flow(gm, gm.cnw, i)

" Template: Constraints which ensure that parallel lines have flow in the same direction "
function constraint_parallel_flow_ne(gm::GenericGasModel{T}, n::Int, idx) where T <: AbstractMIForms
    connection = haskey(ref(gm,n,:connection), idx) ? ref(gm,n,:connection)[idx] : ref(gm,n,:ne_connection)[idx]
    connections = ref(gm,n,:connection)
    ne_connections = ref(gm,n,:ne_connection)

    i = min(connection["f_junction"], connection["t_junction"])
    j = max(connection["f_junction"], connection["t_junction"])

    all_parallel_connections = ref(gm,n,:parallel_ne_connections, (i,j))
    parallel_connections = ref(gm,n,:parallel_connections, (i,j))

    if length(all_parallel_connections) <= 1
        return nothing
    end

    f_connections = filter(i -> connections[i]["f_junction"] == connection["f_junction"], intersect(all_parallel_connections, parallel_connections))
    t_connections = filter(i -> connections[i]["f_junction"] != connection["f_junction"], intersect(all_parallel_connections, parallel_connections))
    f_connections_ne = filter(i -> ne_connections[i]["f_junction"] == connection["f_junction"], setdiff(all_parallel_connections, parallel_connections))
    t_connections_ne = filter(i -> ne_connections[i]["f_junction"] != connection["f_junction"], setdiff(all_parallel_connections, parallel_connections))

    constraint_parallel_flow_ne(gm, n, idx, i, j, f_connections, t_connections, f_connections_ne, t_connections_ne)
end
constraint_parallel_flow_ne(gm::GenericGasModel, i::Int) = constraint_parallel_flow_ne(gm, gm.cnw, i)

" Template: Constraints which ensure that parallel lines have flow in the same direction - customized for ne_pipe"
function constraint_ne_pipe_parallel_flow(gm::GenericGasModel{T}, n::Int, idx) where T <: AbstractMIForms
    pipe = ref(gm,n,:ne_pipe, idx)
    num_connections, aligned_pipes, opposite_pipes, aligned_compressors, opposite_compressors,
           aligned_resistors, opposite_resistors, aligned_short_pipes, opposite_short_pipes, aligned_valves, opposite_valves,
           aligned_control_valves, opposite_control_valves, aligned_ne_pipes, opposite_ne_pipes, aligned_ne_compressors, opposite_ne_compressors =
           calc_parallel_ne_connections(gm, n, pipe)

    if num_connections <= 1
        return nothing
    end

    constraint_ne_pipe_parallel_flow(gm, n, idx, num_connections, aligned_pipes, opposite_pipes, aligned_compressors, opposite_compressors,
                                     aligned_resistors, opposite_resistors, aligned_short_pipes, opposite_short_pipes, aligned_valves, opposite_valves,
                                     aligned_control_valves, opposite_control_valves, aligned_ne_pipes, opposite_ne_pipes, aligned_ne_compressors, opposite_ne_compressors)
end
constraint_ne_pipe_parallel_flow(gm::GenericGasModel, i::Int) = constraint_ne_pipe_parallel_flow(gm, gm.cnw, i)

" Template: Constraints which ensure that parallel lines have flow in the same direction - customized for ne_compressor"
function constraint_ne_compressor_parallel_flow(gm::GenericGasModel{T}, n::Int, idx) where T <: AbstractMIForms
    compressor = ref(gm,n,:ne_compressor, idx)
    num_connections, aligned_pipes, opposite_pipes, aligned_compressors, opposite_compressors,
           aligned_resistors, opposite_resistors, aligned_short_pipes, opposite_short_pipes, aligned_valves, opposite_valves,
           aligned_control_valves, opposite_control_valves, aligned_ne_pipes, opposite_ne_pipes, aligned_ne_compressors, opposite_ne_compressors =
           calc_parallel_ne_connections(gm, n, compressor)

    if num_connections <= 1
        return nothing
    end

    constraint_ne_compressor_parallel_flow(gm, n, idx, num_connections, aligned_pipes, opposite_pipes, aligned_compressors, opposite_compressors,
                                     aligned_resistors, opposite_resistors, aligned_short_pipes, opposite_short_pipes, aligned_valves, opposite_valves,
                                     aligned_control_valves, opposite_control_valves, aligned_ne_pipes, opposite_ne_pipes, aligned_ne_compressors, opposite_ne_compressors)
end
constraint_ne_compressor_parallel_flow(gm::GenericGasModel, i::Int) = constraint_ne_compressor_parallel_flow(gm, gm.cnw, i)
