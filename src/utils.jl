export removetransient, filter_repeats

function removetransient(transient, dt, v...)
    exclude = round(Int, transient / dt)
    if isempty(v)
        return v
    end
    num_elem, min_size = extrema(size(u, ndims(u)) for u in v)
    @assert num_elem == min_size
    if exclude >= num_elem
        return (empty(u) for u ∈ v)
    end
    return (selectdim(u, ndims(u), exclude:size(u, ndims(u))) for u ∈ v)
end


function filter_repeats(vs)
    [vs[i] for i in eachindex(vs) if i == 1 || vs[i-1] != vs[i]]
end