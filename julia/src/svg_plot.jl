_esc(s) = replace(replace(replace(string(s), "&"=>"&amp;"), "<"=>"&lt;"), ">"=>"&gt;")

function _urd_weight_label(weight)
    raw=string(weight)
    length(raw)<=48 && return raw
    weight isa LaurentPolynomial && return "Laurent polynomial ($(length(weight.terms)) terms)"
    weight isa LaurentRational && return "Laurent rational ($(length(weight.numerator.terms))/$(length(weight.denominator.terms)) terms)"
    "symbolic expression ($(length(raw)) chars)"
end

function _rounded_svg_path(points; radius=0.3)
    length(points) < 3 && return "M "*join(("$x $y" for (x,y) in points)," L ")
    commands=String["M $(points[1][1]) $(points[1][2])"]
    for index in 2:length(points)-1
        previous=complex(points[index-1]...); corner=complex(points[index]...); following=complex(points[index+1]...)
        incoming=corner-previous; outgoing=following-corner
        if abs(incoming) == 0 || abs(outgoing) == 0
            push!(commands,"L $(points[index][1]) $(points[index][2])")
            continue
        end
        local_radius=min(radius,0.45abs(incoming),0.45abs(outgoing))
        before=corner-local_radius*incoming/abs(incoming)
        after=corner+local_radius*outgoing/abs(outgoing)
        push!(commands,"L $(real(before)) $(imag(before))")
        push!(commands,"Q $(real(corner)) $(imag(corner)) $(real(after)) $(imag(after))")
    end
    push!(commands,"L $(points[end][1]) $(points[end][2])")
    join(commands," ")
end

function _write_svg(path, body; width=900, height=600)
    open(path, "w") do io
        println(io, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\" viewBox=\"0 0 $width $height\">")
        println(io, "<rect width=\"100%\" height=\"100%\" fill=\"#faf9f7\"/>")
        println(io, body)
        println(io, "</svg>")
    end
    abspath(path)
end

function plot_svg(d::URDiagram, path::AbstractString)
    # The MATLAB URDiagram.plot implementation is intentionally minimal: one
    # left-pointing horizontal quiver of length length(V)+1.
    body=String[]
    push!(body,"<line x1=\"780\" y1=\"300\" x2=\"120\" y2=\"300\" stroke=\"#0000ff\" stroke-width=\"4\"/>")
    push!(body,"<path d=\"M 150 280 L 120 300 L 150 320\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"4\"/>")
    _write_svg(path,join(body,"\n"))
end

"""Draw a readable UR-diagram diagnostic, including signed endpoints and weights.

`plot_svg(::URDiagram, ...)` deliberately preserves MATLAB's one-arrow
`URDiagram.plot`.  That rendering cannot validate a reduction: this companion
view makes the pairing represented by `V` and the current edge weights visible
without changing the source-compatible output.
"""
function plot_urd_structure_svg(d::URDiagram, path::AbstractString)
    validate_data(d)
    width,height=960,480
    margin=80
    count=max(length(d.V),2)
    x(position)=margin+(width-2margin)*(position-1)/(count-1)
    baseline=260
    body=String[
        "<text x=\"$(margin)\" y=\"42\" font-family=\"ui-monospace,monospace\" font-size=\"18\" fill=\"#252525\">UR diagram: signed endpoints, paired edges, weights</text>",
        "<line x1=\"$margin\" y1=\"$baseline\" x2=\"$(width-margin)\" y2=\"$baseline\" stroke=\"#a8a29e\" stroke-width=\"2\"/>"
    ]
    for edge in vec(d.E[1,:])
        positive=vertex_index(d,edge); negative=vertex_index(d,-edge)
        left,right=minmax(x(positive),x(negative))
        lane=42+28*((edge-1)%4)
        control=(left+right)/2
        color=positive<negative ? "#2563eb" : "#9333ea"
        push!(body,"<path d=\"M $left $baseline Q $control $(baseline-lane) $right $baseline\" fill=\"none\" stroke=\"$color\" stroke-width=\"3\"/>")
        label=_urd_weight_label(d.W[edge_index(d,edge)])
        push!(body,"<text x=\"$control\" y=\"$(baseline-lane-10)\" text-anchor=\"middle\" font-family=\"ui-monospace,monospace\" font-size=\"14\" fill=\"#292524\">e$(edge): $(_esc(label))</text>")
    end
    for (position,vertex) in enumerate(d.V)
        color=vertex>0 ? "#1d4ed8" : "#7f1d1d"
        push!(body,"<circle cx=\"$(x(position))\" cy=\"$baseline\" r=\"15\" fill=\"#fff\" stroke=\"$color\" stroke-width=\"3\"/>")
        push!(body,"<text x=\"$(x(position))\" y=\"$(baseline+5)\" text-anchor=\"middle\" font-family=\"ui-monospace,monospace\" font-size=\"14\" fill=\"$color\">$vertex</text>")
    end
    _write_svg(path,join(body,"\n");width=width,height=height)
end

_crossing_sign(slot::Integer) = slot in (1,3) ? -1 : 1
_is_outgoing(slot::Integer,orientation::Integer) =
    slot == 3 || (slot == 4 && orientation == -1) || (slot == 2 && orientation != -1)

function _directed_piece(edge::EdgeLayout,orientation)
    if _is_outgoing(edge.tail_slot,orientation[edge.tail])
        source=_crossing_sign(edge.tail_slot)*edge.tail
        target=_crossing_sign(edge.head_slot)*edge.head
        return (source=source,target=target,points=edge.points,id=edge.id)
    end
    source=_crossing_sign(edge.head_slot)*edge.head
    target=_crossing_sign(edge.tail_slot)*edge.tail
    (source=source,target=target,points=reverse(edge.points),id=edge.id)
end

"""Aggregate planar PD edges into MATLAB REdgeTable arcs between real crossings."""
function _real_arc_layouts(v::VirtualLink,layout::LinkLayout)
    code = v.real_only ? first(planarize_virtual_gauss(v.gauss,v.orientation)) : v.gauss
    pieces=Dict{Tuple{Int,Int},Vector{NamedTuple}}()
    for edge in layout.edges
        piece=_directed_piece(edge,layout.orientation)
        push!(get!(pieces,(piece.source,piece.target),NamedTuple[]),piece)
    end
    arcs=NamedTuple[]
    arc_id=0
    for (component_index,component) in enumerate(code)
        component=filter(!=(0),component)
        isempty(component) && continue
        sequence=NamedTuple[]
        for index in eachindex(component)
            pair=(component[index],component[mod1(index+1,length(component))])
            candidates=get(pieces,pair,NamedTuple[])
            isempty(candidates) && throw(ArgumentError("no directed PD edge for Gauss pair $pair"))
            push!(sequence,popfirst!(candidates))
        end
        real_positions=findall(index -> layout.orientation[abs(component[index])] != 0,eachindex(component))
        for (weight_index,target_position) in enumerate(real_positions)
            previous_position=real_positions[mod1(weight_index-1,length(real_positions))]
            indices=Int[]
            cursor=previous_position
            while cursor != target_position
                push!(indices,cursor)
                cursor=mod1(cursor+1,length(component))
            end
            points=ComplexF64[]
            for index in indices
                append!(points,isempty(points) ? sequence[index].points : sequence[index].points[2:end])
            end
            arc_id += 1
            component_weights=component_index <= length(v.weights) ? v.weights[component_index] : Float64[]
            weight=weight_index <= length(component_weights) ? component_weights[weight_index] : NaN
            push!(arcs,(id=arc_id,points=points,source=component[previous_position],
                        target=component[target_position],weight=weight))
        end
    end
    arcs
end

"""Stable MATLAB `REdgeTable` arc metadata, without invoking the layout solver."""
function _replay_arc_metadata(v::VirtualLink)
    [(id=edge.id,source=edge.crossing[1],target=edge.crossing[2],weight=edge.weight)
     for edge in real_edge_table(v)]
end

"""Render a VirtualLink using either the Julia layout or replayed MATLAB arcs.

`replay_positions`, when supplied, is one polyline per real arc in the stable
order returned by `_real_arc_layouts`.  It is the direct representation of
MATLAB `REdgeTable.Position`; retaining it avoids treating a different optimal
orthogonal embedding as an equivalent visual result.
"""
function plot_svg(v::VirtualLink, path::AbstractString; bending_numbers=nothing,
                  edge_labels=:weight,replay_positions=nothing)
    body=String[]
    if isempty(v.orientation)
        push!(body,"<circle cx=\"450\" cy=\"300\" r=\"145\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"4\"/>")
        push!(body,"<path d=\"M 578 300 L 595 274 L 612 300\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"4\"/>")
        return _write_svg(path,join(body,"\n"))
    end
    layout=nothing
    arcs=if replay_positions === nothing
        layout=orthogonal_layout(v;bending_numbers=bending_numbers)
        _real_arc_layouts(v,layout)
    else
        metadata=_replay_arc_metadata(v)
        length(replay_positions)==length(metadata) || throw(ArgumentError("one MATLAB Position polyline is required for each real arc"))
        [merge(arc,(points=ComplexF64.(collect(points)),)) for (arc,points) in zip(metadata,replay_positions)]
    end
    if replay_positions !== nothing
        all(length(arc.points)>=2 for arc in arcs) || throw(ArgumentError("every MATLAB Position polyline requires at least two points"))
    end
    all_points=reduce(vcat,(arc.points for arc in arcs))
    xmin,xmax=extrema(real.(all_points)); ymin,ymax=extrema(imag.(all_points))
    sx=760/max(xmax-xmin,1); sy=480/max(ymax-ymin,1); scale=min(sx,sy)
    project(z)=((real(z)-xmin)*scale+70,(ymax-imag(z))*scale+60)

    for arc in arcs
        points=copy(arc.points)
        # Match MATLAB's crossing gap: under-strand endpoints are shortened.
        if arc.source < 0 && length(points)>1
            points[1] += 0.20*(points[2]-points[1])/abs(points[2]-points[1])
        end
        if arc.target < 0 && length(points)>1
            points[end] += 0.20*(points[end-1]-points[end])/abs(points[end-1]-points[end])
        end
        projected=[project(point) for point in points]
        data=_rounded_svg_path(projected;radius=0.3scale)
        push!(body,"<path d=\"$data\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"4\" stroke-linejoin=\"round\"/>")

        middle=max(1,min(length(projected)-1,fld(length(projected),2)))
        p1=points[middle]; p2=points[middle+1]; midpoint=(p1+p2)/2
        unit=(p2-p1)/abs(p2-p1); left=midpoint-unit*0.16+unit*im*0.09; right=midpoint-unit*0.16-unit*im*0.09
        arrow=[project(left),project(midpoint),project(right)]
        arrow_data="M $(arrow[1][1]) $(arrow[1][2]) L $(arrow[2][1]) $(arrow[2][2]) L $(arrow[3][1]) $(arrow[3][2])"
        push!(body,"<path d=\"$arrow_data\" fill=\"none\" stroke=\"#0000ff\" stroke-width=\"3\"/>")
        if v.is_weighted || edge_labels in (:id,:all)
            parts=String[]
            edge_labels in (:id,:all) && push!(parts,"E$(arc.id)")
            v.is_weighted && edge_labels in (:weight,:all) && !isnan(arc.weight) &&
                push!(parts,string(arc.weight))
            label=join(parts,", ")
            mx,my=project(midpoint)
            push!(body,"<text x=\"$(mx+5)\" y=\"$(my-5)\" font-family=\"sans-serif\" font-size=\"17\">$(_esc(label))</text>")
        end
    end
    if replay_positions === nothing
        real_index=0
        for (index,point) in enumerate(layout.crossings)
            layout.orientation[index] == 0 && continue
            real_index += 1
            x,y=project(point)
            push!(body,"<circle cx=\"$x\" cy=\"$y\" r=\"5\" fill=\"#0000ff\"/>")
            push!(body,"<text x=\"$(x+8)\" y=\"$(y-8)\" font-family=\"sans-serif\" font-size=\"17\">V$real_index</text>")
        end
    else
        crossings=Dict{Int,ComplexF64}()
        for arc in arcs
            for (crossing,point) in ((abs(arc.source),first(arc.points)),(abs(arc.target),last(arc.points)))
                haskey(crossings,crossing) && !isapprox(crossings[crossing],point;atol=1e-9) &&
                    throw(ArgumentError("MATLAB Position polylines disagree at crossing $crossing"))
                crossings[crossing]=point
            end
        end
        for (real_index,crossing) in enumerate(sort!(collect(keys(crossings))))
            x,y=project(crossings[crossing])
            push!(body,"<circle cx=\"$x\" cy=\"$y\" r=\"5\" fill=\"#0000ff\"/>")
            push!(body,"<text x=\"$(x+8)\" y=\"$(y-8)\" font-family=\"sans-serif\" font-size=\"17\">V$real_index</text>")
        end
    end
    _write_svg(path,join(body,"\n"))
end
