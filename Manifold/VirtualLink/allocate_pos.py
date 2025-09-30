def allocate_positions(self,bendingNumbers):
    # import numpy 
    # import mip
    # pd_code = self.pd_code()
    # color=1
    # if type(color) is not dict:
    #     coloring = {int(i): color for i in set(flatten(pd_code))}
    # else:
    #     from sage.plot.colors import rainbow
    #     ncolors = max([int(i) for i in color.values()]) + 1
    #     arcs = self.arcs()
    #     rainb = rainbow(ncolors)
    #     coloring = {int(i): rainb[color[tuple(j)]] for j in arcs for i in j}
    # comp = self._isolated_components()
    # # Handle isolated components individually
    # if len(comp) > 1:
    #     L1 = Link(comp[0])
    #     L2 = Link(flatten(comp[1:], max_level=1))
    #     P1 = L1.plot(gap, **kwargs)
    #     P2 = L2.plot(gap, **kwargs)
    #     xtra = P1.get_minmax_data()['xmax'] + component_gap - P2.get_minmax_data()['xmin']
    #     for P in P2:
    #         if hasattr(P, 'path'):
    #             for p in P.path[0]:
    #                 p[0] += xtra
    #             for p in P.vertices:
    #                 p[0] += xtra
    #         else:
    #             P.xdata = [p + xtra for p in P.xdata]
    #     return P1 + P2
    # 
    # if 'axes' not in kwargs:
    #     kwargs['axes'] = False
    # if 'aspect_ratio' not in kwargs:
    #     kwargs['aspect_ratio'] = 1
    # 
    # from sage.plot.line import line
    # from sage.plot.bezier_path import bezier_path
    # from sage.plot.circle import circle
    # 
    # # Special case for the unknot
    # if not pd_code:
    #     return circle((0, 0), ZZ.one() / ZZ(2), color=color, **kwargs)
    
    # The idea is the same followed in spherogram, but using MLP instead of
    # network flows.
    # We start by computing a way to bend the edges left or right
    # such that the resulting regions are in fact closed regions
    # with straight angles, and using the minimal number of bends.

    # input interface
    # regions = sorted(self.regions(), key=len)
    # edges = list(set(flatten(pd_code)))
    # edges.sort()

    pd_code=self.pd_code()
    # print(pd_code)
    orientation=self.orientation()
    # print(orientation)
    edges = list(set(flatten(pd_code)))
    regions = sorted(self.regions(), key=len)
    edges.sort()
    MLP = MixedIntegerLinearProgram(maximization=False)
    # v will be the list of variables in the MLP problem. There will be
    # two variables for each edge counting the number of bendings needed.
    # The one with even index corresponds to the flow of this number from
    # the left-hand-side region to the right-hand-side region if the edge
    # is positive oriented. The one with odd index corresponds to the
    # flow in the opposite direction. For a negative oriented edge the
    # same is true but with exchanged directions. At the end, since we
    # are minimizing the total, only one of each will be nonzero.
    v = MLP.new_variable(nonnegative=True, integer=True)
    
    def flow_from_source(e):
        r"""
        Return the flow variable from the source.
        """
        if e > 0:
            return v[2 * edges.index(e)]
        else:
            return v[2 * edges.index(-e) + 1]
    
    def flow_to_sink(e):
        r"""
        Return the flow variable to the sink.
        """
        return flow_from_source(-e)
    
    # one condition for each region
    lr = len(regions)
    for i in range(lr):
        r = regions[i]
        if i < lr - 1:
            # capacity of interior region, sink if positive, source if negative
            capacity = len(r) - 4
        else:
            # capacity of exterior region, only sink (added to fix :issue:`37587`).
            capacity = len(r) + 4
        flow = sum(flow_to_sink(e) - flow_from_source(e) for e in r)
        MLP.add_constraint(flow == capacity)  # exterior region only sink
    
    MLP.set_objective(MLP.sum(v.values()))
    MLP.solve()
    # we store the result in a vector s packing right bends as negative left ones
    # values = MLP.get_values(v, convert=ZZ, tolerance=1e-3)
    values = MLP.get_values(v)
    s = None
    # if bendingNumbers is empty list, we use the MLP solution
    if not bendingNumbers:
        s = [values[2 * i] - values[2 * i + 1] for i in range(len(edges))]
    else:
        s = bendingNumbers
    # print("bending numbers", s)
    
    # segments represents the different parts of the previous edges after bending
    segments = {e: [(e, i) for i in range(abs(s[edges.index(e)]) + 1)]
                for e in edges}
    pieces = {tuple(i): [i] for j in segments.values() for i in j}
    nregions = []
    for r in regions[:-1]:  # interior regions
        nregion = []
        for e in r:
            if e > 0:
                rev = segments[e][:-1]
                sig = sign(s[edges.index(e)])
                nregion += [[a, sig] for a in rev]
                nregion.append([segments[e][-1], 1])
            else:
                rev = segments[-e][1:]
                rev.reverse()
                sig = sign(s[edges.index(-e)])
                nregion += [[a, -sig] for a in rev]
                nregion.append([segments[-e][0], 1])
        nregions.append(nregion)
    N = max(segments) + 1
    segments = [i for j in segments.values() for i in j]
    badregions = [nr for nr in nregions if any(-1 == x[1] for x in nr)]
    while badregions:
        badregion = badregions[0]
        a = 0
        while badregion[a][1] != -1:
            a += 1
        c = -1
        b = a
        while c != 2:
            if b == len(badregion) - 1:
                b = 0
            else:
                b += 1
            c += badregion[b][1]
        otherregion = [nr for nr in nregions
                       if any(badregion[b][0] == x[0] for x in nr)]
        if len(otherregion) == 1:
            otherregion = None
        elif otherregion[0] == badregion:
            otherregion = otherregion[1]
        else:
            otherregion = otherregion[0]
        N1 = N
        N = N + 2
        N2 = N1 + 1
        segments.append(N1)
        segments.append(N2)
        if type(badregion[b][0]) in (int, Integer):
            segmenttoadd = [x for x in pieces
                            if badregion[b][0] in pieces[x]]
            if len(segmenttoadd) > 0:
                pieces[segmenttoadd[0]].append(N2)
        else:
            pieces[tuple(badregion[b][0])].append(N2)
    
        if a < b:
            r1 = badregion[:a] + [[badregion[a][0], 0], [N1, 1]] + badregion[b:]
            r2 = badregion[a + 1:b] + [[N2, 1], [N1, 1]]
        else:
            r1 = badregion[b:a] + [[badregion[a][0], 0], [N1, 1]]
            r2 = badregion[:b] + [[N2, 1], [N1, 1]] + badregion[a + 1:]
    
        if otherregion:
            c = [x for x in otherregion if badregion[b][0] == x[0]]
            c = otherregion.index(c[0])
            otherregion.insert(c + 1, [N2, otherregion[c][1]])
            otherregion[c][1] = 0
        nregions.remove(badregion)
        nregions.append(r1)
        nregions.append(r2)
        badregions = [nr for nr in nregions if any(x[1] == -1 for x in nr)]
    # print(pieces)
    MLP = MixedIntegerLinearProgram(maximization=False)
    v = MLP.new_variable(nonnegative=True, integer=True)
    for e in segments:
        MLP.set_min(v[e], 1)
    for r in nregions:
        horp = []
        horm = []
        verp = []
        verm = []
        direction = 0
        for se in r:
            if direction % 4 == 0:
                horp.append(v[se[0]])
            elif direction == 1:
                verp.append(v[se[0]])
            elif direction == 2:
                horm.append(v[se[0]])
            elif direction == 3:
                verm.append(v[se[0]])
            if se[1] == 1:
                direction += 1
        MLP.add_constraint(MLP.sum(horp) - MLP.sum(horm) == 0)
        MLP.add_constraint(MLP.sum(verp) - MLP.sum(verm) == 0)
    MLP.set_objective(MLP.sum(v.values()))
    MLP.solve()
    v = MLP.get_values(v)
    lengths = {piece: sum(v[a] for a in pieces[piece]) for piece in pieces}
    # image = line([], **kwargs)
    crossings = {tuple(pd_code[0]): (0, 0, 0)}
    availables = pd_code[1:]
    used_edges = []
    coords = []
    edgeID=[]
    while len(used_edges) < len(edges):
        cross_keys = list(crossings.keys())
        i = 0
        j = 0
        while cross_keys[i][j] in used_edges:
            if j < 3:
                j += 1
            else:
                j = 0
                i += 1
        c = cross_keys[i]
        e = c[j]
        used_edges.append(e)
        ce=e
        direction = (crossings[c][2] + c.index(e)) % 4
        orien = orientation[pd_code.index(list(c))]
        if s[edges.index(e)] < 0:
            turn = -1
        else:
            turn = 1
        lengthse = [lengths[(e, k)] for k in range(abs(s[edges.index(e)]) + 1)]
        if c.index(e) == 0 or (c.index(e) == 3 and orien == 1) or (c.index(e) == 1 and orien == -1):
            turn = -turn
            lengthse.reverse()
        tailshort = (c.index(e) % 2 == 0)
        x0 = crossings[c][0]
        y0 = crossings[c][1]
        im = []
        for l in lengthse:
            if direction == 0:
                x1 = x0 + l
                y1 = y0
            elif direction == 1:
                x1 = x0
                y1 = y0 + l
            elif direction == 2:
                x1 = x0 - l
                y1 = y0
            elif direction == 3:
                x1 = x0
                y1 = y0 - l
            im.append(([[x0, y0], [x1, y1]], l, direction))
            direction = (direction + turn) % 4
            x0 = x1
            y0 = y1
        direction = (direction - turn) % 4
        c2 = [ee for ee in availables if e in ee]
        if len(c2) == 1:
            availables.remove(c2[0])
            crossings[tuple(c2[0])] = (x1, y1, (direction - c2[0].index(e) + 2) % 4)
        c2 = [ee for ee in pd_code if e in ee and ee != list(c)]
        if not c2:
            headshort = not tailshort
        else:
            headshort = (c2[0].index(e) % 2 == 0)
        a = deepcopy(im[0][0])
        b = deepcopy(im[-1][0])
        
        c = 0
        
        p = im[0][0][0]
        if len(im) == 4 and max(x[1] for x in im) == 1:
            # bezier_path([[im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]], **kwargs)
            # bezier_pathの代わりに座標を抽出
            bezier_coords = [im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]
            coords.append(bezier_coords)
            edgeID.append(ce)
            p = im[-1][0][1]
        else:
            while c < len(im)-1:
                if im[c][1] > 1:
                    (a, b) = im[c][0]
                    if b[0] > a[0]:
                        e = [b[0] - 1, b[1]]
                    elif b[0] < a[0]:
                        e = [b[0] + 1, b[1]]
                    elif b[1] > a[1]:
                        e = [b[0], b[1] - 1]
                    elif b[1] < a[1]:
                        e = [b[0], b[1] + 1]
                    # line((p, e), **kwargs)の代わりに座標を抽出
                    coords.append([p, e])
                    edgeID.append(ce)
                    p = e
                if im[c+1][1] == 1 and c < len(im) - 2:
                    xr = round(im[c+2][0][1][0])
                    yr = round(im[c+2][0][1][1])
                    xp = xr - im[c+2][0][1][0]
                    yp = yr - im[c+2][0][1][1]
                    q = [p[0] + im[c+1][0][1][0] - im[c+1][0][0][0] - xp,
                            p[1] + im[c+1][0][1][1] - im[c+1][0][0][1] - yp]
                    # bezier_path([[p, im[c+1][0][0], im[c+1][0][1], q]], **kwargs)の代わりに座標を抽出
                    bezier_coords = [p, im[c+1][0][0], im[c+1][0][1], q]
                    coords.append(bezier_coords)
                    edgeID.append(ce)
                    c += 2
                    p = q
                else:
                    if im[c+1][1] == 1:
                        q = im[c+1][0][1]
                    else:
                        q = [im[c+1][0][0][0] + sign(im[c+1][0][1][0] - im[c+1][0][0][0]),
                                im[c+1][0][0][1] + sign(im[c+1][0][1][1] - im[c+1][0][0][1])]
                    # bezier_path([[p, im[c+1][0][0], q]], **kwargs)の代わりに座標を抽出
                    bezier_coords = [p, im[c+1][0][0], q]
                    coords.append(bezier_coords)
                    edgeID.append(ce)
                    p = q
                    c += 1
        # line([p, im[-1][0][1]], **kwargs)の代わりに座標を抽出
        coords.append([p, im[-1][0][1]])
        edgeID.append(ce)
        # 各セグメントの基本座標も追加
        # for a in im:
        #     coords.append(a[0])
    # print(pd_code)
    # print("crossings",crossings)
    
    return [coords,edgeID,crossings]
def to_float_recursive(obj):
    clsname = type(obj).__name__

    # Python の通常数値
    if isinstance(obj, (int, float)):
        return float(obj)

    # Sage の数値型（文字列で判定）
    elif clsname in (
        "RealDoubleElement", "RealDoubleElement_gsl", "RealLiteral",
        "RealFieldElement", "RDFElement"
    ):
        return float(obj)

    # リスト・タプルは再帰
    elif isinstance(obj, (list, tuple)):
        return [to_float_recursive(x) for x in obj]

    else:
        raise TypeError(f"Cannot convert object of type {type(obj)}")
def apply_x_offset(segment_list, xtra):
    """
    Apply an offset to all x-coordinates in a list of segments.

    Parameters:
        segment_list (list): List of segments, where each segment is a list of points [[x, y], ...].
        xtra (float): Offset to add to all x-coordinates.

    Returns:
        list: Updated list of segments with x-coordinates offset by xtra.
    """
    
    return [[[point[0] + xtra, point[1]] for point in segment] for segment in segment_list]

comp = LINK._isolated_components()
component_gap = 0.5  # Gap between components in the plot
# Handle isolated components individually
coords=[]
edgeID=[]
bending_numbers = BENDING_NUMBERS
if bending_numbers is None:
    bending_numbers = [[] for _ in comp]  # Default to empty bending numbers
crossings=dict()
offset_x=0
for i in range(len(comp)):
    link=Link(comp[i])
    ret = allocate_positions(link,bending_numbers[i])
    xs=[x[0] for x in flatten(ret[0],max_level=1)]   
    coord=apply_x_offset(ret[0], offset_x)
    # As a last resort(苦し紛れ), in the case of abalone,
    # direction of edges are exceptionally inverted
    if len(comp[i])==1:
        coord=[e[::-1] for e in coord]
    offset_x+=max(xs)-min(xs)+ component_gap
    coords.extend(coord)
    edgeID.extend(ret[1])
    crossings.update(ret[2])
#  # Special case for the unknot
for i in range(len(crossings)+1,len(crossings)+1+NCIRCLE):
    coords.append([[[offset_x, 0],[offset_x+1,0], [offset_x+1,1],[offset_x,1]]])
    edgeID.append(i)
    offset_x += 1+component_gap
pd_code = LINK.pd_code()
pos_crossings = [crossings[tuple(x)] for x in pd_code]
ret=[coords, edgeID, pos_crossings]
ret=to_float_recursive(ret)
ret