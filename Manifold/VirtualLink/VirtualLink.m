classdef VirtualLink<handle&matlab.mixin.Copyable
    % VirtualLink diagram data tables (structural representation for knots and links)
    % Coordinate system: +x → right, +y → up
    % vertices ID start from 1 to NV, where NV is the number of vertices
    % edges index start from 1 to 2NV+NC, where NC is the number of components
    % But edges ID may not start from 1.
    % Crossing: both virtual and real crossings, vertex: real crossings only

    properties(SetObservable) % Data tables. As for meaning of each table, see the initialization tables
        CrossingTable
        RCrossingTable
        EdgeTable
        REdgeTable
        StrandTable
        RegionTable
        % listener = addlistener(vl, 'REdgeTable', 'PostSet', @(sec, events) arrayfun(@disp,dbstack)));
    end
    properties
        formatFlag % dictionary of format flags
        isWeighted (1,1)=false
        sageLinked =false
        virtualFlag (1,1)=true
    end
    properties
        name
        description
    end
    properties
        GaussCode (1,:) cell % cell array of Gauss code for each component
        % each cell is a vector of integers, positive for over-crossings, negative for under-crossings

        orientation (1,:) double% +1/-1/0 for right/left/virtual type crossings
        ROrientation (1,:) double % +1/-1 for right/left type real crossings
        PDCode (:,4) double% NV by 4 matrix of PD code
        % [[e1,e2,e3,e4],...] anti-clockwise, starting with the incoming undercrossing.

        Ncircle (1,1) double% the number of circle(unknot) components
        DTCode %not supported
        RGaussCode (1,:) cell % Gauss code with only real crossings
        RPDCode (:,4) double % PD code with only real crossings
        RCrossingMap (1,:) double %
        isCut (1,:) logical  % true if the strand is cut
        % headMap: N×2 int matrix of signed vertex IDs
        % Each row = [under, over] → target crossings for under/over strands
        % +/- indicates the edge enters over/under target
        headMap (:,2) double
        RHeadMap (:,2) double % headMap with only real crossings
    end
    properties (Dependent)
        isVertex (1,:) logical % 1 for virtual crossings, 0 for real crossings
        sageName (1,1) string % variable name in Sage
        sageOrientation (1,:) double % orientation in Sage
    end

    properties(Constant)
        formatFlag0  =dictionary(["G" "PD" "VG" "RPD" "TC" "TRC" "TE" "TRE" "TS"],false(1,9))
        CrossingTableInit=table(...
            nan,missing ,[nan nan], nan,nan(1,4),{nan},nan,nan,  ...
            VariableNames = {
            'ID',           % Unique crossing identifier
            'Label',        % Optional symbolic label
            'HeadMap',      % [under,over] vertex IDs for head map
            'Orientation',  % +1/-1 (r/l type) for crossing; 0 if virtual
            'EdgeID',       % [a,b,c,d] edges in PD code
            'Position',     % x+yi coordinate of crossing center
            'Weight',       % Optional weight for the crossing
            'Region'        % Coarse region/zone for layout
            })

        EdgeTableInit =table(...
            nan,missing, [nan,nan], [false,false],nan, nan, {nan}, {nan}, false, ...
            VariableNames = {
            'ID',           % Unique edge identifier
            'Label',        % Optional edg label
            'Crossing',     % Signed crossing ID [source,target](or NaN for free edg)
            'isOver',       % true if the edge goes over [source,target] crossing
            'Weight',       % Optional weight for the edge
            'StrandID',     % Associated strand ID
            'Position',     % array of complex coordinates for drawing
            'Arc'   ,       % signed crossing IDs for the arc
            'IsVirtual'     % true if arc involves only virtual crossings
            })

        StrandTableInit =table(...
            nan,missing, {nan}, {nan}, {nan},{nan},false,{nan}, ...
            VariableNames = {
            'ID',           % Strand ID (component index)
            'Label',        % Optional strand label (e.g., 'L1')
            'Edges',        % Cell array of edge IDs forming the strand
            'Vertices',     % Cell array of vertex IDs forming the strand
            'Crossings',    % Cell array of crossing IDs forming the strand
            'Weights',      % Optional weight for edges
            'isCut',        % true if the strand is cut
            'Color'         % Optional color
            })

        RegionTableInit =table(...
            nan,missing, {nan}, {nan}, {nan}, ...
            VariableNames = {
            'ID',           % Unique region identifier
            'BoundaryEdges',% List of edge IDs forming the boundary
            'Label',        % Optional region label
            'IsOuter',      % true if it is the unbounded face
            'Zone'          % Coarse position or layout group
            })
        % moveType =["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"]% supported move types
    end

    methods
        function set.REdgeTable(obj,arg)
            obj.REdgeTable=arg;
        end
        %% Validation
        function validateStructure(obj)
        end




        %% Moves
        function convert(obj,type,arg)
            arguments
                obj
                type (1,1) string {mustBeMember(type,["mirror","mirrorMfd","reverse","connectedSum","sum","cyclic"])}  % move type
                arg =[]
            end
            switch type
                case "mirror"
                    error("mirror move is not implemented yet")
                case "mirrorMfd"
                    flag=obj.formatFlag(["G" "VG"]);
                    if flag(1)
                        [gc,ori]=obj.getGaussCode;
                    elseif flag(2)
                        [gc,ori]=obj.getRGaussCode;
                    else
                        error("not implemented")
                    end
                    gc=cellfun(@(x){-x},gc);
                    ori=-ori;
                    if flag(1)
                        obj.setData(Gauss=gc,orientation=ori);
                    elseif flag(2)
                        obj.setData(RGauss=gc,orientation=ori);
                    end
                case "reverse"
                    error("reverse move is not implemented yet")
                case "connectedSum"
                    [gc,ori]=obj.getRGaussCode;
                    assert(size(arg,2)==2&&all(ismember(arg(1,1:2),1:length(gc))), ...
                        "arg must be 2 indices of strand of VirtualLink")
                    gc1=gc{arg(1,1)};
                    gc2=gc{arg(1,2)};
                    if size(arg,1)==2
                        idx=find(arg(2,1)==gc1);
                        gc1=gc1([idx:end,1:idx-1]);
                        idx=find(arg(2,2)==gc2);
                        gc2=gc2([idx:end,1:idx-1]);
                    end
                    NV=length(ori);
                    arr=[-2,-1,1,4,-4,-3,3,2];
                    arr=sign(arr)*NV+arr;
                    gc{arg(1,1)}=[arr(1:4),gc2,arr(5:8),gc1];
                    gc(arg(1,2))=[];
                    ori=[ori,-1,1,-1,1];
                    obj.setData(Gauss=gc,orientation=ori);
                case "sum"
                    assert(isa(arg,"VirtualLink"))
                    [gc1,ori1]=obj.getRGaussCode;
                    [gc2,ori2]=arg.getRGaussCode;
                    NV1=length(ori1);
                    NV2=length(ori2);
                    map=@(x)sign(x)*NV1+x;
                    gc2=cellfun(map,gc2,"UniformOutput",false);
                    gc=[gc1,gc2];
                    ori=[ori1,ori2];
                    obj.setData(Gauss=gc,orientation=ori);
                case "cyclic"
                    % cyclic shift of gauss code
                    % specify argument in the format [strand idx,move num]
                    assert(isnumeric(arg)&&numel(arg)==2)
                    if obj.virtualFlag
                        [gc,ori]=obj.getRGaussCode();
                    else
                        [gc,ori]=obj.getGaussCode();
                    end
                    gc{arg(1)}=circshift(gc{arg(1)},arg(2));
                    if obj.virtualFlag
                        obj.setData(RGaussCode=gc,orientation=ori);
                    else
                        obj.setData(GaussCode=gc,orientation=ori);
                    end
            end
        end
        function convertKnotCompl(obj)
            [vgc,ori]=obj.getRGaussCode;
            NV=length(ori);
            dic = dictionary();
            for i = 1:NV
                dic{i} = 4*(i-1) + [1,2];
                if ori(i)>0
                    dic{-i} = -4*(i-1) - [2,3];
                else
                    dic{-i} = -4*(i-1) - [4,1];
                end
            end
            vgc2=cellfun(@(x) {horzcat(dic{x})},vgc);
            for i = 1:NV
                dic{i} = 4*(i-1) + [3,4];
                if ori(i)>0
                    dic{-i} = -4*(i-1) - [4,1];
                else
                    dic{-i} = -4*(i-1) - [2,3];
                end
            end
            vgc2=[vgc2;cellfun(@(x) {horzcat(dic{fliplr(x)})},vgc)];
            ori2=repmat([-1,1,-1,1],1,NV);
            obj.setData(Gauss=vgc2,orientation=ori2);
        end
        function tbl=movable(obj,type,argstr,arg)
            % MOVABLE Get a table of valid moves for the specified type
            % See `move` for the meaning of arguments and move types.
            % if the table returned is empty, no valid moves are found.
            % tbl: table with columns "param" and "vertex" for the move
            arguments
                obj
                type (1,1) string {mustBeMember(type,["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"])}  % move type
                argstr struct =struct.empty;
                arg.param (1,1) =nan % parameter for the move
                arg.v (1,:) double =[]% vertex indices for move
                arg.e (1,:) double =[]% edge indices for move
                arg.strict (1,1) logical =false % strict mode, fix the order of vertices
            end
            if ~isempty(argstr)
                arg=argstr;
            end
            tbl=table;
            if isfield(arg,"v")
                V=arg.v;
            else
                V=[];
            end
            if isfield(arg,"e"), E=arg.e; end
            strict=arg.strict;
            % if isfield(arg,"e")&&~isempty(arg.e), warning("arg.e is not supported yet"); end
            switch type
                case "R1"
                    assert(length(V)==1||length(arg.e)==1)
                case "R2"
                    assert(length(V)==2||length(arg.e)==2) % 2 vertices or 1 edge for R2 move
                case "R3"
                    assert(length(arg.e)==3||length(V)==3)
                case "MP"
                    assert(any(length(V)==[2,3])) % 2 or 3 vertices for Pachner move
                    V=arg.v;
                    N=length(V);
                    [gc,ori]=obj.getRGaussCode;
                    if N==2
                        T=load("MPmoveData.mat").T_MP_L;
                    else
                        if isequal(V,[0,0,0])
                            arr=1:length(ori);
                            C=binom(arr,3);
                            tbl=table;
                            for j=1:size(C,1)
                                tbl=[tbl;obj.movable("MP",v=C(j,:))];
                            end
                            return
                        end
                        T=load("MPmoveData.mat").T_MP_R;
                    end
                    if ~isa(arg.param,"double")
                        T=T(find(T.param==arg.param),:);
                    end
                    if length(ori)<max(V)
                        return; % no valid moves if vertex IDs exceed
                    end
                    cur=ori(V);
                    X=searchMoves(T,gc,V,cur);
                    % disp(X)
                    if ~isempty(X)
                        tbl=X(:,["param","vertex"]);%消す必要無い？
                    end
                    % tbl.Properties.VariableNames=["name","vertex"];
                    % disp(T)
                case "PS"
                    % 2 vertices or 2 edge for Pure Sliding move (inverse)
                    assert(length(V)==2||length(arg.e)==2)
                    [gc,ori]=obj.getRGaussCode;
                    if length(V)==2
                        assert(V(1)~=V(2)&&V(1)*V(2)>0)
                        [sj,vj]=findC(V(1),gc);
                        cgc=gc{sj}([end,1:end,1]);
                        isadj_p=cgc(vj)==V(2)|| cgc(vj+2)==V(2);
                        [sj,vj]=findC(-V(1),gc);
                        cgc=gc{sj}([end,1:end,1]);
                        isadj_n=cgc(vj)==-V(2)|| cgc(vj+2)==-V(2);
                        if isadj_p&&isadj_n
                            tbl=table("PS",V,VariableNames=["param","vertex"]);
                        end
                    else
                        error("PS move is not implemented yet")
                    end
                    N=length(V);
                    load("PSmoveData.mat","T_PS_R");
                    cur=ori(V);
                    X=searchMoves(T_PS_R,gc,V,cur);
                    if ~isempty(X)
                        tbl=X(:,["param","vertex"]);
                    end
                case "CP"
                    assert(any(length(V)==[3,5]))
                    load("CPmoveData.mat","T_CP_R","T_CP_L");
                    [gc,ori]=obj.getRGaussCode;
                    if length(V)==3
                        X=searchMoves(T_CP_L,gc,V,ori(V));
                    else
                        X=searchMoves(T_CP_R,gc,V,ori(V));
                    end
                    if ~isempty(X)
                        tbl=X(:,["param","vertex"]);
                    end
                case "02" % move 0-2
                    assert(length(V)==2||length(arg.e)==1) % 2 vertices or 1 edge for 0-2 move
                    load("M02moveData.mat","T_02_R");
                    % inv
                    assert(ismember(arg.param,["inv","dir"]),'parameter must be "inv" or "dir"')
                    if arg.param=="inv" % number of vertices: 2 to 0
                        T=T_02_R;
                        [gc,ori]=obj.getRGaussCode;
                        cur=ori(V);
                        X=searchMoves(T,gc,V,cur);
                        if ~isempty(X)
                            tbl=X(:,["param","vertex"]);
                        end
                    else % dir: number of vertices: 0 to 2
                        gc=obj.getRGaussCode;
                        [sj1,vj1]=findC(V(1),gc);
                        [sj2,vj2]=findC(V(2),gc);
                        if((sj1==sj2)&&(vj1+1==vj2||vj1==length(gc{sj1})&&vj2==1))
                            % the vertices must be adjacent
                            tbl=table("dir",arg.v,VariableNames=["param","vertex"]);
                        end
                    end
                case "H"
                    assert(length(V)==1) % 2 vertices or 1 edge for handle move
                    gc=obj.getRGaussCode;
                    NV=max(horzcat(gc{:}));
                    assert(V(1)<=NV,"vertex index exceeds the number of vertices")
                    tbl=table(0,V,VariableNames=["param","vertex"]);
                case "BMP" % Bumping move
                    %: 4通りのBMP moveを網羅できていない
                    assert(length(V)==2||length(V)==3) % 2 or 3 vertices for Bumping move
                    if length(V)==2
                        assert((V(1)~=V(2))&&(V(1)*V(2)>0))
                        isPositiveV=V(1)>0;
                        V=abs(V(:).');
                        [gc,ori]=obj.getRGaussCode;
                        [sj,vj]=arrayfun(@(x)findC(x,gc),[V,-V]);
                        acceptable=[
                            sj(1)==sj(2)&&vj(1)+1==vj(2)&&all(ori(V)==[-1,1]);
                            sj(3)==sj(4)&&vj(3)+1==vj(4)&&all(ori(V)==[1,-1]);
                            sj(2)==sj(1)&&vj(2)+1==vj(1)&&all(ori(V)==[1,-1]);
                            sj(4)==sj(3)&&vj(4)+1==vj(3)&&all(ori(V)==[-1,1])];
                        tbl=table(strings(4,1),[V;-V;V([2,1]);-V([2,1])],VariableNames=["param","vertex"]);
                        if strict, acceptable=acceptable(2-isPositiveV); end
                        tbl=tbl(acceptable,:);
                    else
                        [gc,ori]=obj.getRGaussCode();
                        cur=ori(V);
                        load("moveData/BMPmoveData.mat","T_BMP_R")
                        X=searchMoves(T_BMP_R,gc,V,cur);
                        if ~isempty(X)
                            tbl=X(:,["param","vertex"]);
                        end
                    end
                case "B02" % Bumping 0-2 move
                    assert(length(V)==2||isequal(size(E),[2,2]),"invalid parameter")
                    if length(V)==2
                        assert(V(1)~=V(2))
                        [gc,ori]=obj.getRGaussCode;
                        sgn=ori(V);
                        if sgn(1)*sgn(2)~=-1 % opposite orientation
                            return; % no valid moves
                        elseif sgn(1)==-1
                            V=V([2,1]); % make sure V(1) is 1
                        end
                        % to check if a strand made of -V is included
                        % if any(cellfun(@(x)isequal(x,-V),gc))
                        if any(cellfun(@(x)isequal(x,-V)||isequal(x,V),gc))
                            tbl=table("",V,VariableNames=["param","vertex"]);
                            % elseif any(cellfun(@(x)isequal(x,-V([2,1])),gc))
                        elseif any(cellfun(@(x)isequal(x,-V([2,1]))||isequal(x,V([2,1])),gc))
                            tbl=table("",V([2,1]),VariableNames=["param","vertex"]);
                        end
                    else
                        % not implemented
                        tbl=table("",E(:)',VariableNames=["param","edge"]);
                    end
                otherwise
            end
            % issue: 同一データでも適用位置によって複数解が出てくる可能性を考慮する必要
            tbl=unique(tbl,'rows');
            if nargout==0
                disp(tbl)
            end
            function ret=flatGC(gc,n)
                % Flatten the Gauss code
                Ns=length(gc);
                for ii=1:Ns
                    gc{ii}=[gc{ii},gc{ii}(mod(0:n-1,Ns)+1),0];
                end
                ret=[horzcat(gc{:}),zeros(1,n)];
            end
            function X=searchMoves(T,gc,V,ori_c)
                % Search the permutations of vertices for all possible moves
                % dic: dictionary T.gc->gc
                % X:solution table
                % ori_c: orientation of vertex used to replace
                % V: vertex index corresponds to 1:N
                N_=length(V);

                dic=dictionary([(1:N_),-N_:-1],1:2*N_);
                X=table;
                if strict
                    P=1:N_;
                else
                    P=perms(1:N_); % consider all permutations of vertices V
                end
                for iP=1:size(P,1)
                    % iP: current permutation index of V
                    % Tj: move type index with matching orientation
                    Tj=find(all(ori_c(P(iP,:))==T.ori,2));
                    vertex=V(P(iP,:));
                    dic(1:N_)=vertex;
                    dic(-1:-1:-N_)=-vertex;
                    for iTj=1:size(Tj,1)
                        % iTj: current move type index
                        [sIdx,vIdx]=arrayfun(@(x)findC(x,gc),dic(T.gcFirst(Tj(iTj),:)));
                        gc_r=T.gc(Tj(iTj),:);
                        NS=size(T.gc,2);
                        flag=false(1,NS);
                        for iS = 1:NS
                            % Check if the fragment of Gauss code matches the segment
                            vIdxMod = mod(vIdx(iS) - 1:vIdx(iS) + length(gc_r{iS}) - 2, length(gc{sIdx(iS)})) + 1;
                            flag(iS) = isequal(dic(gc_r{iS}), gc{sIdx(iS)}(vIdxMod));
                        end
                        if all(flag)
                            % If all fragments match, add to the result
                            X=[X;[table(vertex),T(Tj(iTj),:)]];
                        end
                    end
                end
            end
        end
        function move(obj,type,argstr,arg)
            % MOVE performs a move on the virtual link diagram
            % arguments:
            % type: move type, one of ["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"]
            % argstr: structure with move parameters
            % arg: move parameters
            % Move list:
            % - R1, R2, R3: Reidemeister moves
            % - MP: Pachner 2-3 move
            % - PS: Pure Sliding move
            % - CP: Combinatorial Pontryagin move
            % - 02: 0-2 move
            % - H: H move
            % - BMP: Bumping MP move
            % - B02: Bumping 0-2 move
            % parameters:
            % arg.v: vertex indices for the move
            % arg.e: edge indices for the move
            % arg.param: parameter for the move
            % for MP move, specify A1~D4 types of MP moves
            % for BMP move, specify param=0 or param=1 for the orientation of the new circle
            arguments
                obj
                type (1,1) string {mustBeMember(type,["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"])}  % move type
                argstr struct =struct.empty;
                arg.param (1,1)=nan % parameter for the move
                arg.v (1,:) double % vertex indices for move
                arg.e double % edge indices for move
            end
            if ~isempty(argstr)
                arg=argstr;
            end
            arg.strict=true;
            tbl=obj.movable(type,arg);
            if height(tbl)==0
                error("No valid moves found.");
            elseif height(tbl)>1
                warning("moves may be ambiguous. check the specification");
            end
            if isfield(arg,"v"), V=arg.v; end
            if isfield(arg,"e"), E=arg.e; end

            switch type
                case "R1"
                    error("R1 move is not implemented yet")
                case "R2"
                    error("R2 move is not implemented yet")
                case "R3"
                    error("R3 move is not implemented yet")
                case "MP"
                    N=length(V);
                    [gc,ori]=obj.getRGaussCode;
                    load("MPmoveData.mat","T_MP_L","T_MP_R");
                    param=string(arg.param);
                    if N==2
                        idx=find(param==T_MP_L.param,1);
                        T0=T_MP_L(idx,:);
                        T1=T_MP_R(idx,:);
                    else
                        idx=find(param==T_MP_R.param,1);
                        T0=T_MP_R(idx,:);
                        T1=T_MP_L(idx,:);
                    end

                    [gc,ori,Vnew,dic2]=convertGaussCode(gc,ori,T0,T1,V);
                    TE0=obj.REdgeTable;
                    obj.setData(RGauss=gc,orientation=ori);

                    if obj.isWeighted
                        convertWeight(obj,V,TE0,T0,T1)
                    end

                case "PS"
                    if ~isempty(V)
                        [gc,ori]=obj.getRGaussCode;
                        idx=double(tbl.param(1));
                        load("PSmoveData.mat","T_PS_R");
                        T0=T_PS_R(idx,:);
                        load("PSmoveData.mat","T_PS_L");
                        [gc,ori]=convertGaussCode_(gc,ori,T0,T_PS_L,V);
                        obj.setData(RGauss=gc,orientation=ori);
                    else

                        error("PS move is not implemented yet")
                    end


                case "CP"
                    load("CPmoveData.mat","T_CP_R","T_CP_L");
                    [gc,ori]=obj.getRGaussCode;
                    if length(V)==3
                        T0=T_CP_L; T1=T_CP_R;
                    else
                        T0=T_CP_R; T1=T_CP_L;
                    end
                    [gc,ori]=convertGaussCode_(gc,ori,T0,T1,V);
                    obj.setData(RGauss=gc,orientation=ori);
                    if obj.isWeighted
                        error("weighted CP move is not implemented yet")
                    end
                case "02"
                    assert(length(V)==2||length(arg.e)==1) % 2 vertices or 1 edge for 0-2 move
                    load("M02moveData.mat","T_02_R","T_02_L");
                    if arg.param=="inv" % number of vertices: 2 to 0
                        [gc,ori]=obj.getRGaussCode;
                        [gc,ori]=convertGaussCode_(gc,ori,T_02_R,T_02_L,V);
                        TE0=obj.REdgeTable;
                        obj.setData(RGauss=gc,orientation=ori);
                        if obj.isWeighted
                            convertWeight(obj,V,TE0,T_02_R,T_02_L)
                        end
                    else % dir: number of vertices: 0 to 2
                        [gc,ori]=obj.getRGaussCode;
                        [sj,vj]=findC(V(1),gc);
                        NV=max(horzcat(gc{:}));
                        dic=dictionary([-2,-1,0,1,2],[-NV-2,-NV-1,0,NV+1,NV+2]);
                        gc{sj}=[gc{sj}(1:vj),dic(T_02_R.gc{1}),gc{sj}(vj+1:end)];
                        ori=[ori,T_02_R.ori];
                        obj.setData(RGauss=gc,orientation=ori);
                        if obj.isWeighted
                            TE0=obj.REdgeTable;
                            [~,widx]=ismember(V,TE0.Crossing,"rows");
                            wdata=T_02_R.weight{1};
                            TE10_=repmat(obj.EdgeTableInit,5,1);
                            TE10_.Crossing=dic(wdata(:,1:2));
                            TE10_.Weight=wdata(:,3);
                            TE10_.Weight(1)=TE0.Weight(widx);
                            TE10_.Crossing(1,1)=V(1);
                            TE10_.Crossing(2,2)=V(2);
                            TE0(widx,:)=[];
                            obj.REdgeTable=[TE0;TE10_];
                            obj.formatFlag("TRE")=true;
                        end
                    end
                case "H"
                    assert(~isempty(obj.movable("H",arg)),"invalid arguments")
                    assert(obj.isWeighted&&~isnan(arg.param),"must be weighted")
                    gc=obj.getRGaussCode;
                    [sj,vj]=findC(V,gc);
                    cgc=[gc{sj}(end),gc{sj},gc{sj}(1)];
                    edge=cgc(vj+[0,1;1,2]);
                    [sj,vj]=findC(-V,gc);
                    cgc=[gc{sj}(end),gc{sj},gc{sj}(1)];
                    edge=[edge;cgc(vj+[0,1;1,2])];
                    [~,idx]=ismember(edge,obj.REdgeTable.Crossing,"rows");
                    warray=[-1,1,-1,1];
                    for j=1:length(idx)
                        obj.REdgeTable.Weight(idx(j),:)=obj.REdgeTable.Weight(idx(j),:)+ ...
                            warray(j)*arg.param;
                    end
                case "BMP"
                    if length(V)==2
                        if ~ismember(arg.param,[0,1])
                            error("specify the parameter param=0,or param=1")
                        else
                            orientationOfNewCircle=logical(double(arg.param));
                        end
                        % issue: mirrored version
                        tbl=obj.movable("BMP",v=V);
                        assert(height(tbl)==1 && isequal(tbl.vertex, V), "BMP moves are not unique or invalid.");
                        [gc,ori]=obj.getRGaussCode;
                        S=sign(V(1));
                        [sj1,vj1]=findC(-V(1),gc);
                        gc{sj1}=gc{sj1}([vj1+1:end,1:vj1-1]);
                        [sj2,~]=findC(-V(2),gc);
                        [sj0,vj0]=findC(V(1),gc);
                        Vnew=S*(max([gc{:}])+1);
                        gc{sj0}(vj0)=Vnew;
                        vj0=find(gc{sj0}==V(2),1);
                        gc{sj0}(vj0)=[];
                        vj2=find(gc{sj2}==-V(2),1);
                        if orientationOfNewCircle
                            gc_circle=[-V(1:2),-Vnew];
                            ori(abs([V,Vnew]))=-[-1,1,1]*S;
                        else
                            gc_circle=[-V([2,1]),-Vnew];
                            ori(abs([V,Vnew]))=[-1,1,1]*S;
                        end
                        if sj1==sj2 % if the vertices V are in the same strand, split into two strands
                            gc(length(gc)+(1:2))={[V(1),gc{sj2}(1:vj2-1)],[V(2),gc{sj2}(vj2+1:end)]};
                            gc{sj2}=gc_circle;
                        else % if the vertices are in different strands, form one strand
                            gc{sj2}=[V(2),gc{sj2}([vj2+1:end,1:vj2-1]),V(1),gc{sj1}];
                            gc{sj1}=gc_circle;
                        end
                        obj.setData(RGauss=gc,orientation=ori);
                    else
                        assert(ismember(arg.param,["E1";"F1";"E2";"F2"]))
                        [gc,ori]=obj.getRGaussCode;
                        load("moveData/BMPmoveData.mat","T_BMP_L","T_BMP_R")
                        T_BMP_L=T_BMP_L(arg.param==T_BMP_L.param,:);
                        T_BMP_R=T_BMP_R(arg.param==T_BMP_R.param,:);
                        [gc,ori]=convertGaussCode(gc,ori,T_BMP_R,T_BMP_L,arg.v);
                        obj.setData(RGauss=gc,orientation=ori);
                    end
                case "B02"
                    if isfield(arg,"e")&&isequal(size(arg.e),[2,2])
                        assert(ismember(arg.param,[1,2]),"specify the parameter param=1,or param=2")
                        [gc,ori]=obj.getRGaussCode;
                        load("moveData/B02moveData.mat","T_B02_L","T_B02_R")
                        T_B02_R=T_B02_R(double(arg.param)==double(T_B02_R.param),:);% orientationOfNewCircle
                        [gc,ori]=convertGaussCodeFromEdge(gc,ori,T_B02_L,T_B02_R,arg.e);
                        obj.setData(RGauss=gc,orientation=ori);
                        % [sj,vj]=arrayfun(@(x)findC(x,gc),E(:,1));
                        % error("not implemented")
                    elseif length(V)==2
                        assert(all(V>0))
                        tbl=obj.movable("B02",v=V);
                        assert(height(tbl)==1,"No valid BMP move found");
                        [gc,ori]=obj.getRGaussCode;
                        NV=max(cellfun(@(x)max(x),gc));
                        dic=dictionary(setdiff(1:NV,V),1:NV-2);
                        dic(-dic.keys)=-dic.values;
                        dic(0)=0;
                        isVpositive=any(cellfun(@(x)isequal(x,V)||isequal(x,V([2,1])),gc));
                        if isVpositive, V=-V; end
                        [sj0,~]=findC(-V(1),gc);
                        [sj2,vj2]=findC(V(2),gc);
                        gc{sj2}=gc{sj2}([vj2+1:end,1:vj2-1]);
                        [sj1,vj1]=findC(V(1),gc);
                        if sj1==sj2 % if the vertices V are in the same strand, split into two strands
                            gc([sj0,sj1])={gc{sj1}(1:vj1-1),gc{sj1}(vj1+1:end)};
                        else % if the vertices are in different strands, form one strand
                            gc{sj0}=[gc{sj2},gc{sj1}([vj1+1:end,1:vj1-1])];
                            gc([sj1,sj2])=[];
                        end
                        gc=cellfun(@(x){dic(x)},gc);
                        ori(abs(V))=[];
                        obj.setData(RGauss=gc,orientation=ori);
                    else
                        error("invalid parameter or B02 move is not implemented yet")
                    end
            end
            function [gc,ori,Vnew,map]=convertGaussCode_(gc,ori,T0,T1,V)
                % Convert Gauss code based on the move(legacy)
                % T0: original Gauss code data
                % T1: new Gauss code data
                % V: vertices involved in the move
                NV_=max([gc{:}]);
                map=struct;
                if T0.NV<=T1.NV
                    % if the number of vertices increases after the move
                    V(T0.NV+1:T1.NV)=NV_+(1:T1.NV-T0.NV);
                    Mgc01=@(x)x;
                    ori(V)=T1.ori;
                    Vnew=V(1:T1.NV);
                elseif T0.NV>T1.NV
                    Vsort=sort(V);
                    vrem=setdiff(1:NV_,Vsort(end-(0:T0.NV-T1.NV-1)));
                    Mgc01=dictionary(vrem,1:NV_-(T0.NV-T1.NV));
                    Mgc01(Vsort(1:T1.NV))=Vsort(1:T1.NV);
                    Mgc01(-Mgc01.keys)=-Mgc01.values;
                    Mgc01(0)=0;
                    ori(Vsort(1:T1.NV))=T1.ori;
                    ori=ori(vrem);
                    Vnew=Vsort(1:T1.NV);
                else
                    Mgc01=@(x)x;
                    ori(V)=T0.ori;
                    Vnew=V;
                end
                % issue: don't we replace T1.gc by Mt2gc1?and Mtgc0?
                Mt2gc1=dictionary([1:length(V),-1:-1:-length(V)],[V,-V]);
                % position of start point of each fragment
                [startStrandIdx,startVertexIdx]=arrayfun(@(x)findC(x,gc),Mt2gc1(T0.gcFirst));

                for s=1:length(gc)
                    % s:current component
                    % update Gauss code for each component
                    k=find(startStrandIdx==s);
                    if isempty(k)
                        gc{s}=Mgc01(gc{s}); % no related crossings in this component
                        continue; % no related crossings in this component
                    end
                    [~,sortidx]=sort(startVertexIdx(k));
                    k=k(sortidx);
                    % gcpart: temporary Gauss code [gc,T1,gc,T1,...,gc] T1 start at vj_(k)
                    % dic1: dictionary to convert T1.gc to V
                    % gc01: delete and shift indices of crossings
                    gcpart=gc{s}(1:startVertexIdx(k(1))-1);
                    for i=1:length(k)-1
                        gcpart=[gcpart,Mt2gc1(T1.gc{k(i)}), ...
                            gc{s}(startVertexIdx(k(i))+length(T0.gc{k(i)}):startVertexIdx(k(i+1))-1)];
                    end
                    vjend=startVertexIdx(k(end))+length(T0.gc{k(end)});
                    gcpart=[gcpart,Mt2gc1(T1.gc{k(end)}),gc{s}(vjend:end)];
                    if vjend>length(gc{s})+1
                        % if at the end of the Gauss code, restart from the beginning
                        gcpart(1:(vjend-length(gc{s})-1))=[];
                    end
                    gc{s}=Mgc01(gcpart);
                end
                if length(T1.gc)>length(T0.gc)
                    % if the number of strands increases after the move, add new strands
                    gc=[gc,cellfun(@(x){Mt2gc1(x)},T1.gc(length(T0.gc)+1:end))];
                elseif length(T1.gc)<length(T0.gc)
                    % if the number of strands decreases after the move, remove strands
                    sj_remove=startStrandIdx(length(T1.gc)+1:end);
                    % strand to remove must be empty
                    assert(all(cellfun(@(x)isempty(x),gc(sj_remove))))
                    gc(sj_remove)=[]; % remove empty strands
                end
            end
            function [gc,ori,Vnew,map]=convertGaussCodeFromEdge(gc,ori,T0,T1,E)
                % Convert Gauss code based on the move 
                % brash up to handle bumping moves but not completed for
                % all moves
                % T0: original Gauss code data
                % T1: new Gauss code data
                % V: vertices involved in the move
                V_=length(ori)+1;
                T0.gc={1,-1};
                T0.gcFirst=[1,-1];
                T0.NV=1;
                after=E(:,1)~=0;
                for i=1:2
                    [sj_,vj_]=findC(E(i,2-after(i)),gc);
                    % insert ±V at the edges specified
                    gc{sj_}=gc{sj_}([1:vj_,vj_:end]);
                    gc{sj_}(vj_+after(i))=V_*(-1)^(i-1);
                end
                ori(end+1)=0;
                 [gc,ori,Vnew,map]=convertGaussCode(gc,ori,T0,T1,V_);
            end
            function [gc,ori,Vnew,map]=convertGaussCode(gc,ori,T0,T1,V)
                % 
                % issue: when seed data contain empty array?
                NV_=length(ori);
                map=struct;
                if T0.NV<=T1.NV
                    % if the number of vertices increases after the move
                    V(T0.NV+1:T1.NV)=NV_+(1:T1.NV-T0.NV);
                    Mgc01=@(x)x;
                    ori(V)=T1.ori;
                    Vnew=V(1:T1.NV);
                elseif T0.NV>T1.NV
                    Vsort=sort(V);
                    vrem=setdiff(1:NV_,Vsort(end-(0:T0.NV-T1.NV-1)));
                    Mgc01=dictionary(vrem,1:NV_-(T0.NV-T1.NV));
                    Mgc01(Vsort(1:T1.NV))=Vsort(1:T1.NV);
                    Mgc01(-Mgc01.keys)=-Mgc01.values;
                    Mgc01(0)=0;
                    ori(Vsort(1:T1.NV))=T1.ori;
                    ori=ori(vrem);
                    Vnew=Vsort(1:T1.NV);
                end
                Mt2gc0=dictionary([1;-1]*(1:length(V)),[1;-1]*V);
                Mt2gc1=dictionary([1;-1]*(1:length(Vnew)),[1;-1]*Vnew);
                % position of start point of each fragment
                [startStrandIdx,startVertexIdx]=arrayfun(@(x)findC(x,gc),Mt2gc0(T0.gcFirst));
                cap=zeros(2,0);
                frag={};
                if length(T1.gc)<length(T0.gc)
                    % if the number of strands decreases after the move, remove strands
                    sj_remove=startStrandIdx(length(T1.gc)+1:end);
                    gc_remove=gc(sj_remove);
                    gc_remove2=cellfun(@(x){Mt2gc0(x)},T0.gc(length(T1.gc)+1:length(T0.gc)));
                    cond=cellfun(@(A,B)length(A)==length(B) && contains(join(string([B B])), join(string(A))) ...
                        ,gc_remove,gc_remove2);
                    % strand to remove must be empty
                    assert(all(cond),"cannot apply move")
                else
                    sj_remove=[];
                end
                Ncap=max([T0.cap,T1.cap]);
                % divide Gauss code into fragment
                for si=1:length(gc)
                    % s:current component
                    % k:current strand to devide into fragment
                    % update Gauss code for each component
                    si0s=find(startStrandIdx==si);
                    
                    if ismember(si,sj_remove)
                        continue; % remove strand to delete
                    elseif isempty(si0s)% no related crossings in this component
                        cap(:,end+1)=[1;1]*(Ncap+si);
                        frag{end+1}=Mgc01(gc{si});
                        continue; % no related crossings in this component
                    end
                    iscut=gc{si}(end)==0;
                    [~,sortidx]=sort(startVertexIdx(si0s));
                    si0s=si0s(sortidx);
                    delimeter=[startVertexIdx(si0s),length(gc{si})+1];
                    % Since the Gauss code is cyclic, if a change wraps from the end to the start,
                    % remove the first Ndel elements to make it work correctly.
                    Ndel=delimeter(end-1)+length(T0.gc{si0s(end)})-(length(gc{si})+1);
                    frag{end+1}=Mgc01(gc{si}(max(Ndel,0)+1:delimeter(1)-1));
                    % gcpart: temporary Gauss code [gc,T1,gc,T1,...,gc] T1 start at vj_(k)
                    % dic1: dictionary to convert T1.gc to V
                    % gc01: delete and shift indices of crossings
                    Nsi0s=length(si0s);
                    fi=length(frag);
                    frag(end+(1:2*Nsi0s))={[]};
                    cap(1:2,end+(1:2*Nsi0s+1))=(Ncap+si)*(~iscut); % when not cut, connect fragment at the end
                    for i=1:Nsi0s
                        si0=si0s(i);
                        frag{fi+2*i-1}=Mt2gc1(T1.gc{si0});
                        frag{fi+2*i}=Mgc01(gc{si}(delimeter(i)+length(T0.gc{si0}):delimeter(i+1)-1));
                        cap(2*fi+4*(i-1)+(0:3))=[-si0,-si0,T1.cap(si0),T0.cap(si0)];
                    end
                end
                % table(cap',frag')
                gc={};
                f_used=false(1,length(frag));
                % merge fragments to form Gauss code
                while ~all(f_used)
                    cur=find(~f_used,1);
                    next=cur;
                    f_used(next)=true;
                    gc{end+1}=frag{cur};
                    while cap(1,next)~=0
                        next=cap(1,next)==cap(2,:);
                        if f_used(next), break; end
                        f_used(next)=true;
                        gc{end}=[frag{next},gc{end}];
                    end
                    prev=cur;
                    while cap(2,prev)~=0
                        prev=cap(2,prev)==cap(1,:);
                        if f_used(prev), break; end
                        f_used(prev)=true;
                        gc{end}=[gc{end},frag{prev}];
                    end  
                end    
                % gc=cellfun(@(x){Mgc01(x)},gc);
                if length(T1.gc)>length(T0.gc)
                    % if the number of strands increases after the move, add new strands
                    gc=[gc,cellfun(@(x){Mt2gc0(x)},T1.gc(length(T0.gc)+1:end))];
                end
            end
            function convertWeight(obj,V,TE0,T0,T1)
                % update the weight of the crossing
                % 現状，MP moveなどに対してadhocな実装にしている
                Mt2gc0=dictionary(-N:N,[-fliplr(V),0,V]);
                wdata0=T0.weight{1};
                wdata1=T1.weight{1};
                assert(~isequal(wdata0,nan),"not supported yet");

                jmiddle=~(wdata0(:,4));
                jE=zeros(size(wdata0,1),1);
                for i=1:size(wdata0,1) % find the edge indices from data(:,1:2)
                    iodir=wdata0(i,4); %=0,1,2
                    if iodir==0
                        jE(i)=find(TE0.Crossing(:,1)==Mt2gc0(wdata0(i,1)),1);
                    else
                        jE(i)=find(TE0.Crossing(:,iodir)==Mt2gc0(wdata0(i,iodir)),1);
                        vother=TE0.Crossing(jE(i),3-iodir);
                        % wdata0(i,3-iodir)=vother;
                        if ismember(vother,[V,-V])
                            wdata1(i,3-iodir)=wdata1(find(Mt2gc0(wdata0(:,3-iodir))==vother),3-iodir);
                        else
                            wdata1(i,3-iodir)=vother;
                        end
                    end
                end
                % issue: adjust weight by applying H-move
                wmatch=TE0.Weight(jE(jmiddle))==wdata0(jmiddle,3);
                if ~all(wmatch)
                    Idx=find(jmiddle);
                    Idx=Idx(~wmatch);
                    str=sprintf(join(join(string(wdata0(Idx,1:2)),",",2),"]:%d ,[")+"]:%d",wdata0(Idx,3).');
                    error("weight mismatch:\n weight must be %s",str);
                end
                CrossingArr=[];
                WeightArr=[];
                jEuA=unique(jE(~jmiddle).');
                for jEu=jEuA
                    WeightArr(end+1)=sum([-wdata0(jEu==jE,3);wdata1(jEu==jE,3);TE0.Weight(jEu)]);
                    widx=find(jEu==jE,1);
                    CrossingArr(end+1,:)=wdata1(widx,1:2);
                end
                jm2=~wdata1(:,4);
                WeightArr=[WeightArr.';wdata1(jm2,3)];
                CrossingArr=[CrossingArr;wdata1(jm2,1:2)];
                TE10=repmat(obj.EdgeTableInit,length(jEuA)+sum(jm2),1);
                TE10.Weight=WeightArr(:);
                TE10.Crossing=CrossingArr;
                TE0(jE,:)=[];
                obj.REdgeTable=[TE0;TE10];
                obj.formatFlag("TRE")=true;
            end
        end

    end
    methods(Static)
        function [gc,ori]=calcReverseStrand(gc,ori,reverseFlag)
            % Calculate the reverse strand Gauss code and orientation
            % gc: cell array of Gauss code for each component
            % ori: orientation vector
            v=[];
            assert(length(gc)==length(reverseFlag))
            for i=1:length(gc)
                if reverseFlag(i)
                    gc{i}=fliplr(gc{i});
                    v=[v,gc{i}];
                end
            end
            vp=v(v>0);
            vpm=vp(ismember(-vp,v));
            v=setdiff(abs(v),vpm);
            v0s=v(ori(v)==0);
            ori(v)=-ori(v);
            for v0=v0s(:)'
                % set orientation of virtual crossings to 0
                [idxC,idx]=findC(v0,gc);
                gc{idxC}(idx)=-v0;
                [idxC,idx]=findC(-v0,gc);
                gc{idxC}(idx)=v0;
            end
        end
    end
    %% fundamental calculations
    methods

        function buildFromGaussCode(obj)
        end

        function buildFromPdCode(obj)
        end


        function updateRegions(obj)
        end

        function ret=get.sageName(obj)
            ret="Link"+string(mod(keyHash(obj),1000));
        end
        function ret=get.isVertex(obj)
            % Return a logical array indicating virtual crossings
            ret=~logical(obj.orientation);
        end
        function set.isVertex(obj,arg)
            % isVertex is read-only property
            obj.orientation(logical(arg))=0; % set orientation to 0 for virtual crossings
        end
        function [ret,criterion]=isClosed(obj)
            % isClosed : determine obj represents closed normal o-graph
            % C1 頂点を取り除いて, 対応する辺同士をつなげるとただ 1 つの単純閉曲線となる.
            % C2 [BP] の Figure 1.2 で 3 価グラフに置き換えたものは連結である.
            % C3 [BP] の Figure 1.3 で置き換えてできる単純閉曲線を考える. このとき連結成分の個数は元の
            % グラフの頂点 +1 と等しい.

            criterion=false(1,3);
            [gc,ori]=obj.getRGaussCode;
            criterion(1)=length(gc)==1; % C1

            G=obj.getGraphA;
            criterion(2)=max(G.conncomp)==1; % C2
            % issue: calc vg2g isn't needed
            tbl=obj.getDiskTable;
            NV=length(ori);
            criterion(3)=height(tbl)==NV+1; %C3

            ret=all(criterion);
        end
    end
    %% setter/getter
    % いちいち説明をつけない
    methods
        function set.GaussCode(obj,arg)
            obj.StrandTable.Crossings=arg(:);
        end
        function ret=get.GaussCode(obj)
            ret=obj.StrandTable.Crossings;
        end
        function set.PDCode(obj,arg)
            obj.CrossingTable.EdgeID=arg;
        end
        function ret=get.PDCode(obj)
            ret=obj.CrossingTable.EdgeID;
        end
        function set.RGaussCode(obj,arg)
            obj.StrandTable.Vertices=arg(:);
        end
        function ret=get.RGaussCode(obj)
            ret=obj.StrandTable.Vertices;
        end
        function set.RPDCode(obj,arg)
            obj.RCrossingTable.EdgeID=arg;
        end
        function ret=get.RPDCode(obj)
            ret=obj.RCrossingTable.EdgeID;
        end
        function set.orientation(obj,arg)
            obj.CrossingTable.Orientation=arg(:);
        end
        function ret=get.orientation(obj)
            ret=obj.CrossingTable.Orientation';
        end
        function set.ROrientation(obj,arg)
            obj.RCrossingTable.Orientation=arg(:);
        end
        function ret=get.ROrientation(obj)
            ret=obj.RCrossingTable.Orientation';
        end
        function set.isCut(obj,arg)
            obj.StrandTable.isCut(:)=arg;
        end
        function ret=get.isCut(obj)
            ret=obj.StrandTable.isCut;
        end

    end
    %% system
    methods
        function ret=get.sageOrientation(obj)
            % Return orientation in Sage format
            % convert 0 to 1 for Sage compatibility
            ret=obj.orientation;
            ret(ret==0)=+1;
        end
        function set.sageOrientation(obj,arg)
            % Set orientation from Sage format
            jreal=~obj.isVertex;
            obj.orientation(jreal)=arg(jreal);
        end
        function obj=VirtualLink()
            obj.formatFlag=obj.formatFlag0;

        end
        function delete(obj)
            if obj.sageLinked
                SageWrapper.H.exec(sprintf("del %s",obj.sageName));
            end
        end

        function setSageLink(obj)
            SW=SageWrapper.H;
            Pname=obj.sageName;
            if obj.formatFlag('G')
                data=sprintf("[%s,%s]", ...
                    SW.toStr(obj.GaussCode), ...
                    SW.toStr(obj.sageOrientation,1));
            elseif obj.formatFlag('PD')
                data=SW.toStr(obj.PDCode,2);
                if ismissing(data)
                    data="[]";
                end
            else
                error('Invalid format flag. Please check the input data.');
            end
            SW.exec(sprintf("%s=Link(%s)",Pname,data));
            obj.sageLinked=true;
        end
    end
    %% graphics
    methods
        function [PE,V,E]=calcPositions(obj,strct)
            % calcPositions determines the positions of vertices, edges
            % issue: Ncircle are not implemented
            % issue: plot fails when bending number is not appropriate
            %  (also in case of auto layouting).
            %
            SW=SageWrapper.H;
            try
                SW.exec(obj.sageName)
            catch ME
                obj.setSageLink; % if not linked, link to Sage
            end
            gap=strct.componentGap;
            bendingNumbers=strct.bendingNumbers;
            if isempty(bendingNumbers)
                bendingNumbers="None";
            else
                bendingNumbers=SW.toStr(bendingNumbers);
            end
            Pname=obj.sageName;
            cmd=readlines(fullfile(TopologyConfig.ProjectPath, ...
                'Manifold','VirtualLink','allocate_pos.py'));
            cmd=strrep(cmd,"LINK",Pname);
            cmd=strrep(cmd,"BENDING_NUMBERS",bendingNumbers);
            cmd=strrep(cmd,"COMPONENT_GAP",string(gap));
            cmd=strrep(cmd,"NCIRCLE",string(obj.Ncircle));
            try
                C=SW.exec(cmd);
            catch ME
                % if contains(ME.message, 'MIPSolverException') && ...
                %    contains(ME.message, 'proven infeasible')
                %     error("plot algorithm error")
                %     % 必要ならばリカバリ処理をここに書く
                % else
                %     rethrow(ME);  % 他のエラーは通常通り投げる
                % end
                rethrow(ME);
            end
            pxy=@(x)complex(double(x(1))+double(x(2))*1i);
            P=arrayfun(@(x){cell(x{1})},C{1});
            P=cellfun(@(x){cellfun(pxy,x)},P);
            E=double(C{2});
            V=cellfun(pxy,cell(C{3})).';
            PE=arrayfun(@(e){P(e==E)},unique(E).');
            PE=cellfun(@connectPolylines,PE,UniformOutput=false);
            obj.CrossingTable.Position=V;
            obj.calcEdgeTable;
            obj.calcEdgeDirection;
            obj.calcStrandTable;

            TS=obj.StrandTable;
            TE=obj.EdgeTable;
            E=unique(E);
            E2idx=dictionary(E,1:length(E));
            for si=1:height(TS)
                edges=E2idx(TS.Edges{si});
                % reverse edge's segments coordinates array if the ends are not matching
                % or the tangent directions of both ends are not matching
                % vec=diff(TE.Position{edges(end-[0,1])});
                if isscalar(edges)
                    continue
                end
                edges(end+(1:2))=edges(1:2); % close the loop

                for ei=2:length(edges)
                    flag=0;
                    segs=PE{edges(ei)};
                    vend=abs(TE.Crossing(edges(ei),:));
                    if segs(end)~=V(vend(2))
                        % if the end of the segment does not match the coordinate of the crossing
                        flag=1;
                    elseif vend(1)==vend(2)
                        % if the edge is a loop, check the tangent direction
                        if ei<length(edges)
                            % need to shrink loop segment to avoid collision
                            closedLoopSize=2/3; % size of the closed loop segment
                            segs=segs(1)+closedLoopSize*(segs-segs(1));
                        end
                        segs2=PE{edges(ei-1)};
                        z=(segs2(end)-segs2(end-1))/(segs(2)-segs(1));
                        assert(z~=0)
                        if abs(angle(z))>0.1
                            flag=1;
                        end
                        if length(segs)==4
                            mp=segs(2)+segs(3)-segs(1);
                            segs=[segs(1:2) mp segs(3:4)];
                        end
                    end
                    if flag
                        % reverse the segment if flag is set and cast double to complex
                        PE{edges(ei)}=complex(flip(segs));
                    else
                        PE{edges(ei)}=complex(segs);
                    end
                end
            end
            TE.Position=PE;
            TE.ID=E';
            obj.EdgeTable=TE;
            TRE=obj.calcREdgeTable;
            for ei=1:height(TRE)
                arc=TRE.Arc{ei};
                [~,idx]=ismember(arc([1:end-1;2:end].'),TE.Crossing,'rows');
                [~,idx2]=ismember(arc([2:end;1:end-1].'),TE.Crossing,'rows');
                % issue: calc Positionにおいてedgeの順序が反転している
                idx(idx==0)=idx2(idx==0);
                TRE.Position{ei}=connectPolylines(TE.Position(idx));
            end
            obj.REdgeTable=TRE;
            % [~,ori]=obj.getPDCode;
            % isVertex=logical(ori);
            % TRE=
            % TRE.Position=cellfun(@(e){connectPolylines(TE.Position(E2idx(e)))},TRE.arc);

            function C = connectPolylines(PE)
                %CONNECTPOLYLINES Reorder and concatenate polylines in cell array PE
                %   C is a single row vector of complex points, with duplicate
                %   connection points removed.
                % issue: fail if coordinates of a segment are flipped
                pidx=cellfun(@(e)length(e)==2&&e(1)==e(2),PE);
                pcoord=[PE{pidx}];
                PE(pidx)=[]; % delete a point
                if length(PE)==1&&PE{1}(1)==PE{1}(end)
                    C=PE{1};
                    return
                end
                starts=cellfun(@(s)s(1),PE);
                ends=cellfun(@(s)s(end),PE);

                p_is_included=ismember(pcoord,[starts,ends]);
                [~,next]=ismember(ends,starts);
                e_cur=1;
                eidx=[];
                if all(next~=0)
                    eidx=1:length(PE);
                    %issue: loopのとき，vertexの位置がわからない
                else
                    while e_cur~=0
                        eidx=[e_cur,eidx];
                        [~,e_cur]=ismember(e_cur,next);
                    end
                    e_cur=next(1);
                    while e_cur~=0
                        eidx(end+1)=e_cur;
                        e_cur=next(e_cur);
                    end
                end
                assert(length(eidx)==length(PE)&&all(p_is_included), ...
                    'connectPolylines:NoMatch', 'No matching single polyline found')
                PEsorted=PE(eidx);
                C=cellfun(@(e)e(1:end-1),PEsorted,UniformOutput=false);
                C=[C{:},PEsorted{end}(end)];


                % % https://chatgpt.com/share/686d8571-c200-8003-8380-463bf957ed03
                % % issue: if arcs are not simple curve (self cross), it fail
                % % Filter out cells with fewer than 2 points
                % mask = cellfun(@(v) numel(v) > 1, PE);
                % polys = PE(mask);
                % n = numel(polys);
                %
                % % Precompute start and end of each polyline
                % starts = cellfun(@(v) v(1), polys);
                % ends   = cellfun(@(v) v(end), polys);
                % used = false(1, n);
                %
                % % Initialize with the first polyline
                % idx0 = 1;
                % used(idx0) = true;
                % C = polys{idx0};
                % current_end = C(end);
                % current_start = C(1);
                %
                % % Iteratively find and append matching polylines
                % for k = 2:n
                %     found = false;
                %     for j = find(~used)
                %         if starts(j) == current_end          % matches start
                %             next_poly = polys{j};
                %         elseif ends(j)==current_start
                %
                %         elseif ends(j) == current_end      % matches end, reverse needed
                %             next_poly = flip(polys{j});
                %             disp("flipped")
                %         else
                %             continue;
                %         end
                %         used(j) = true;
                %         % Remove duplicate connection point
                %         next_poly = next_poly(2:end);
                %         % Append
                %         C = [C, next_poly];
                %         current_end = C(end);
                %         found = true;
                %         break;
                %     end
                %     if ~found
                %         error('connectPolylines:NoMatch', 'No matching polyline found at step %d.', k);
                %     end
                % end
                % if C(end)==C(end-1)
                %     C=C(1:end-1);
                % end
            end

        end
        function [C3,P,Cr,C4]=assignPositions2(obj,gap)
            % legacy code
            if nargin < 2
                gap = 0.1;
            end
            Pname=obj.sageName;
            cmd=readlines(fullfile(TopologyConfig.ProjectPath,'Manifold','VirtualLink','plotfun.py'));
            cmd=strrep(cmd,"Link000",Pname);
            cmd=strrep(cmd,"GAP",string(gap));
            SageWrapper.H.exec(cmd(1:end-1))
            C=cell(SageWrapper.H.exec(cmd(end)));
            P=C{1};
            Cr=SageWrapper.toDouble(C{3});
            C4=cell(C{4});
            C=cell(C{2});
            C2=cellfun(@cell,C,UniformOutput=false);
            pxy=@(x)complex(double(x(1))+double(x(2))*1i);
            C3=cellfun(@(x){cellfun(pxy,x)},C2);

        end
    end
    methods
        function [PDCode,orientation,Ncircle]=getPDCode(obj)
            if obj.virtualFlag
                obj.calcVG2G;
            end
            if obj.formatFlag("PD")
                PDCode=obj.PDCode;
                orientation=obj.orientation;
            else
                SW=SageWrapper.H;
                obj.setSageLink;
                PDCode=SW.exec(sprintf("%s.pd_code()",obj.sageName));
                PDCode=cellfun(@(x){double(x)},cell(PDCode));
                PDCode=vertcat(PDCode{:});
                orientation=double(SW.exec(sprintf("%s.orientation()",obj.sageName)));
                obj.PDCode=PDCode;
                % obj.orientation=orientation; %isVirtualCrossingを経由する必要
            end
            Ncircle=obj.Ncircle;
        end
        function [GaussCode,orientation]=getGaussCode(obj,cutsufix)
            if nargin==1, cutsufix=true; end
            if obj.virtualFlag
                obj.calcVG2G;
            end
            if obj.formatFlag("G")
                GaussCode=obj.GaussCode.';
                orientation=obj.orientation;
            else
                SW=SageWrapper.H;
                obj.setSageLink;
                GaussCode=SW.exec(sprintf("%s.gauss_code()",obj.sageName));
                GaussCode=cellfun(@(x){double(x)},cell(GaussCode));
                orientation=double(SW.exec(sprintf("%s.orientation()",obj.sageName)));
                if ~obj.formatFlag("TS")
                    obj.StrandTable=repmat(obj.StrandTableInit,length(GaussCode),1);
                end
                obj.GaussCode=GaussCode;
                obj.orientation=orientation;
                obj.formatFlag("G")=true;
            end
            if ~isempty(obj.isCut)&&cutsufix
                % If there are cut edges, append 0 to the end of each Gauss code
                GaussCode(obj.isCut)=cellfun(@(x){[x,0]},GaussCode(obj.isCut));
            end
        end
        function [GaussCode,orientation]=getRGaussCode(obj,skip)
            if obj.formatFlag("VG")
                GaussCode=obj.RGaussCode;
                orientation=obj.ROrientation;
            else
                if obj.formatFlag("G")
                    GaussCode=obj.GaussCode;
                    orientation=obj.orientation;
                else
                    [GaussCode,orientation]=obj.getGaussCode();
                end
                rc=find(orientation);
                obj.RCrossingMap=rc; %?
                dic=dictionary(1:length(orientation),0);
                dic(rc)=1:length(rc);
                dic(-dic.keys)=-dic.values;
                for i = 1:length(GaussCode)
                    gc=dic(GaussCode{i});
                    GaussCode{i} = gc(gc~=0); % Remove R crossings
                end
                orientation=orientation(rc);
                obj.RGaussCode=GaussCode;
                if nargin<2||~skip
                    obj.calcRCrossingTable; % duplicate?
                    obj.ROrientation=orientation;
                end
                obj.formatFlag("VG")=true;
            end
            if ~isempty(obj.isCut)
                % If there are cut edges, append 0 to the end of each Gauss code
                GaussCode(obj.isCut)=cellfun(@(x){[x,0]},GaussCode(obj.isCut));
            end
        end
        function [RPD,vori,Nc]=getRPDCode(obj)
            % Get the virtual PD code and orientation
            if obj.formatFlag("RPD")
                RPD=obj.RPDCode;
                vori=obj.ROrientation;
                Nc=obj.Ncircle;
            else
                obj.calcRCrossingTable;
                [vgc,vori]=obj.getRGaussCode();
                % Convert virtual Gauss code to virtual PD code
                if isempty(vori)
                    % Handle empty case (unknot components only)
                    RPD = [];
                    Nc = length(vgc);
                else
                    % Get all crossing vertices
                    NV=length(vori);

                    % Count unknot components (empty Gauss codes)
                    Nc = sum(cellfun(@isempty, vgc));

                    % Initialize PD code matrix
                    pd2 = nan(NV, 4);
                    edgeCounter = 0;
                    for strandIdx = 1:length(vgc)
                        if isempty(vgc{strandIdx})
                            continue; % Skip unknot components
                        end

                        strand = vgc{strandIdx};
                        strandEdges = edgeCounter + [1:length(strand),1];
                        edgeCounter = edgeCounter + length(strand);
                        strand(end+1)= strand(1); % Close the loop for the strand
                        % Process each crossing in this strand
                        for i = 1:length(strand)-1
                            v = strand(i);
                            pd2(abs(v),[1,3]+(v>0)) = strandEdges(i+[0,1]);
                        end
                    end
                    if any(isnan(pd2(:)))
                        % Handle case where some crossings might not be properly connected
                        warning('possible issue with virtual Gauss code');
                    end
                    % Convert to PD code format
                    j=find(vori>0);
                    RPD=pd2;
                    RPD(j,[4,2])=RPD(j,[2,4]);

                end
                obj.RPDCode=RPD;
                obj.ROrientation=vori;
                obj.formatFlag("RPD")=true;
            end
        end
        function G=getGraphA(obj)
            G=graph;
            [PD,ori,Ncircle]=obj.getRPDCode();
            edges=unique(PD(:)).';
            NV=length(ori);
            arr=[4,1,2,3,4,1];
            for e=edges
                idx=find(PD==e,2);
                v=mod(idx-1,NV)+1;
                d=(idx-v)/NV;
                G=G.addedge(3*(v(1)-1)+arr(d(1)+(1:3)), ...
                    3*(v(2)-1)+arr(d(2)+(1:3)),e);
            end
        end
        function Hmat=getMatrixHMove(obj)
            %  H-move is a move that changes weights around a crossing
            %  Hmat is (Nedge,Nvertex) integer matrix
            obj.calcRCrossingTable;
            obj.calcREdgeTable;
            rhm=obj.RCrossingTable.HeadMap;
            NV=size(rhm,1);
            v=(1:NV).';
            rhm2=pred(rhm,[-v,v]);
            vend=obj.REdgeTable.Crossing;
            NE=length(vend);
            idx=zeros(NV,4);
            [~,idx(:)]=ismember([rhm2(:),[-v;v];[-v;v],rhm(:)],vend,"rows");
            idx=idx+NE*(v-1);
            Hmat=zeros(NE,NV);
            Hmat(idx(:,1:2))=-1;
            Hmat(idx(:,3:4))=Hmat(idx(:,3:4))+1;
        end
        function G=getDSDiagram(obj)
            % 3-valent directed graph PDにしたほうがよい
            G=digraph;
            [PD,ori,Ncircle]=obj.getRPDCode();
            edges=unique(PD(:)).';
            NV=length(ori);
            arr=[4,1,2,3,4,1];
            for e=edges
                idx=find(PD==e,2);
                v=mod(idx-1,NV)+1;
                d=(idx-v)/NV;
                if d(1)==1||d(1)==2&&ori(v(1))>0||d(1)==4&&ori(v(1))<0
                    % if the edge is a virtual crossing, add it as a directed edge
                    G=G.addedge(3*(v(1)-1)+arr(d(1)+(1:3)), ...
                        3*(v(2)-1)+arr(d(2)+(1:3)),e);
                else

                end
                % if the edge is a real crossing, add it as an undirected edge
                G=G.addedge(3*(v(1)-1)+arr(d(1)+(1:3)), ...
                    3*(v(2)-1)+arr(d(2)+(1:3)),e);
            end
            % [~,vjp]=ismember(PD2(:,3),PD2(:,2));
            % vertexの接続関係を示すPDのデータ構造があれば便利かもしれない
            % [~,vjm]=ismember(PD2(:,3),PD2(:,1));
            % PD3(:,3)=vjp-vjm;
            % [~,vjp]=ismember(PD2(:,4),PD2(:,2));
            % [~,vjm]=ismember(PD2(:,4),PD2(:,1));
        end
        function [tbl,G]=getDiskTable(obj)
            % getDiskTable returns a table of disks

            % dic: the orientation of a vertex -> the edge correspondence of the dual tetrahedron
            persistent dic
            if isempty(dic) %setup dictionary for connection of disk boundaries
                dp([1,2,4,5,9,12,7,6,8,11,10,3])=[1:6,1:6];
                dm([1,2,4,5,9,12,7,10,3,11,6,8])=[1:6,1:6];

                dic=dictionary([1,-1],{dp,dm});
            end
            T=obj.calcREdgeTable;
            G=digraph();
            [PD,ori,Ncircle]=obj.getRPDCode();
            edges=unique(PD(:));
            NE=length(edges); % ignoring Ncircle
            vec=dictionary([edges;-edges], ...
                mat2cell([eye(NE);-eye(NE)],ones(2*NE,1),NE));
            tmp=[0;1;0;0;1;0]*(ori>0)+[0;0;1;0;0;1]*(ori<0);
            % G: Node: an edge of a tetrahedron, Edge: a face of a tetrahedron
            G=G.addnode(table(tmp(:),Var="Ndot"));
            for ei=1:height(T) % connect edges of G
                v=T.Crossing(ei,:);
                tmp=6*(abs(v)-1)+ ...
                    [dic{ori(abs(v(1)))}((1:3)+6+3*(v(1)>0));
                    dic{ori(abs(v(2)))}((1:3)+3*(v(2)>0))].';
                G=G.addedge(tmp([1,2,6]),tmp([4,5,3]),ei*[1,1,-1]);
            end
            [cycles,edges_u]=G.allcycles; % find all cycles (disk) in G
            edges_s=cellfun(@(c){G.Edges.Weight(c)'},edges_u);
            % nodes=cellfun(@(c){G.Edges.EndNodes(c,1)'},cycles);
            Ndots=cellfun(@(c)sum(G.Nodes.Ndot(c)),cycles);
            tbl=table(cycles,edges_s,Ndots);
            tmp=cellfun(@(c){sum(vertcat(vec{c}),1)},edges_s);
            % dx: coboundary of disk?  cp: Euler cochain
            tbl.delta=vertcat(tmp{:});
            tbl.cp=1-tbl.Ndots/2;
        end
        function calcVG2G(obj)
            %  calcVG2G converts the virtual Gauss code to the real Gauss code
            % add virtual crossings in order to make the objects valid in visualization
            % issue: just after VG2G, must set sage link (for consistency)
            vgc=obj.RGaussCode;
            vori=obj.ROrientation;
            [gc,ori]=vg2g(vgc,vori);
            obj.setData(reset=false,GaussCode=gc,orientation=ori)
            obj.setSageLink;
        end
        function calcG2VG(obj)
            % calcG2VG update information of vertices
            obj.RCrossingTable
        end
        function calcStrandTable(obj)
            % Calculate the strand table based on edge connections

            [gc,ori_]=getGaussCode(obj,false);
            ori_(ori_==0)=1;
            PD=obj.getPDCode();
            strands=cell(length(gc),1); % Cell array to hold strands
            for i=1:length(gc)
                %  {-nl: 3, -nr:3, +nr:2, +nl:4}
                % strand = arrayfun(@(v) PD(abs(v), 3-(v > 0)*(ori(abs(v)) >= 0)), gc{i});
                strand = arrayfun(@(v) PD(abs(v), 3-(v > 0)*ori_(abs(v))), gc{i});
                strands{i} = strand; % Add the strand to the list
            end
            if height(obj.StrandTable) ~= length(strands)
                % If the number of strands has changed, update the table
                obj.StrandTable = cell2table(num2cell(nan(length(strands), length(obj.StrandTableVariableNames))), ...
                    VariableNames=obj.StrandTableVariableNames);
            end
            obj.StrandTable.Edges = strands; % Assign the strands to the table
            obj.StrandTable.Crossings =gc.'; % Assign the vertices to the table
        end
        function [weight0,weight1,Hperp]=calcWeight(obj,set)
            % calcWeight calculates the integral weight of o-graph
            % Solve equation cp=-dx, where cp is the Euler 2-cochain
            % and dx is the coboundary of 1-cochain of the polyhedron
            % Axis 2 of the arrays correspond to edge cochains.
            % wei0: one particular solution of weight
            % wei1: kernel of coboundary (homogeneous solution of weight)
            % wei0+wei1: affine space of solutions
            % hperp: generators in wei1 perpendicular to H-moves
            % hmoves: generators of H-move
            arguments
                obj
                set (1,1) logical=false % if true, set the special solution as weight
            end
            tbl=obj.getDiskTable;
            SW=SageWrapper.H;
            dx=SW.toStr(tbl.delta,2);
            cp=SW.toStr(-tbl.cp,1);
            hmoves=obj.getMatrixHMove';
            HM=SW.toStr(hmoves,2); % H-move matrix
            cmd=[
                sprintf("bndM = Matrix(ZZ, %s)", dx)
                sprintf("eulerChain = vector(ZZ, %s)", cp)
                sprintf("HMoves = Matrix(ZZ, %s)", HM)
                "weightSpace = ZZ**bndM.ncols()"
                "HMoves = weightSpace.submodule(HMoves)"
                "weight0 = bndM.solve_right(eulerChain) # solve cp = -dx"
                "weight1_ = bndM.right_kernel_matrix()"
                "kerM = weightSpace.submodule(weight1_)"
                "quotM = kerM.quotient(HMoves)"
                "invariants = str(list(quotM.invariants()))"
                "# print(""invariant"",invariants)"
                "Hperp = [list(gen.lift()) for gen in quotM.gens()]"
                "Hperp=Matrix(ZZ,Hperp)"
                "[str(list(weight0)), str(weight1_), invariants, str(Hperp)]"];
            try
                C=SW.exec(cmd);
            catch ME
                if contains(ME.message,'no solution')
                    warning("No solution found for the weight calculation. Returning empty arrays.");
                    weight0=[];weight1=[];Hperp=[];
                    return
                else
                    rethrow(ME);
                end
            end
            weight0=str2num(string(C{1}));
            weight1=str2num(string(C{2}));
            D=str2num(string(C{3}));
            Hperp=str2num(string(C{4}));

            if ~isempty(D)&&D(1)~=0
                warning("weight solution space has torsion: %s",mat2str(D));
            end
            if isempty(weight1)
                weight1=zeros(0,size(tbl.delta,2));
            end
            if set
                obj.setWeight(weight0,"REdge")
            end
            % kerMc=null(kerM,"rational")'; % kerMc: complement space of kerM
            % V_=null([HMoveM;kerMc],"rational"); % V_:complement space of H in kerM
            %     [~,den]=numden(sym(V_(:,i)));
            %     den=lcm(den);
        end
        function ret=calcSpec(obj,arg,param)
            arguments
                obj
                arg string {mustBeMember(arg,["b1","pi1","h1"])}
                param =[]
            end
            switch arg
                case "b1"
                    assert(obj.isClosed,"not impl for non-closed manifold")
                    [~,~,hp]=obj.calcWeight;
                    ret=size(hp,1);
                case "pi1"
                    % MST2021 page19 fundamental group calculus by Koda
                    % generator: e1,...,e_NE
                    % relation: g=l, hg=k ,at a vertex: edges=[g,k,l,h]
                    % The output is equal to the order.
                    % parameter=variable name in sage
                    if isempty(param), param="group"; end
                    vname=param;
                    [RPD,ori,NC]=obj.getRPDCode;
                    NE=max(RPD,[],"all")+NC;
                    gen=join("e"+(1:NE),",");
                    rel1=sprintf("e%d/e%d,",RPD(:,[1,3]).');
                    rel2=sprintf("e%d*e%d/e%d,",RPD(:,[4,1,2])');
                    SW=SageWrapper.H;
                    cmd(1)=sprintf("%s = FreeGroup(%d)",vname,NE);
                    cmd(2)=sprintf("%s=%s.gens()",gen,vname);
                    cmd(3)=sprintf("%s=%s/[%s%s]",vname,vname,rel1,rel2);
                    SW.exec(cmd);
                    ret=string(SW.exec(sprintf("str(%s.order())",vname)));
                    if ret=="+Infinity", ret=inf; end
                    ret=double(ret);
                case "h1"
                    if isempty(param), param="group"; end
                    calcSpec(obj,"pi1",param);
                    SageWrapper.H.exec(sprintf( ...
                        "%s=%s.abelianization_map().codomain()",param,param));
            end
        end
        function ret=calcInvariant(obj,hopfalg)
            % must execute in Execution repository

            if isa(hopfalg,'VectAlg')
                assert(obj.isWeighted)
                % assume that Tp,Tm,S are already calc right after declaration
                Tp=hopfalg.getSC('Tp');
                Tm=hopfalg.getSC('Tm');
                S=hopfalg.getSC('antipode');
                obj.calcREdgeTable;
                [~,ori]=obj.getRGaussCode;
                TRE=obj.REdgeTable;
                expr1=repmat("Tp",1,length(ori));
                expr1(ori<0)="Tm";
                for i=1:length(ori)
                    expr1(i)=expr1(i)+sprintf("{%s}",join(string(4*(i-1)+(1:4)),","));
                end
                expr2=repmat("(S^%d){%d,%d}",1,height(TRE));
                M=abs(TRE.Crossing)*4+(TRE.Crossing>0)-[1,3];
                W=2*TRE.Weight;
                for i=1:height(TRE)
                    expr2(i)=sprintf(expr2(i),W(i),M(i,1),M(i,2));
                end
                expr=join([expr1,expr2]);
                ret=calcTensorExpression(expr,[]);
            elseif strcmp(hopfalg,"UR")
                D=VL2URD(obj);
                b1=obj.calcSpec("b1");
                % D.simplify2;
                % vars=symvar([D.W,D.C]);
                % expr=simplify(subs([D.C,D.W],vars,eps*ones(size(vars))));
                ret=D.trace;
            end
        end
        function [hi,ti]=calcEdgeDirection(obj)
            % Calculate edge directions based on PD code on Sage
            pname=obj.sageName;
            cmd=[...
                "t,h="+pname+"._directions_of_edges()";
                "pd="+pname+".pd_code()";
                "ti = [(1 + pd.index(t[e]))*(1-2*(t[e][2]==e))" + ...
                "                  for e in sorted(t.keys())]";
                "hi = [(1 + pd.index(h[e]))*(1-2*(h[e][0]==e))" + ...
                "                  for e in sorted(h.keys())]";
                "ti,hi"];
            ret=SageWrapper.H.exec(cmd);
            hi=double(ret{1});
            ti=double(ret{2});
            obj.EdgeTable.Crossing=[hi;ti].';
        end
        function calcEdgeTable(obj)
            hm=obj.headMap;
            NC=obj.Ncircle;
            NV=size(hm,1);
            if obj.formatFlag("TE")
                TE=obj.EdgeTable;
            else
                TE=repmat(obj.EdgeTableInit,NC+2*NV,1);
            end
            TE.Crossing=[[[-(1:NV),1:NV].',hm(:)];nan(NC,2)];
            obj.EdgeTable=TE;
        end
        function calcRCrossingTable(obj)
            % Calculate the crossing table for real crossings
            if obj.formatFlag("TRC")
                return;
            end
            [gc,ori]=obj.getRGaussCode(true);
            TRC=repmat(obj.CrossingTableInit,length(ori),1);
            rhm=nan(length(ori),2);
            for si=1:length(gc)
                cgc=gc{si};
                if isempty(cgc)
                    continue; % skip empty strands
                end
                cgc=[cgc,cgc(1)]; % close the loop
                for vi=1:length(cgc)-1
                    rhm(abs(cgc(vi)),1+(cgc(vi) > 0))=cgc(vi+1);
                end
            end

            TRC.Orientation=ori(:);
            TRC.HeadMap=rhm;
            obj.RCrossingTable=TRC;
            obj.formatFlag("TRC")=1;
        end
        function TRE=calcREdgeTable(obj)
            % real/virtual
            % if decide arc= if obj is real
            % issue: getGaussCodeだと計算量が無駄な場合がある．isVirtualか確認

            [gc,ori]=obj.getGaussCode;
            emptyIdx=cellfun(@isempty,gc);
            gc=gc(~emptyIdx);
            aj=0;
            ej=0;
            if obj.formatFlag("TRE")
                TRE=obj.REdgeTable;
            else
                NE=sum(cellfun(@length,gc))+length(emptyIdx);
                TRE=repmat(obj.EdgeTableInit,NE,1);
            end
            TRE.Arc(:)={[]};
            iscut=obj.isCut;
            for si=1:length(gc)
                % issue: if the strand is cut, need to add edge
                cgc=gc{si}(1:end-iscut(si));
                firstVertex=find(ori(abs(cgc)),1);
                cgc=cgc([firstVertex:end,1:firstVertex]);
                cori=ori(abs(cgc));
                for vi=1:length(cgc)-1
                    ej=ej+1;
                    if cori(vi)~=0
                        % if the vertex is real
                        aj=aj+1;
                        TRE.Arc{aj}=cgc(vi);
                        Crossing(aj,1)=cgc(vi);
                    end
                    TRE.Arc{aj}(end+1)=cgc(vi+1);
                    Crossing(aj,2)=cgc(vi+1);
                end
            end
            vj=find(ori);
            TRE=TRE(1:aj,:);
            dic=dictionary(vj,1:length(vj));
            dic(-dic.keys)=-dic.values;

            [~,widx]=ismember(Crossing,TRE.Crossing,'rows');
            weight=[TRE.Weight;0];
            widx(widx==0)=height(TRE)+1;
            TRE.Crossing=dic(Crossing);
            TRE.Weight=weight(widx); % sort weight array
            obj.REdgeTable=TRE;
            obj.formatFlag("TRE")=1;
            obj.RCrossingMap=vj;
            % if ~obj.formatFlag("TRC")
            %     obj.CrossingTable=repmat(obj.CrossingTableInit,length(vj),1);

            % end
        end
    end
    %% plot, display, and disp methods
    methods
        function p=plot(obj,arg)
            % Plot the VirtualLink diagram
            % overview: first, calculate position and .
            arguments
                obj
                arg.virtual=true
                arg.cut (1,:) = [] % cut specified edge indices of strands
                arg.smooth (1,1) logical = true % smooth the curve
                arg.color  = 'b' % color of the knot
                arg.lineWidth (1,1) double = 1 % line width for the knot
                arg.markerSize (1,1) double = 20 % size of markers for crossings
                arg.coordinate (1,1) logical =false % show coordinates and axis grid
                arg.edgeLabel {mustBeMember(arg.edgeLabel,["none","ID","weight","all"])} = "weight" % edge label type
                arg.vertexLabel {mustBeMember(arg.vertexLabel,["none","ID","ori","all"])} = "all" % vertex label type
                arg.vertexStyle {mustBeMember(arg.vertexStyle,["none","real","imag","all"])} = "real" % vertex dot style
                arg.gap (1,1) double = 0.2 % gap between edges
                arg.componentGap (1,1) double = 1 % gap between components
                % how many times to bend the edge to the left
                arg.bendingNumbers (1,:) cell = {} %WARNING: this may freeze the program
            end
            if obj.virtualFlag
                obj.calcVG2G;
            end
            obj.calcREdgeTable;

            [PE,V]=obj.calcPositions(arg);


            ori=obj.orientation;
            obj.calcEdgeDirection;
            if arg.virtual
                obj.calcRCrossingTable;
                TE=obj.REdgeTable;
                TC=obj.RCrossingTable;

            else
                TE=obj.EdgeTable;
                TC=obj.CrossingTable;
            end

            obj.calcStrandTable;
            TS=obj.StrandTable;
            gap=0.2; % gap between edges

            edgeSet=unique(horzcat(TS.Edges{:}));
            dic=dictionary(edgeSet,1:length(edgeSet));
            if isempty(arg.cut)
                arg.cut=0+obj.isCut;
                for si=1:height(TS)
                    if arg.cut(si)&&arg.virtual
                        [~,arg.cut(si)]=ismember(TS.Vertices{si}([end,1]),TE.Crossing,"rows");
                    elseif arg.cut(si)~arg.virtual
                        [~,arg.cut(si)]=ismember(TS.Crossings{si}([end,1]),TE.Crossing,"rows");
                    end
                end
            end
            figure
            hold on
            for ei=1:height(TE)

                segs=TE.Position{ei};
                vend=TE.Crossing(ei,:);
                label=TE.Label(ei);
                if ismissing(label), label=""; end
                if arg.edgeLabel=="ID"||arg.edgeLabel=="all"
                    label=sprintf('E%d',ei);
                end
                if obj.isWeighted&&(arg.edgeLabel=="weight"||arg.edgeLabel=="all")
                    weight=TE.Weight(ei);
                    label=join([label,num2str(weight)],", ") ;
                end
                % calculate middle point of the segment to place the label
                if ~ismember(ei,arg.cut)
                    seg=segs(floor(end/2)+(0:1));
                    mp=mean(seg);
                    arrowhead=uvec(seg(1:2))*[-1+.5i,0,-1-.5i]/7+mp;
                    plot(arrowhead,LineWidth=arg.lineWidth,Color=arg.color)
                    % Place the label at the middle point
                    text(real(mp),imag(mp),label,'VerticalAlignment','top');
                else
                    arrowhead=segs(1)+uvec(segs(1:2))/3+uvec(segs(1:2))*[-1+.5i,0,-1-.5i]/7;
                    plot(arrowhead,LineWidth=arg.lineWidth,Color=arg.color)
                    segs=[segs(1)+[0,1/3]*uvec(segs(1:2)) nan segs(end)-[1/3,0]*uvec(segs(end-1:end))];
                end
                % Why doubles may be casted into double here?
                segs(1)=segs(1)+gap*(vend(1)<=0&&TC.Orientation(abs(vend(1))))*uvec(segs(1:2));
                segs(end)=segs(end)+gap*(vend(2)<=0&&TC.Orientation(abs(vend(2))))*uvec(segs(end-(0:1)));
                if arg.smooth
                    segs=smoothCurve(segs,0.3);
                end
                plot(complex(segs),LineWidth=arg.lineWidth,Color=arg.color)

            end


            isVertex=logical(ori);
            Vr=V(isVertex);
            scatter(real(Vr),imag(Vr),arg.markerSize,arg.color,"filled")
            text(real(Vr)+arg.gap/2,imag(Vr)+arg.gap/2,"V"+string(1:length(Vr)),'VerticalAlignment','bottom')
            xlim padded
            ylim padded
            axis equal
            if ~arg.coordinate, axis off; end
            hold off
            set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0))

            function ret=uvec(v)
                ret=(v(2)-v(1))/abs(v(2)-v(1));
            end
            function ret = smoothCurve(pts, radius)
                % Smooth a curve using spline interpolation for smoother turns
                % pts: complex array of points defining the polyline
                % radius: smoothing radius, controls interpolation density
                if nargin < 2
                    radius = 0.3; % Default smoothing radius
                end

                if length(pts) < 3
                    ret = pts;
                    return;
                end

                % Separate real and imaginary parts
                x = real(pts);
                y = imag(pts);

                if length(pts) < 3
                    ret = pts;
                    return;
                end

                % Create smoothed points with circular arcs at corners
                smooth_pts = pts(1); % Start with first point

                for i = 2:length(pts)-1
                    p1 = pts(i-1);
                    p2 = pts(i);
                    p3 = pts(i+1);

                    % Calculate vectors
                    v1 = p2 - p1;
                    v2 = p3 - p2;

                    % Calculate distances
                    d1 = abs(v1);
                    d2 = abs(v2);

                    if d1 > 0 && d2 > 0
                        % Unit vectors
                        u1 = v1 / d1;
                        u2 = v2 / d2;

                        % Points before and after corner
                        start_pt = p2 - u1 * radius;
                        end_pt = p2 + u2 * radius;

                        % Add line segment to start of arc
                        smooth_pts(end+1) = start_pt;

                        % Create circular arc
                        % Calculate angle between vectors
                        angle_dif=angle(u2)-angle(u1)+[-2*pi,0,2*pi];
                        [~,min_idx]=min(abs(angle_dif));
                        angle_dif=angle_dif(min_idx);
                        if abs(angle_dif) > 0.1 % Only add arc if angle is significant
                            arc_angles = linspace(0,angle_dif,30);
                            % Rotate u1 by intermediate angle and scale by radius
                            arc_pt = p2 + u1/1i*radius*((exp(1i * arc_angles)-1)*sign(angle_dif)-1i);
                            smooth_pts=[smooth_pts, arc_pt];
                        end
                        % Add end point of arc
                        smooth_pts(end+1) = end_pt;
                    else
                        smooth_pts(end+1) = p2;
                    end
                end

                % Add final point
                smooth_pts(end+1) = pts(end);

                ret = smooth_pts;
            end
        end

        function p=plotSage(obj)
            SW=SageWrapper.H;
            fname=fullfile(TopologyConfig.ProjectPath, ...
                "Manifold","VirtualLink","tmpknotimage.png");
            SW.exec(sprintf("img=%s.plot()",obj.sageName))
            SW.exec("img.save_image("""+fname+""")")
            imshow(fname)
            delete(fname)
        end

        function disp(obj)
            feval("disp"+ TopologyConfig.H.VL.disp,obj);
        end
        function disp0(obj)
            % Default display
            builtin('disp',obj);
        end
        function disp1(obj,type)
            % Short text display
            arguments
                obj
                type string {mustBeMember(type,["","G","VG","PD"])} = ""
            end
            disp('<a href="matlab:help VirtualLink">Virtual Link object</a>:')
            if type==""&&obj.formatFlag("G")||type=="G"
                fprintf("Gauss Code:")
                dispC(obj.GaussCode);
                fprintf("Orientation:")
                disp("["+join(string(obj.orientation),", ") + "]")
            elseif type==""&&obj.formatFlag("VG")||type=="VG"
                fprintf("Virtual Gauss Code:")
                dispC(obj.RGaussCode);
                fprintf("Orientation:")
                disp("["+join(string(obj.ROrientation),", ") + "]")
            elseif type==""&&obj.formatFlag("PD")||type=="PD"
                fprintf("PD Code:")
                disp(mat2str(obj.PDCode))
            end
        end
        function disp2(obj)
            % Long text display
            if isnan(obj.C)
                strC="nan";
            else
                strC=obj.C;
                if isempty(symvar(strC))
                    [n,d]=rat(strC,1e-4);
                    strC=sym(n)./sym(d);
                end
                strC=string(strC);
            end
            fprintf("C=%s, V=%s\n",strC,mat2str(obj.V))
            disp(obj.tbl)
        end
        function disp3(obj)
            disp(obj.unicodeArt)
        end
        function disp4(obj)
            T=obj.unicodeArt(true);
            T.weight=double(T.weight);
            disp(T)
        end
        function disp5(obj)
            % Short text display
            if obj.formatFlag("G")
                fprintf("Gauss=")
                fprintf(stringC(obj.GaussCode));
                fprintf(", Orientation=")
                disp("["+join(string(obj.orientation),", ") + "]")
            elseif obj.formatFlag("VG")
                fprintf("RGauss=")
                fprintf(stringC(obj.RGaussCode));
                fprintf(", Orientation=")
                disp("["+join(string(obj.ROrientation),", ") + "]")
            elseif obj.formatFlag("PD")
                fprintf("PD Code:")
                disp(mat2str(obj.PDCode))
            end
        end
        function ret=toString(obj)
            % issue: format flagごとに処理を分岐する
            gc=obj.RGaussCode;
            gc=join(cellfun(@(arr)"["+join(string(arr),", ") + "]",gc),",");
            ori="["+join(string(obj.ROrientation),", ") + "]";
            ret=sprintf("[%s;%s]",gc,ori);
        end
    end
    %% setData
    methods
        function setData(obj,arg)
            % SetData Set data for the VirtualLink object.
            % Input can be (real)GaussCode, (real)PDCode, or ostring.
            % real data means that extra virtual crossings have to be added to be a (virtual) link.
            %
            arguments
                obj
                arg.headMap
                arg.RHeadMap
                arg.orientation
                arg.GaussCode
                arg.RGaussCode
                arg.cut
                arg.PDCode
                arg.Ncircle
                arg.isVertex
                arg.ostring
                arg.DTCode
                arg.table
                arg.odata
                % arg.reset: true if the mathematical structure of the diagram
                % (e.g. crossings, connections, topology) has changed.
                % This forces to initialize tables.
                arg.reset = true
            end
            
            if arg.reset
                obj.formatFlag=obj.formatFlag0;
                obj.sageLinked=false;
            end
            mode=0;
            if isfield(arg,"headMap")
                mode=1;
                assert(length(arg.headMap)==2)
                obj.headMap=arg.headMap;
            elseif isfield(arg,"RHeadMap")
                mode=2;
                assert(length(arg.RHeadMap)==2)
                obj.RHeadMap=arg.RHeadMap;
            elseif isfield(arg,"GaussCode")
                mode=3;
                assert(all(isfield(arg,["orientation"])))
                assert(isa(arg.GaussCode,'cell'),"Gauss code must be cell array of double array")
                NC=sum(cellfun(@length,arg.GaussCode)==0);
                NV=length(arg.orientation);
                if arg.reset
                    obj.StrandTable=repmat(obj.StrandTableInit,length(arg.GaussCode),1);
                end
                obj.CrossingTable=repmat(obj.CrossingTableInit,NV,1);
                obj.GaussCode=arg.GaussCode;
                obj.orientation=arg.orientation;
                obj.formatFlag("G")=true;
                obj.formatFlag("TS")=true;
                obj.formatFlag("TC")=true;
                obj.virtualFlag=false;
                obj.headMap=gc2hm(arg.GaussCode,NV);
                obj.CrossingTable.HeadMap=obj.headMap;
                % issue: cut
            elseif isfield(arg,"RGaussCode")
                mode=4;
                assert(all(isfield(arg,"orientation")))
                if isnumeric(arg.RGaussCode)&&isvector(arg.RGaussCode)
                    arg.RGaussCode={arg.RGaussCode};
                end
                NC=sum(cellfun(@length,arg.RGaussCode)==0);
                NV=length(arg.orientation);
                obj.StrandTable=repmat(obj.StrandTableInit,length(arg.RGaussCode),1);
                obj.RCrossingTable=repmat(obj.CrossingTableInit,NV,1);
                obj.RGaussCode =helperCut( arg.RGaussCode,obj);
                % obj.orientation=arg.orientation;
                obj.ROrientation=arg.orientation;
                obj.formatFlag("VG")=true;
                obj.formatFlag("TS")=true;
                obj.formatFlag("TRC")=true;
                % obj.calcVG2G();
                obj.virtualFlag=true;
                obj.headMap=gc2hm(obj.RGaussCode,NV);
                obj.RCrossingTable.HeadMap=obj.headMap;
            elseif isfield(arg,"DTCode")
                error("DTCode is not supported yet")
            elseif isfield(arg,"PDCode")
                mode=5;
                % assert(all(isfield(arg,["isVertex"])))
                if isfield(arg,["Ncircle"])
                    NC=arg.Ncircle;
                elseif ~isempty(arg.PDCode)
                    NC=0; % if not specified, assume that no unknots are included
                else
                    NC=1; % unknot
                end

                NV=size(arg.PDCode,1);
                obj.CrossingTable=repmat(obj.CrossingTableInit,NV,1);
                obj.PDCode=arg.PDCode;
                % obj.orientation=arg.orientation;
                % obj.orientation=ones(1,NV);
                if isfield(arg,"isVertex")
                    obj.isVertex=arg.isVertex;
                else
                    obj.isVertex=false(1,NV);
                end
                obj.virtualFlag=false;
                obj.calcEdgeTable;
                if ~isfield(arg,"orientation") %isVertex,orientationの統合
                    try
                        obj.formatFlag("PD")=true;
                        obj.setSageLink;
                        arg.orientation=double(SageWrapper.H.exec( ...
                            sprintf("%s.orientation()",obj.sageName)));
                    catch ME
                        obj.formatFlag("PD")=false;
                        rethrow(ME);
                    end
                end
                obj.PDCode=arg.PDCode;
                obj.orientation=arg.orientation;
                obj.headMap=pd2hm(arg.PDCode,arg.orientation);
                obj.CrossingTable.HeadMap=obj.headMap;
                obj.formatFlag("PD")=true;
                obj.formatFlag("TC")=true;
            elseif isfield(arg,"table")
                assert(numel(arg.table)==2)
                data=join(string(arg.table),",");
                data=SageWrapper.H.exec(sprintf( ...
                    "Knots().from_table(%s).pd_code()",data));
                data=cellfun(@(x){double(x).'},cell(data));
                data=[data{:}].';
                obj.setData(PDCode=data);
                % disp(data)
                return
            elseif isfield(arg,"odata")
                try
                    % 入力: O-data形式 {{[vertices]},orientation,{[weight]}}
                    % 出力: obj.GaussCode, obj.orientation を復元

                    obj.setData(RGaussCode=arg.odata{1}, ...
                        orientation=arg.odata{2});
                catch ME
                    assert(isa(arg.odata,"cell")&& ...
                        all(cellfun(@(x)isa(x,"double"),arg.odata{1})) ...
                        && any(length(arg.odata)==[2,3]), ...
                        "odata must be format {{[vertices]},orientation,{[weight]}}")
                    rethrow(ME);
                end
                if length(arg.odata)==3
                    obj.setWeight(arg.odata{3});
                end
                return
            elseif isfield(arg,"ostring")
                % 入力: O-data文字列 'v_{1}^{+r}v_{2}^{-l}⊗v_{3}^{+r}'
                % 出力: obj.V, obj.W を復元
                gc = {}; ori = [];
                pattern = 'v_\{(\d+)\}\^\{([+-][rl0])\}|⊗|,';
                tokens = regexp([char(arg.ostring),'⊗'], pattern, 'tokens');
                cnt=1;
                cgc=[];
                for i = 1:length(tokens)
                    token = tokens{i};
                    if isempty(token)
                        gc{cnt}=cgc;
                        cnt=cnt+1;
                        cgc=[];
                    else
                        idx = str2double(token{1});
                        sign = token{2}(1);
                        dir = token{2}(2);
                        if sign == '+'
                            cgc(end+1) = idx;
                        else
                            cgc(end+1) = -idx;
                        end
                        if dir == 'r'
                            ori(idx) = +1;
                        elseif dir == 'l'
                            ori(idx) = -1;
                        else
                            ori(idx) = 0;
                        end
                    end
                end
                % gc=fliplr(gc); %反転を有効にする設定を用意するか
                obj.setData(RGauss=gc,orientation=ori);
                return
            end
            if arg.reset&&~obj.virtualFlag
                obj.setSageLink()
            end
            obj.Ncircle=NC;

            function gc=helperCut(gc,vlobj)
                % if 0 is inserted in the Gauss code, it is a cut edge
                % and gc code is stored after removing 0
                % cut flag is stored in obj.isCut
                [~,cutIdx]=cellfun(@(x)ismember(0,x),gc);
                if any(cutIdx)
                    % If there are cut edges, append 0 to the end of each Gauss code
                    % issue: where to put 0 (end or start?)
                    isCut=logical(cutIdx);
                    gc(isCut)=arrayfun(@(x,i){x{1}([i+1:end,1:i-1])},gc(isCut),cutIdx(isCut));
                    vlobj.isCut=isCut;
                end
            end
            function hm=gc2hm(gc,NV)
                hm=nan(NV,2);
                for i=1:length(gc)
                    cgc=[gc{i},gc{i}(1)];
                    cgc(isinf(cgc))=[];
                    iscut=cgc(end-1)==0;
                    for j=1:length(cgc)-1-iscut
                        hm(abs(cgc(j)),1+(cgc(j)>0))=cgc(j+1);
                    end
                end
            end
            function hm=pd2hm(pd,ori)
                % Convert PD code to head map
                idx=find(ori>=0);
                pd(idx,[4,2])=pd(idx,[2,4]);
                [~,vm]=ismember(pd(:,3:4),pd(:,1));
                [~,vp]=ismember(pd(:,3:4),pd(:,2));
                hm=[-vm(:,1)+vp(:,1),-vm(:,2)+vp(:,2)];
            end

        end
        function setWeight(obj,weight,type)
            % Set the weight of the strands
            % arguments
            % obj: VirtualLink object
            % weight: weight of the strands, can be a table or a double array
            % type: type of the weight, can be "REdge", "Edge", "Crossing", "RCrossing", or "Strand"
            % As for weight of odata, type= "REdge" and weight is a table
            %  from edge={[source, target]} to weight value, or array of weights of each edge in REdgeTable
            arguments
                obj
                weight % can be a table or a double array
                type {mustBeMember(type,["REdge","Edge","Crossing","RCrossing","Strand"])} = "REdge"
            end
            switch type
                case "REdge"
                    obj.calcREdgeTable;
                    if isa(weight,"table")
                        [~,idx]=ismember(weight.Crossing,obj.REdgeTable.Crossing,'rows');
                        obj.REdgeTable.Weight(idx)=weight.Weight;
                    elseif isa(weight,"double")
                        obj.REdgeTable.Weight=weight(:);
                    end
                otherwise
                    error("Weight for %s is not supported",type)
            end
            obj.isWeighted=true;
        end
    end
end
function ret=succ(hm,v)
    ret=hm(abs(v)+(v>0)*size(hm,1));
end
function ret=pred(hm,v)
    [~,vm]=ismember(v,hm(:,1));
    [~,vp]=ismember(v,hm(:,2));
    ret=vp-vm;
end
function ret=stringC(arg)
    % Convert cell array of strings to a single string
    if isempty(arg)
        ret="{}";
    else
        for i=1:length(arg)
            if isempty(arg{i})
                str(i)="";
            else
                str(i)=join(string(arg{i}),",");
            end
        end
        ret="{["+join(str,"],[")+"]}";
    end
end
function dispC(arg)
    if isempty(arg)
        disp("[[]]")
    else
        str=strings(1,length(arg));
        for i=1:length(arg)
            if isempty(arg{i})
                str(i)="";
            else
                str(i)=join(string(arg{i}),",");
            end
        end
        disp("[["+join(str,"],[")+"]]")
    end
end
function [idxC,idx]=findC(key,arr)
    % find the index of key in arr={[...],...,[...]}]}
    try
        idxC=find(cellfun(@(C)ismember(key,C),arr));
        idx=find(key==arr{idxC});
    catch
        assert(isscalar(key))
        if isempty(idxC)
            error("cannot find the key`%s` in the array`%s` ",string(key), ...
                "[["+join(cellfun(@(x)join(string(x),","),arr),"], [") + "]]")
        end
        assert(isscalar(idxC))

    end
end


function [gc,ori]=vg2g(rgc,ro)


    if isempty(ro)
        assert(~any(cellfun(@numel,rgc)))
        gc=rgc;
        ori=ro;
        return
    end
    % store information about the regions of the virtual graph
    % using vertex/edge/face (regV/regE/regF)
    % regF is a graph object representing the adjacency of regions
    % and its Weights are the edge IDs connecting the regions
    regV={};
    regE={};
    gc=rgc;
    ge=cell(size(gc));
    ori=ro;
    usedV=0;
    cntE=0;
    cntV=max(cellfun(@max,rgc))+1;
    regF=graph;
    for si=1:length(rgc)
        % loop through each strand
        cgc=[rgc{si} 0];
        if isscalar(cgc)
            rgc{si}=[];
            continue
        end
        cntE=cntE+1;
        ge{si}(end+1)=cntE;
        v1=abs(cgc(1));
        if ismember(v1,usedV)
            % start strand sweeping from a used vertex
            [f2,vidx2]=findC(-v1,regV);
            regV{f2}=[cntV+1,regV{f2}([vidx2:end,1:vidx2])];
            regE{f2}=[-cntE,regE{f2}(vidx2:end),regE{f2}(1:vidx2-1),cntE];
            regF=regF.addedge(f2,[f2 f2],[0 cntE]);
            [f1,vidx]=findC(v1,regV);
            regV{f1}=[regV{f1}(1:vidx-1),-v1,0,-v1,regV{f1}(vidx+1:end)];
            regE{f1}=[regE{f1}(1:vidx-1),0,0,regE{f1}(vidx:end)];
            % regF=regF.addedge(f1,f1,0);
            f1=f2;
        else
            % start strand sweeping from a new vertex
            f1=length(regV)+1;
            regV{f1}=[0,-1,0,1]*ro(v1)*cgc(1)+[cntV+1,0,0,0];%?
            regE{f1}=[-cntE,0,0,cntE];
            regF=regF.addedge(f1,[f1 f1],[0 cntE]);
            usedV(end+1)=v1;
        end
        for vi=2:length(cgc)
            % loop through each vertex in the strand
            fd()
            cntE=cntE+1;
            ge{si}(end+1)=cntE;
            v1=abs(cgc(vi));
            if ismember(v1,usedV)
                % pattern 1,2,3,4,6
                [f2,vidx2]=findC(v1,regV);
                fnew=findC(-v1,regV);
                % Find the shortest path from f1 to f2 (and fnew) in the graph
                path = [shortestpath(regF, f1, f2) fnew];
                len=length(path);
                if len==1
                    % pattern 1: connect to another component
                    % f1⋯>f2->fnew
                    [fnew,vidx3]=findC(-v1,regV);
                    regV{f2}=[regV{f2}(1:vidx2-1),-v1,regV{f1}(2:end),-v1,regV{f2}(vidx2+1:end)];
                    regV{fnew}=[cntV+1,-v1,regV{fnew}([vidx3+1:end,1:vidx3-1]),-v1];
                    regV(f1)=[];
                    regE{f2}=[regE{f2}(1:vidx2-1),regE{f1},regE{f2}(vidx2:end)];
                    regE{fnew}=[-cntE,regE{fnew}([vidx3:end,1:vidx3-1]),cntE];
                    regE(f1)=[];
                    % contract the faces f1<->f2
                    ns=regF.neighbors(f1);
                    for n=ns(ns~=f1)'
                        bnd=regF.Edges.Weight(regF.findedge(f1,n));
                        regF=regF.addedge(f2,n*ones(size(bnd)),bnd);
                    end
                    bnd=regF.Edges.Weight(regF.findedge(f1,f1));
                    regF=regF.addedge(f2,f2*ones(size(bnd)),bnd);
                    regF=regF.addedge(fnew,fnew,cntE); % rmnode後に実行していたので，ズレていた
                    regF=regF.rmnode(f1);
                    f1=fnew-(f1<fnew); % issue:faceが消えるとき，face indexが1つずれる対処
                elseif path(1)==path(2)&&vi<length(cgc)
                    % pattern 2,3: connect to vertex
                    vidx=vidx2;
                    [f2,vidx2]=findC(-v1,regV);
                    fnew=length(regV)+1;
                    if vidx2<vidx
                        % pattern 2:
                        regV{fnew}=regV{f1}([vidx2,vidx+1:end]);
                        regV{f1}=regV{f1}([1,vidx2:vidx-1,vidx2,2:vidx2]);
                        regE{fnew}=regE{f1}(vidx:end);
                        regE{f1}=[-cntE regE{f1}([vidx2:vidx-1,1:vidx2-1]) cntE];
                    elseif vidx<vidx2
                        % pattern 3:
                        regV{fnew}=regV{f1}([vidx2,2:vidx-1]);
                        regV{f1}=regV{f1}([1,vidx2:end,vidx2,vidx+1:vidx2]);
                        regE{fnew}=regE{f1}(1:vidx-1);
                        regE{f1}=[-cntE regE{f1}([vidx2:end,vidx:vidx2-1]) cntE];
                    end
                    newbnd=abs(regE{fnew});
                    regF=updateFace(regF,newbnd,f1,f2,fnew,cntE);
                    f1=f2;
                else
                    % pattern 4,6
                    ge{si}(end)=[];
                    cntE=cntE-1;
                    for fj=1:len-1
                        f1=path(fj);
                        f2=path(fj+1);
                        if fj<len-1
                            % add virtual crossing if we have to cut through edges
                            cntE=cntE+1;
                            bnds=regF.Edges.Weight(findedge(regF,f1,f2));
                            bnd=bnds(1); %ここが選択の余地あり
                            vidx=find(bnd==abs(regE{f1}),1)+1; % put new v in the middle of regV{f1}(vidx-[0,1])
                            if TopologyConfig.H.VL_vg2g_disp
                                fprintf("cut through the edge between:[%d,%d]\n",regV{f1}(vidx-[0,1]))
                            end
                            eo=regE{f1}(vidx-1)>0; % orientation of edge
                            es=2*eo-1;
                            regV{f1}=[regV{f1}(1:vidx-1),cntV,regV{f1}(vidx:end)];
                            regE{f1}=[regE{f1}(1:vidx-1-~eo),es*cntE,regE{f1}(vidx-1+eo:end)];
                            vidx2=find(bnd==abs(regE{f2}),1)+1;
                            regV{f2}=[regV{f2}(1:vidx2-1),-cntV,regV{f2}(vidx2:end)];
                            regE{f2}=[regE{f2}(1:vidx2-1-eo),-es*cntE,regE{f2}(vidx2-1+~eo:end)];
                            regF=regF.addedge(f1,f2,cntE);
                            vs=cntV;
                            [~,gcidx]=ismember(cgc(vi),[gc{si} 0]);
                            gc{si}=[gc{si}(1:gcidx-1),-es*vs,gc{si}(gcidx:end)];
                            [sidx,geidx]=findC(bnd,ge);
                            ge{sidx}=[ge{sidx}(1:geidx), cntE, ge{sidx}(geidx+1:end)];
                            gc{sidx}=[gc{sidx}(1:geidx), es*vs, gc{sidx}(geidx+1:end)];
                            ori(cntV)=0;
                            v1=cntV;
                            cntV=cntV+1;
                        elseif vi==length(cgc)
                            % connect ends
                            break
                        else
                            v1=abs(cgc(vi));
                            [~,vidx]=ismember(v1,regV{f1});
                            [f2,vidx2]=findC(-v1,regV);
                        end
                        cntE=cntE+1;
                        fnew=length(regV)+1;
                        % pattern 4: cut through edges and add virtual crossing
                        regV{fnew}=[-v1 regV{f1}(vidx+1:end)];
                        regV{f1}=[-v1 regV{f1}(2:vidx-1)]; %ここの2はなぜ？ 1にしたらどうなる？
                        regV{f2}=[cntV+1 regV{f2}([vidx2:end,1:vidx2])];
                        regE{fnew}=regE{f1}(vidx:end);
                        regE{f1}=regE{f1}(1:vidx-1);
                        regE{f2}=[-cntE regE{f2}([vidx2:end,1:vidx2-1]) cntE];
                        regF=updateFace(regF,abs(regE{fnew}),f1,f2,fnew,cntE);
                        f1=f2;
                        ge{si}(end+1)=cntE;
                    end
                end
            else
                % pattern 5: add new crossing to temporary graph
                usedV(end+1)=v1;
                vs=ro(v1)*cgc(vi);
                regV{f1}=[cntV+1 -vs regV{f1}(2:end) vs];
                regE{f1}=[-cntE regE{f1} cntE];
                regF=regF.addedge(f1,f1,cntE);
            end

        end
        % pattern 6: connect ends of the strand
        vidx=find(regV{f1}==0,1);
        fnew=length(regV)+1;
        regV{fnew}=regV{f1}(vidx+1:end);
        regE{fnew}=regE{f1}(vidx+1:end);
        regV{f1}=regV{f1}(2:vidx-1);
        regE{f1}=regE{f1}([2:vidx-2,1]);
        regF=updateFace(regF,abs(regE{fnew}),f1,f2,fnew,cntE);
        regF=regF.rmedge(f2,f2);
        fd()
        % regF.plot
        % dispC(ge)
        % dispC(gc)
    end
    %
    % function [idxC,idx]=findC(key,arr)
    %     % find the index of key in arr={[...],...,[...]}]}
    %     idxC=find(cellfun(@(C)ismember(key,C),arr));
    %     assert(isscalar(idxC))
    %     idx=find(key==arr{idxC});
    % end
    function regF=updateFace(regF,bnd,f1,f2,fnew,cntE)
        % update the face in regF with the new boundary bnd
        % bnd=unique(bnd);
        [~,notSelfAdj]=ismember(bnd,bnd);
        bnd2=bnd;
        bnd2(notSelfAdj)=[];
        bnd=setdiff(bnd,bnd2);
        [~,eidx_]=ismember(bnd,regF.Edges.Weight);
        assert(all(eidx_~=0))
        BT=regF.Edges.EndNodes(eidx_,:);
        bidx=(BT(:,1)==f1).'+1;
        fbnd=arrayfun(@(x,y)BT(x,y),1:length(eidx_),bidx);
        regF=regF.rmedge(eidx_);
        regF=regF.addedge(fnew,fbnd,bnd);
        regF=regF.addedge(f2,f2,cntE);
        if ~isempty(bnd2)
            [~,eidx_]=ismember(bnd2,regF.Edges.Weight);
            regF=regF.rmedge(eidx_);
            regF=regF.addedge(fnew,fnew*ones(size(bnd2)),bnd2);
        end
        % disp(regF.Edges)
    end
    function fd()
        % for debugging
        if ~TopologyConfig.H.VL_vg2g_disp, return;end
        try
            fprintf("f1=%d, ",f1)
            fprintf("si=%d, vi=%d, v_gc=%d",si,vi,cgc(vi))
            fprintf(", s_gc=%d",ro(abs(cgc(vi))))
        catch
        end
        disp(" ")
        fprintf("edges of strands=    ")
        dispC(ge)
        fprintf("vertices of regions= ")
        dispC(regV)
        fprintf("edges of regions=    ")
        dispC(regE)
        fprintf("gauss code=          ")
        dispC(gc)
        disp("connections of face: [f_s;f_t;edge]")
        T=regF.Edges{:,:}';
        for ii=1:size(T,2)
            if T(3,ii)==0, continue; end
            [si_,vi_]=findC(T(3,ii),regE);
            T(4:5,ii)=regV{si_}(mod(vi_-[1,0],length(regV{si_}))+1);
        end
        disp(T)

    end
end