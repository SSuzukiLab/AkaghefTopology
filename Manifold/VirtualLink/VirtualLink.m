classdef VirtualLink<handle&matlab.mixin.Copyable
    % VirtualLink diagram data tables (structural representation for knots and links)
    % Coordinate system: +x → right, +y → up
    % vertices ID start from 1 to NV, where NV is the number of vertices
    % edges index start from 1 to 2NV+NC, where NC is the number of components
    % But edges ID may not start from 1.

    properties

        CrossingTable
        EdgeTable
        StrandTable
        RegionTable
    end
    properties
        formatFlag  =dictionary(["G" "PD" "VG" "VPD"],false(1,4))
        sageLinked =false
        virtualFlag (1,1)=true
    end
    properties
        name
        description
    end
    properties
        GaussCode (1,:) cell % cell array of Gauss code for each component
        orientation (1,:) double% +1/-1 for right/left type crossings
        virtualOrientation (1,:) double % +1/-1 for right/left type real crossings
        PDCode (:,4) double% NV by 4 matrix of PD code
        Ncircle (1,1) double% the number of circle(unknot) components
        DTCode %not supported
        virtualGaussCode (1,:) cell % Gauss code not including virtual crossings
        virtualPDCode (:,4) double % PD code not including virtual crossings
        cutFlag (1,:) logical  % true if the strand is cut
    end
    properties (Dependent)
        isVirtual (1,:) logical % 1 for virtual crossings, 0 for real crossings
        sageName (1,1) string % variable name in Sage
        sageOrientation (1,:) double % orientation in Sage
    end

    properties(Constant)
        formatFlag0  =dictionary(["G" "PD" "VG" "VPD"],false(1,4))
        CrossingTableVariableNames = {
            'ID',           % Unique crossing identifier
            'Label',        % Optional symbolic label
            'Orientation',  % +1/-1 (r/l type) for crossing; 0 if virtual
            'EdgeID',       % [a,b,c,d] edges in PD code
            'Position',     % x+yi coordinate of crossing center
            'Region',       % Coarse region/zone for layout
            }

        EdgeTableVariableNames = {
            'ID',           % Unique edge identifier
            'Label',        % Optional edg label
            'CrossingID',       % Source crossing ID [source,target](or NaN for free edg)
            'StrandID',     % Associated strand ID
            'Position',     % array of complex coordinates for drawing
            'IsVirtual'     % true if arc involves only virtual crossings
            }

        StrandTableVariableNames = {
            'ID',           % Strand ID (component index)
            'Label',         % Optional strand label (e.g., 'L1')
            'Edges',      % Cell array of edge IDs forming the strand
            'Vertices',  % Cell array of vertex IDs forming the strand
            'Color'         % Optional color
            }

        RegionTableVariableNames = {
            'ID',           % Unique region identifier
            'BoundaryEdges', % List of edge IDs forming the boundary
            'Label',        % Optional region label
            'IsOuter',      % true if it is the unbounded face
            'Zone'          % Coarse position or layout group
            }
        % moveType =["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"]% supported move types
    end

    methods
        %% Validation
        function validateStructure(obj)
        end

        function isClosed(obj)
        end


        %% Moves
        function convert(obj,type,arg)
            arguments
                obj
                type (1,1) string {mustBeMember(type,["mirror","mirrorMfd","reverse"])}  % move type
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
                        [gc,ori]=obj.getVirtualGaussCode;
                    else
                        error("not implemented")
                    end
                    gc=cellfun(@(x){-x},gc);
                    ori=-ori;
                    if flag(1)
                        obj.setData(Gauss=gc,orientation=ori);
                    elseif flag(2)
                        obj.setData(virtualGauss=gc,orientation=ori);
                    end
                case "reverse"
                    error("reverse move is not implemented yet")
            end
        end
        function tbl=movable(obj,type,argstr,arg)
            % MOVABLE Get a table of valid moves for the specified type
            % See `move` for the meaning of arguments and move types.
            % if the table returned is empty, no valid moves are found.
            % tbl: table with columns "detail" and "vertex" for the move
            arguments
                obj
                type (1,1) string {mustBeMember(type,["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"])}  % move type
                argstr struct =struct.empty;
                arg.detail (1,1) string="" % detailed move type
                arg.v (1,:) double =[]% vertex indices for move
                arg.e (1,:) double =[]% edge indices for move
                arg.strict (1,1) logical =false % strict mode, fix the order of vertices
            end
            if ~isempty(argstr)
                arg=argstr;
            end
            tbl=table;
            V=arg.v;
            strict=arg.strict;
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
                    [gc,ori]=obj.getVirtualGaussCode;
                    if N==2
                        T=load("MPmoveData.mat").T_MP_L;
                    else
                        T=load("MPmoveData.mat").T_MP_R;
                    end
                    if strlength(arg.detail)>0
                        T=T(find(T.detail==arg.detail),:);
                    end
                    if length(ori)<max(V)
                        return; % no valid moves if vertex IDs exceed
                    end
                    cur=ori(V);
                    X=searchMoves(T,gc,V,cur);
                    % disp(X)
                    if ~isempty(X)
                        tbl=X(:,["detail","vertex"]);%消す必要無い？
                    end
                    % tbl.Properties.VariableNames=["name","vertex"];
                    % disp(T)
                case "PS"
                    % 2 vertices or 2 edge for Pure Sliding move (inverse)
                    assert(length(V)==2||length(arg.e)==2)
                    [gc,ori]=obj.getVirtualGaussCode;
                    if length(V)==2
                        assert(V(1)~=V(2)&&V(1)*V(2)>0)
                        [sj,vj]=findC(V(1),gc);
                        cgc=gc{sj}([end,1:end,1]);
                        isadj_p=cgc(vj)==V(2)|| cgc(vj+2)==V(2);
                        [sj,vj]=findC(-V(1),gc);
                        cgc=gc{sj}([end,1:end,1]);
                        isadj_n=cgc(vj)==-V(2)|| cgc(vj+2)==-V(2);
                        if isadj_p&&isadj_n
                            tbl=table("PS",V,VariableNames=["detail","vertex"]);
                        end
                    else
                        error("PS move is not implemented yet")
                    end
                    N=length(V);
                    load("PSmoveData.mat","T_PS_R");
                    cur=ori(V);
                    X=searchMoves(T_PS_R,gc,V,cur);
                    if ~isempty(X)
                        tbl=X(:,["detail","vertex"]);
                    end
                case "CP"
                    assert(any(length(V)==[3,5]))
                case "02" % move 0-2
                    assert(length(V)==2||length(arg.e)==1) % 2 vertices or 1 edge for 0-2 move
                case "H" % handle move
                    assert(length(V)==1) % 2 vertices or 1 edge for handle move
                case "BMP" % Bumping move
                    assert(length(V)==2||length(V)==1) % 2 vertices for Bumping move
                    if length(V)==2
                        isPositiveV=V(1)>0;
                        assert((V(1)~=V(2))&&(V(2)>0==isPositiveV))
                        [gc,ori]=obj.getVirtualGaussCode;
                        [sj1,vj1]=findC(V(1),gc);
                        [sj2,vj2]=findC(V(2),gc);
                        sgn=ori(abs(V));
                        if sj1==sj2 && (vj1==vj2+1|| vj1==vj2-1)&& sgn(1)*sgn(2)==-1
                            % check if the vertices are adjacent
                            V=V((3+sgn*(2*isPositiveV-1))/2); % V has orientation [1,-1]
                            tbl=table("",V,VariableNames=["detail","vertex"]);
                        end
                    else
                        error("BMP move is not implemented yet")
                    end
                case "B02" % Bumping 0-2 move
                    assert(length(V)==2)
                    if length(V)==2
                        assert(V(1)~=V(2))
                        [gc,ori]=obj.getVirtualGaussCode;
                        sgn=ori(V);
                        if sgn(1)*sgn(2)~=-1 % opposite orientation
                            return; % no valid moves
                        elseif sgn(1)==-1
                            V=V([2,1]); % make sure V(1) is 1
                        end
                        % to check if a strand made of -V is included
                        % if any(cellfun(@(x)isequal(x,-V),gc))
                        if any(cellfun(@(x)isequal(x,-V)||isequal(x,V),gc))
                            tbl=table("",V,VariableNames=["detail","vertex"]);
                            % elseif any(cellfun(@(x)isequal(x,-V([2,1])),gc))
                        elseif any(cellfun(@(x)isequal(x,-V([2,1]))||isequal(x,V([2,1])),gc))
                            tbl=table("",V([2,1]),VariableNames=["detail","vertex"]);
                        end
                    else
                        error("B02 move is not implemented yet")
                    end
                otherwise
            end
            % issue: 同一データでも適用位置によって複数解が出てくる可能性を考慮する必要
            tbl=unique(tbl,'rows');
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
            % arg.detail: detailed move type 
            % for MP move, specify A1~D4 types of MP moves
            % for BMP move, specify detail=0 or detail=1 for the orientation of the new circle
            arguments
                obj
                type (1,1) string {mustBeMember(type,["R1","R2","R3","MP","PS","CP","02","CP","H","BMP","B02"])}  % move type
                argstr struct =struct.empty;
                arg.detail (1,1) string="" % detailed move type
                arg.v (1,:) double % vertex indices for move
                arg.e (1,:) double % edge indices for move
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
            V=arg.v;
            switch type
                case "R1"
                    error("R1 move is not implemented yet")
                case "R2"
                    error("R2 move is not implemented yet")
                case "R3"
                    error("R3 move is not implemented yet")
                case "MP"
                    N=length(V);
                    [gc,ori]=obj.getVirtualGaussCode;
                    load("MPmoveData.mat","T_MP_L","T_MP_R");
                    if N==2
                        idx=find(T_MP_L.detail==arg.detail,1);
                        T0=T_MP_L(idx,:);
                        T1=T_MP_R(idx,:);
                    else
                        idx=find(T_MP_R.detail==arg.detail,1);
                        T0=T_MP_R(idx,:);
                        T1=T_MP_L(idx,:);
                    end
                    [gc,ori]=convertGaussCode(gc,ori,T0,T1,V);
                    obj.setData(virtualGauss=gc,orientation=ori);

                case "PS"
                    if ~isempty(V)
                        [gc,ori]=obj.getVirtualGaussCode;
                        idx=double(tbl.detail(1));
                        load("PSmoveData.mat","T_PS_R");
                        T0=T_PS_R(idx,:);
                        load("PSmoveData.mat","T_PS_L");
                        [gc,ori]=convertGaussCode(gc,ori,T0,T_PS_L,V);
                        obj.setData(virtualGauss=gc,orientation=ori);
                    else

                        error("PS move is not implemented yet")
                    end


                case "CP"
                    error("CP move is not implemented yet")
                case "02"
                    error("02 move is not implemented yet")
                case "H"
                    error("H move is not implemented yet")
                case "BMP"
                    if length(V)==2
                        if strlength(arg.detail)==0
                            error("specify the parameter detail=0,or detail=1")
                        else
                            orientationOfNewCircle=double(arg.detail);
                        end
                        % issue: B02 move
                        tbl=obj.movable("BMP",v=V);
                        assert(height(tbl)==1&&isequal(tbl.vertex,V),"No valid BMP move found");
                        [gc,ori]=obj.getVirtualGaussCode;
                        s0=double(tbl.detail);
                        [sj1,vj1]=findC(-V(1),gc);
                        gc{sj1}=gc{sj1}([vj1+1:end,1:vj1-1]);
                        [sj2,~]=findC(-V(2),gc);
                        [sj0,vj0]=findC(V(1),gc);
                        vnew=max([gc{:}])+1;
                        gc{sj0}(vj0)=-vnew;
                        vj0=find(gc{sj0}==V(2),1);
                        gc{sj0}(vj0)=[];
                        vj2=find(gc{sj2}==-V(2),1);
                        if orientationOfNewCircle
                            gc_circle=[-V(1:2),vnew];
                            ori([abs(V),vnew])=[-1,1,1];
                        else
                            gc_circle=[-V([2,1]),vnew];
                            ori([abs(V),vnew])=-[-1,1,1];
                        end
                        if sj1==sj2 % if the vertices V are in the same strand, split into two strands
                            gc(length(gc)+(1:2))={[V(1),gc{sj2}(1:vj2-1)],[V(2),gc{sj2}(vj2+1:end)]};
                            gc{sj2}=gc_circle;
                        else % if the vertices are in different strands, form one strand
                            gc{sj2}=[V(2),gc{sj2}([vj2+1:end,1:vj2-1]),V(1),gc{sj1}];
                            gc{sj1}=gc_circle;
                        end
                        obj.setData(virtualGauss=gc,orientation=ori);
                    else
                        error("BMP move inverse is not implemented yet")
                    end
                case "B02"
                    if length(V)~=2
                        error("B02 move is not implemented yet")
                    else
                        assert(all(V>0))
                        tbl=obj.movable("B02",v=V);
                        assert(height(tbl)==1,"No valid BMP move found");
                        [gc,ori]=obj.getVirtualGaussCode;
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
                        obj.setData(virtualGauss=gc,orientation=ori);
                    end
            end
            function [gc,ori]=convertGaussCode(gc,ori,T0,T1,V)
                % Convert Gauss code based on the move
                % T0: original Gauss code data
                % T1: new Gauss code data
                % V: vertices involved in the move
                NV_=max(cellfun(@(x)max(x),gc));
                if T0.NV<T1.NV
                    % if the number of vertices increases after the move
                    V(T0.NV+1:T1.NV)=NV_+(1:T1.NV-T0.NV);
                    dic2=@(x)x;
                    ori(V)=T1.ori;
                elseif T0.NV>T1.NV
                    vrem=setdiff(1:NV_,V(end-(0:T0.NV-T1.NV-1)));
                    dic2=dictionary(vrem,1:NV_-(T0.NV-T1.NV));
                    dic2(-dic2.keys)=-dic2.values;
                    dic2(0)=0;
                    ori(V(1:T1.NV))=T1.ori;
                    ori=ori(vrem);
                else
                    dic2=@(x)x;
                    ori(V)=T0.ori;
                end
                dic1=dictionary([1:length(V),-1:-1:-length(V)],[V,-V]);
                % position of start point of each fragment
                [sj_,vj_]=arrayfun(@(x)findC(dic1(x),gc),T0.gcFirst);

                for s=1:length(gc)
                    % s:current component
                    % update Gauss code for each component
                    k=find(sj_==s);
                    if isempty(k)
                        gc{s}=dic2(gc{s}); % no related crossings in this component
                        continue; % no related crossings in this component
                    end
                    [~,sortidx]=sort(vj_(k));
                    k=k(sortidx);
                    % tmp: temporary Gauss code [gc,T1,gc,T1,...,gc] T1 start at vj_(k)
                    % dic1: dictionary to convert T1.gc to V
                    % dic2: delete and shift indices of crossings, dic: convert to new indices T1->gc
                    tmp=gc{s}(1:vj_(k(1))-1);
                    for i=1:length(k)-1
                        tmp=[tmp,dic1(T1.gc{k(i)}),gc{s}(vj_(k(i))+length(T0.gc{k(i)}):vj_(k(i+1))-1)];
                    end
                    vjend=vj_(k(end))+length(T0.gc{k(end)});
                    tmp=[tmp,dic1(T1.gc{k(end)}),gc{s}(vjend:end)];
                    if vjend>length(gc{s})+1
                        % if at the end of the Gauss code, restart from the beginning
                        tmp(1:(vjend-length(gc{s})-1))=[];
                    end
                    gc{s}=dic2(tmp);
                end
                if length(T1.gc)>length(T0.gc)
                    % if the number of strands increases after the move, add new strands
                    gc=[gc,cellfun(@(x){dic1(x)},T1.gc(length(T0.gc)+1:end))];
                elseif length(T1.gc)<length(T0.gc)
                    % if the number of strands decreases after the move, remove strands
                    sj_remove=sj_(length(T1.gc)+1:end);
                    % strand to remove must be empty
                    assert(all(cellfun(@(x)isempty(x),gc(sj_remove))))
                    gc(sj_remove)=[]; % remove empty strands
                end
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
        function ret=get.isVirtual(obj)
            % Return a logical array indicating virtual crossings
            ret=~logical(obj.orientation);
        end
        function set.isVirtual(obj,arg)
            % isVirtual is read-only property
            obj.orientation(logical(arg))=0; % set orientation to 0 for virtual crossings
        end
        function ret=get.sageOrientation(obj)
            % Return orientation in Sage format
            % convert 0 to 1 for Sage compatibility
            ret=obj.orientation;
            ret(ret==0)=+1;
        end
        function set.sageOrientation(obj,arg)
            % Set orientation from Sage format
            jreal=~obj.isVirtual;
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
        function [PE,V,E]=calcPositions(obj,gap)
            if nargin < 2
                gap = 1; % Default gap between components
            end
            Pname=obj.sageName;
            cmd=readlines(fullfile(TopologyConfig.ProjectPath, ...
                'Manifold','VirtualLink','allocate_pos.py'));
            cmd=strrep(cmd,"LINK",Pname);
            cmd=strrep(cmd,"COMPONENT_GAP",string(gap));
            cmd=strrep(cmd,"NCIRCLE",string(obj.Ncircle));
            try
                C=SageWrapper.H.exec(cmd);
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
            obj.calcEdgeDirection;
            obj.calcStrandTable;
            TS=obj.StrandTable;
            TE=obj.EdgeTable;
            E=unique(E);

            for si=1:height(TS)
                [~,edges]=ismember(TS.Edges{si},E);
                % reverse edge's segments coordinates array if the ends are not matching
                % or the tangent directions of both ends are not matching
                % vec=diff(TE.Position{edges(end-[0,1])});
                if isscalar(edges)
                    continue
                end
                edges(end+(1:2))=edges(1:2); % close the loop
                for ei=1:length(edges)-1
                    flag=0;
                    segs=PE{edges(ei)};
                    vend=TE.CrossingID(edges(ei),:);
                    if segs(end)~=V(vend(2))
                        flag=1;
                    elseif vend(1)==vend(2)
                        segs=segs(1)+2/3*(segs-segs(1));
                        segs2=PE{edges(ei+1)};
                        z=(segs2(2)-segs2(1))/(segs(end)-segs(end-1));
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
                        PE{edges(ei)}=flip(segs);
                    else
                        PE{edges(ei)}=segs;
                    end
                end
            end
            obj.EdgeTable.Position=PE;
            E=unique(E);
            obj.EdgeTable.ID=E';

            function C = connectPolylines(PE)
                %CONNECTPOLYLINES Reorder and concatenate polylines in cell array PE
                %   C is a single row vector of complex points, with duplicate
                %   connection points removed.
                % https://chatgpt.com/share/686d8571-c200-8003-8380-463bf957ed03

                % Filter out cells with fewer than 2 points
                mask = cellfun(@(v) numel(v) > 1, PE);
                polys = PE(mask);
                n = numel(polys);

                % Precompute start and end of each polyline
                starts = cellfun(@(v) v(1), polys);
                ends   = cellfun(@(v) v(end), polys);
                used = false(1, n);

                % Initialize with the first polyline
                idx0 = 1;
                used(idx0) = true;
                C = polys{idx0};
                current_end = C(end);

                % Iteratively find and append matching polylines
                for k = 2:n
                    found = false;
                    for j = find(~used)
                        if starts(j) == current_end          % matches start
                            next_poly = polys{j};
                        elseif ends(j) == current_end      % matches end, reverse needed
                            next_poly = flip(polys{j});
                            disp("flipped")
                        else
                            continue;
                        end
                        used(j) = true;
                        % Remove duplicate connection point
                        next_poly = next_poly(2:end);
                        % Append
                        C = [C, next_poly];
                        current_end = C(end);
                        found = true;
                        break;
                    end
                    if ~found
                        error('connectPolylines:NoMatch', 'No matching polyline found at step %d.', k);
                    end
                end
                if C(end)==C(end-1)
                    C=C(1:end-1);
                end
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
        function [PDCode,orient,Ncircle]=getPDCode(obj)
            SW=SageWrapper.H;
            obj.setSageLink;
            PDCode=SW.exec(sprintf("%s.pd_code()",obj.sageName));
            PDCode=cellfun(@(x){double(x)},cell(PDCode));
            PDCode=vertcat(PDCode{:});
            orient=double(SW.exec(sprintf("%s.orientation()",obj.sageName)));

        end
        function [GaussCode,orientation]=getGaussCode(obj)
            if obj.formatFlag("G")
                GaussCode=obj.GaussCode;
                orientation=obj.orientation;
            else
                SW=SageWrapper.H;
                obj.setSageLink;
                GaussCode=SW.exec(sprintf("%s.gauss_code()",obj.sageName));
                GaussCode=cellfun(@(x){double(x)},cell(GaussCode));
                orientation=double(SW.exec(sprintf("%s.orientation()",obj.sageName)));
                obj.GaussCode=GaussCode;
                obj.orientation=orientation;
                obj.formatFlag("G")=true;
            end
            if ~isempty(obj.cutFlag)
                % If there are cut edges, append 0 to the end of each Gauss code
                GaussCode(obj.cutFlag)=cellfun(@(x){[x,0]},GaussCode(obj.cutFlag));
            end
        end
        function [GaussCode,orientation]=getVirtualGaussCode(obj)
            if obj.formatFlag("VG")
                GaussCode=obj.virtualGaussCode;
                orientation=obj.virtualOrientation;
            else
                if obj.formatFlag("G")
                    GaussCode=obj.GaussCode;
                    orientation=obj.orientation;
                else
                    [GaussCode,orientation]=obj.getGaussCode();
                end
                vv=find(orientation==0);
                for i = 1:length(GaussCode)
                    idx=ismember(GaussCode{i}, [vv, -vv]);
                    GaussCode{i} = GaussCode{i}(~idx); % Remove virtual crossings
                end
                orientation(vv)=[];
                obj.virtualGaussCode=GaussCode;
                obj.virtualOrientation=orientation;
                obj.formatFlag("VG")=true;
            end
            if ~isempty(obj.cutFlag)
                % If there are cut edges, append 0 to the end of each Gauss code
                GaussCode(obj.cutFlag)=cellfun(@(x){[x,0]},GaussCode(obj.cutFlag));
            end
        end
        function calcVG2G(obj)
            vgc=obj.virtualGaussCode;
            vsgn=obj.orientation;
            [gc,sgn]=vg2g(vgc,vsgn);
            obj.GaussCode=gc;
            obj.orientation=sgn;
            obj.formatFlag("G")=true;
        end
        function calcStrandTable(obj)
            % Calculate the strand table based on edge connections

            [gc,ori_]=getGaussCode(obj);
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
            obj.StrandTable.Vertices =gc.'; % Assign the vertices to the table
        end
        function [hi,ti]=calcEdgeDirection(obj)
            % Calculate edge directions based on PD code on Sage
            pname=obj.sageName;
            cmd=[...
                "t,h="+pname+"._directions_of_edges()";
                "pd="+pname+".pd_code()";
                "ti = [1 + pd.index(t[e])  for e in sorted(t.keys())]";
                "hi = [1 + pd.index(h[e])  for e in sorted(h.keys())]"
                "ti,hi"];
            ret=SageWrapper.H.exec(cmd);
            hi=double(ret{1});
            ti=double(ret{2});
            obj.EdgeTable.CrossingID=[hi;ti].';
        end


    end
    %% plot, display, and disp methods
    methods
        function p=plot(obj,arg)
            arguments
                obj
                arg.cut (1,:) = [] % cut specified edge indices of strands
                arg.color  = 'b' % color of the knot
                arg.lineWidth (1,1) double = 2 % line width for the knot
                arg.markerSize (1,1) double = 100 % size of markers for crossings
                arg.edgeLabel {mustBeMember(arg.edgeLabel,["none","ID","weight","all"])} = "all" % edge label type
                arg.vertexLabel {mustBeMember(arg.vertexLabel,["none","ID","ori","all"])} = "all" % vertex label type
                arg.vertexStyle {mustBeMember(arg.vertexStyle,["none","real","imag","all"])} = "real" % vertex dot style
                arg.gap (1,1) double = 0.2 % gap between edges
                arg.componentGap (1,1) double = 1 % gap between components
            end
            % Plot the VirtualLink diagram
            [PE,V]=obj.calcPositions;
            PD=obj.PDCode;
            orientation=obj.orientation;
            obj.calcEdgeDirection;
            TE=obj.EdgeTable;
            obj.calcStrandTable;

            gap=0.2; % gap between edges
            gappos=[nan,nan];
            TS=obj.StrandTable;
            if isempty(arg.cut)
                arg.cut=false(1,height(TS));
            end
            assert(length(arg.cut)==height(TS))
            arg.cut(arg.cut==Inf)=cellfun(@length,TS.Edges(arg.cut==Inf));

            figure
            hold on
            for si=1:height(TS) % over strands
                edgeList=TS.Edges{si};
                sgc=sign(TS.Vertices{si});
                short=and(sgc<0,orientation(abs(TS.Vertices{si})));
                % short=[short([end,1:end-1]);short];%1ずれてる可能性あり
                short=[short;short([2:end,1])];%1ずれてる可能性あり
                strand=[];
                for ei=1:length(edgeList) % over edges
                    edge=edgeList(ei);
                    segs=TE.Position{edge};
                    % calculate middle point of the segment to place the label
                    if ei~=arg.cut(si)
                        seg=segs(floor(end/2)+(0:1));
                        mp=mean(seg);
                        arrowhead=uvec(seg(1:2))*[-1+.5i,0,-1-.5i]/7+mp;
                        plot(arrowhead,LineWidth=2,Color='b')
                        % Place the label at the middle point
                        text(real(mp),imag(mp),sprintf('E%d',edge),'VerticalAlignment','top');
                    else
                        segs=[segs(1)+[0,1/3]*uvec(segs(1:2)) nan segs(end)-[1/3,0]*uvec(segs(end-1:end))];
                    end
                    % short=~mod(find(PD(TE.CrossingID(edge,:),:).'==edge),2);
                    segs(1)=segs(1)+gap*short(1,ei)*uvec(segs(1:2));
                    segs(end)=segs(end)+gap*short(2,ei)*uvec(segs([end,end-1]));
                    strand=[strand,gappos(short(1,ei)),segs,gappos(short(2,ei))];
                end

                plot(strand,LineWidth=arg.lineWidth,Color=arg.color)
            end
            isVirtual=logical(orientation);
            Vr=V(isVirtual);
            scatter(real(Vr),imag(Vr),arg.markerSize,arg.color,"filled")
            text(real(Vr),imag(Vr),"V"+string(1:length(Vr)),'VerticalAlignment','bottom')
            xlim padded
            ylim padded
            axis equal
            hold off

            set(gca, 'LooseInset', max(get(gca, 'TightInset'), 0))

            function ret=uvec(v)
                ret=(v(2)-v(1))/abs(v(2)-v(1));
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
        function disp1(obj)
            % Short text display
            disp('<a href="matlab:help VirtualLink">Virtual Link object</a>:')
            if obj.formatFlag("G")
                fprintf("Gauss Code:")
                dispC(obj.GaussCode);
                fprintf("Orientation:")
                disp("["+join(string(obj.orientation),", ") + "]")
            elseif obj.formatFlag("VG")
                fprintf("Virtual Gauss Code:")
                dispC(obj.virtualGaussCode);
                fprintf("Orientation:")
                disp("["+join(string(obj.virtualOrientation),", ") + "]")
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
    end
    %% setData
    methods
        function setData(obj,arg)
            % SetData Set data for the VirtualLink object.
            % Input can be (virtual)GaussCode, (virtual)PDCode, or ostring.
            % virtual data means that extra virtual crossings have to be added to be a (virtual) link.
            %
            arguments
                obj
                arg.orientation
                arg.GaussCode
                arg.virtualGaussCode
                arg.PDCode
                arg.Ncircle
                arg.isVirtual
                arg.ostring
                arg.DTCode
                arg.table
            end
            obj.formatFlag=obj.formatFlag0;
            if isfield(arg,"GaussCode")
                assert(all(isfield(arg,["orientation"])))
                assert(isa(arg.GaussCode,'cell'))
                obj.GaussCode=arg.GaussCode;
                obj.orientation=arg.orientation;
                obj.formatFlag("G")=true;
                NC=sum(cellfun(@length,arg.GaussCode)==0);
                NV=length(arg.orientation);
                obj.virtualFlag=false;
            elseif isfield(arg,"virtualGaussCode")
                assert(all(isfield(arg,["orientation"])))
                obj.virtualGaussCode =helperCut( arg.virtualGaussCode,obj);
                % obj.orientation=arg.orientation;
                obj.virtualOrientation=arg.orientation;
                obj.formatFlag("VG")=true;
                % obj.calcVG2G();
                NC=sum(cellfun(@length,arg.virtualGaussCode)==0);
                NV=length(arg.orientation);
                obj.virtualFlag=true;
            elseif isfield(arg,"PDCode")
                % assert(all(isfield(arg,["isVirtual"])))
                if isfield(arg,["Ncircle"])
                    NC=arg.Ncircle;
                elseif ~isempty(arg.PDCode)
                    NC=0; % if not specified, assume that no unknots are included
                else
                    NC=1; % unknot
                end
                obj.formatFlag("PD")=true;
                obj.PDCode=arg.PDCode;
                NV=size(arg.PDCode,1);
                obj.orientation=ones(1,NV);
                if isfield(arg,"isVirtual")
                    obj.isVirtual=arg.isVirtual;
                else
                    obj.isVirtual=false(1,NV);
                end
                obj.virtualFlag=false;
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
                obj.setData(virtualGauss=gc,orientation=ori);
                return
            end
            if ~obj.virtualFlag
                obj.setSageLink()
                if obj.formatFlag("PD")
                    [~,ori]=obj.getPDCode;
                    obj.sageOrientation=ori;
                end
            end
            obj.Ncircle=NC;
            obj.CrossingTable=cell2table(num2cell(nan(NV,length(obj.CrossingTableVariableNames))), ...
                VariableNames=obj.CrossingTableVariableNames);
            obj.EdgeTable=cell2table(num2cell(nan(NV*2+NC,length(obj.EdgeTableVariableNames))), ...
                VariableNames=obj.EdgeTableVariableNames);

            function gc=helperCut(gc,vlobj)
                % if 0 is inserted in the Gauss code, it is a cut edge
                % and gc code is stored after removing 0
                % cut flag is stored in obj.cutFlag
                [~,cutIdx]=cellfun(@(x)ismember(0,x),gc);
                if any(cutIdx)
                    % If there are cut edges, append 0 to the end of each Gauss code
                    cutFlag=logical(cutIdx);
                    gc(cutFlag)=arrayfun(@(x,i){x{1}([i+1:end,1:i-1])},gc(cutFlag),cutIdx(cutFlag));
                    vlobj.cutFlag=cutFlag;
                end
            end
        end
    end
end


function dispC(arg)
    if isempty(arg)
        disp("[[]]")
    else
        disp("[["+join(cellfun(@(x)join(string(x),","),arg),"], [") + "]]")
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
function [gc,sgn]=vg2g(vgc,vsgn)


    if isempty(vsgn)
        assert(~any(cellfun(@numel,vgc)))
        gc=vgc;
        sgn=vsgn;
        return
    end
    % store information about the regions of the virtual graph
    % using vertex/edge/face (regV/regE/regF)
    % regF is a graph object representing the adjacency of regions
    % and its Weights are the edge IDs connecting the regions
    regV={};
    regE={};
    gc=vgc;
    ge=cell(size(gc));
    sgn=vsgn;
    usedV=0;
    cntE=0;
    cntV=max(cellfun(@max,vgc))+1;
    regF=graph;
    for si=1:length(vgc)
        % loop through each strand
        cgc=[vgc{si} 0];
        if isscalar(cgc)
            vgc{si}=[];
            continue
        end
        CU
        CI
        v1=abs(cgc(1));
        if ismember(v1,usedV)
            % start strand sweeping from a used vertex
            [f2,vidx2]=findC(-v1,regV);
            regV{f2}=[cntV+1,regV{f2}([vidx2:end,1:vidx2])];
            regE{f2}=[-cntE,regE{f2}(1:vidx2-1),regE{f2}(vidx2:end),cntE];
            regF=regF.addedge(f2,[f2 f2],[0 cntE]);
            [f1,vidx]=findC(v1,regV);
            regV{f1}=[regV{f1}(1:vidx-1),-v1,0,-v1,regV{f1}(vidx+1:end)];
            regE{f1}=[regE{f1}(1:vidx-1),0,0,regE{f1}(vidx:end)];
            % regF=regF.addedge(f1,f1,0);
            f1=f2;
        else
            % start strand sweeping from a new vertex
            f1=length(regV)+1;
            regV{f1}=[0,-1,0,1]*vsgn(v1)*cgc(1)+[cntV+1,0,0,0];%?
            regE{f1}=[-cntE,0,0,cntE];
            regF=regF.addedge(f1,[f1 f1],[0 cntE]);
            usedV(end+1)=v1;
        end
        for vi=2:length(cgc)
            % loop through each vertex in the strand
            fd
            CU
            CI
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
                    regE{f2}=[regE{f2}(1:vidx2-1),regE{f1},regE{f2}(vidx2+1:end)];
                    regE{fnew}=[-cntE,regE{fnew}([vidx3:end,1:vidx3-1]),cntE];
                    % contract the faces f1<->f2
                    ns=regF.neighbors(f1);
                    for n=ns(ns~=f1)'
                        bnd=regF.Edges.Weight(regF.findedge(f1,n));
                        regF.addedge(f2,n*ones(size(bnd)),bnd);
                    end
                    bnd=regF.Edges.Weight(regF.findedge(f1,f1));
                    regF.addedge(f2,f2*ones(size(bnd)),bnd);
                    regF.rmnode(f1);
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
                            CU
                            bnds=regF.Edges.Weight(findedge(regF,f1,f2));
                            bnd=bnds(1); %ここが選択の余地あり
                            vidx=find(bnd==abs(regE{f1}),1)+1;
                            eo=regE{f1}(vidx-1)>0;
                            regV{f1}=[regV{f1}(1:vidx-1),cntV,regV{f1}(vidx:end)];
                            regE{f1}=[regE{f1}(1:vidx-1-~eo),(2*eo-1)*cntE,regE{f1}(vidx-1+eo:end)];
                            vidx2=find(bnd==abs(regE{f2}),1)+1;
                            regV{f2}=[regV{f2}(1:vidx2-1),-cntV,regV{f2}(vidx2:end)];
                            regE{f2}=[regE{f2}(1:vidx2-1-eo),-(2*eo-1)*cntE,regE{f2}(vidx2-1+~eo:end)];
                            regF=regF.addedge(f1,f2,cntE);
                            vs=cntV;
                            [sidx,geidx]=findC(bnd,ge);
                            ge{sidx}=[ge{sidx}(1:geidx), cntE, ge{sidx}(geidx+1:end)];
                            gc{sidx}=[gc{sidx}(1:geidx), vs, gc{sidx}(geidx+1:end)];
                            [~,gcidx]=ismember(abs(cgc(vi)),[abs(gc{si}) 0]);
                            gc{si}=[gc{si}(1:gcidx-1),-vs,gc{si}(gcidx:end)];
                            sgn(cntV)=0;
                            v1=cntV;
                            cntV=cntV+1;
                        elseif vi==length(cgc)
                            % connect ends
                            break
                        else
                            v1=abs(cgc(vi));
                            vidx=vidx2;
                            [f2,vidx2]=findC(-v1,regV);
                        end
                        CU
                        fnew=length(regV)+1;
                        % pattern 4: cut through edges and add virtual crossing
                        regV{fnew}=[-v1 regV{f1}(vidx+1:end)];
                        regV{f1}=[-v1 regV{f1}(2:vidx-1)];
                        regV{f2}=[cntV+1 regV{f2}([vidx2:end,1:vidx2])];
                        regE{fnew}=regE{f1}(vidx:end);
                        regE{f1}=regE{f1}(1:vidx-1);
                        regE{f2}=[-cntE regE{f2}([vidx2:end,1:vidx2-1]) cntE];
                        regF=updateFace(regF,abs(regE{fnew}),f1,f2,fnew,cntE);
                        f1=f2;
                        CI
                    end
                end
            else
                % pattern 5: add new crossing to temporary graph
                usedV(end+1)=v1;
                vs=vsgn(v1)*cgc(vi);
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
        fd
        % regF.plot
        dispC(ge)
        dispC(gc)
    end


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
        try
            fprintf("f1=%d\n",f1)
            fprintf("si=%d, vi=%d, v_gc=%d, s_gc=%d\n",si,vi,cgc(vi),vsgn(abs(cgc(vi))))
        catch
        end
        dispC(ge)
        dispC(regV)
        dispC(regE)
        disp(regF.Edges{:,:}')

    end
    function CU()
        cntE=cntE+1;
    end
    function CI()
        ge{si}(end+1)=cntE;
    end
end
