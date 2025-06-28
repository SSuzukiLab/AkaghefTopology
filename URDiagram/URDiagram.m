classdef URDiagram<handle
    % URDiagram - A class to represent and manipulate UR Diagrams.

    properties
        C (1,1) sym=1
        V (1,:) double % V must contain ±pair of indices of edges
        E %(prop,chord idx)
        W (1,:) sym
        section
    end
    properties(Dependent)
        tbl
        rank
    end
    properties(Constant,Hidden)
        tblname=["chord","notdirty","dilation","dilinv","empty2"];
    end

    methods
        %% Relation
        function ret=nextIdx(obj)
            ret=max(obj.V)+1;
        end
        function [ret,idx]=getW(obj,arg)
            [~,idx]=ismember(abs(arg),obj.E(1,:));
            ret=obj.W(idx);
        end
        function ret=getEidx(obj,arg)
            [~,ret]=ismember(arg,obj.E(1,:));
        end
        function ret=getVidx(obj,arg)
            [~,ret]=ismember(arg,obj.V);
        end
        function ret=getTidx(obj,arg)
            idx=obj.getVidx(arg);
            ret=arg;
            for i=1:numel(arg)
                ret(i)=sum(obj.V(1:idx(i))==0)+1;
            end
        end
        function obj=put0(obj,argP,argN)
            newIdx=obj.nextIdx();
            V=[obj.V newIdx -newIdx];
            if argP<=argN
                obj.V=V([1:argP end-1 argP+1:argN end argN+1:end-2]);
            else
                obj.V=V([1:argN end-1 argN+1:argP end argP+1:end-2]);
            end
            obj.E(1,end+1)=newIdx;
            obj.W(end+1)=0;

        end
        function swap(obj,arg,skip)
            arguments
                obj
                arg
                skip  =false
            end

            V=obj.V;
            pos=obj.getVidx(arg);
            [~,pos]=ismember(arg,V);
            if ~skip
                assert(numel(arg)==2,"must specify 2adjacent indices")
                if any(pos==0)
                    error('Invalid: argument must be existing vertices.');
                elseif abs(pos(1)-pos(2))~=1
                    error('Invalid: the vertices must be adjacent.');
                end
            end
            sgn=sign(V(pos));
            V(pos)=V(pos([2 1]));
            if arg(1)==-arg(2)
                % dilation swap
                [W,Widx]=obj.getW(arg(1));
                if ~xor(pos(1)<pos(2),sgn(1)<0)
                    % ∂x->x∂
                    obj.W(Widx)=simplify(W/(1-W));
                    obj.C=simplify(obj.C/(1-W));
                else
                    % x∂->∂x
                    obj.W(Widx)=simplify(W/(1+W));
                    obj.C=simplify(obj.C/(1+W));
                end
            elseif sgn(1)==sgn(2)
                % swap, H×H, H^*×H^*
            else
                % swap, H×H^*, H×H^*
                weight=obj.getW(arg);
                [~,pos2]=ismember(-arg,V);
                if xor(sgn(1)>0,pos(2)<pos(1))
                    newWeight=-weight(1)*weight(2);
                else
                    newWeight=weight(1)*weight(2);
                end
                newIdx=obj.nextIdx();
                if pos2(1)<pos2(2)
                    V=[V(1:pos2(1)-1) -sgn(1)*newIdx V(pos2(1):pos2(2)-1) -sgn(2)*newIdx V(pos2(2):end)];
                else
                    V=[V(1:pos2(2)-1) -sgn(2)*newIdx V(pos2(2):pos2(1)-1) -sgn(1)*newIdx V(pos2(1):end)];
                end
                obj.W(end+1)=newWeight;
                obj.E(1,end+1)=newIdx;
            end
            obj.V=V;
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
            try
            catch
                error('Invalid index for swapping vertices.');
            end
        end
        function add(obj,arg,skip)
            arguments
                obj
                arg (1,2) double
                skip=false
            end
            V=obj.V;
            arg=abs(arg);
            [~,posPlus]=ismember(arg,V);
            [~,posMinus]=ismember(-arg,V);
            if ~skip
                if any(posPlus==0)||any(posMinus==0)
                    error('Invalid: argument must be existing vertices.');
                elseif abs(posPlus(1)-posPlus(2))~=1||abs(posMinus(1)-posMinus(2))~=1
                    error('Invalid: the vertices must be adjacent.');
                end
            end
            eidx=obj.getEidx(arg);
            obj.V([posPlus(2) posMinus(2)])=[];
            obj.W(eidx(1))=obj.W(eidx(1))+obj.W(eidx(2));
            obj.W(eidx(2))=[];
            obj.E(:,eidx(2))=[];
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end

        end
        function obj=dilation(obj,dil,arg,skip)
            arguments
                obj
                dil (1,1) double
                arg (1,1) double
                skip=false
            end
            dil=abs(dil);
            pos=obj.getVidx([dil -dil arg]);
            [dist,sgn1]=min(abs(pos(1:2)-pos(3)));
            e1=pos(2)-pos(1);
            e2=sign(arg);
            e3=sign(pos(3)-pos(1));
            idx1=obj.getEidx(dil);
            idx2=obj.getEidx(abs(arg));
            if ~skip
                if any(pos==0)
                    error('Invalid: specify the edge indices')
                elseif ~isempty(setxor(pos-min(pos),0:2))
                    error('Invalid: must be adjacent indices')
                elseif 0
                    error('Invalid: 2nd arguments must be a dilation operator')
                end
            end
            if e1<0
                mul=obj.W(idx1)/(1-obj.W(idx1));
            else
                mul=obj.W(idx1)+1;
            end
            obj.W(idx2)=obj.W(idx2)*mul^(e2*e3);
            obj.V(pos+e3*[1 1 -2])=obj.V(pos);
        end
        function obj=composeDilation(obj,arg,skip)
            arguments
                obj
                arg (1,2) double
                skip=false
            end
            arg=abs(arg);
            pos=obj.getVidx([arg -arg]);
            if ~skip
                if any(pos==0)
                    error('Invalid: argument must be existing vertices.');
                elseif abs(pos(1)-pos(3))~=1||abs(pos(2)-pos(4))~=1
                    error('Invalid: the vertices must be dilation indices.');
                elseif ~isempty(setxor(pos-min(pos),0:3))
                    error('Invalid: must be adjacent dilation indices.');
                end
            end
            if pos(1)>pos(3)
                % obj.V()を追加する必要がある
                obj.swap(obj.V(pos([1 3])),skip);
            elseif pos(2)>pos(4)
                obj.swap(obj.V(pos([2 4])),skip);
            end

            [w,idx]=obj.getW(arg);
            obj.W(idx(1))=(1+w(1))*(1+w(2))-1;
            obj.E(:,idx(2))=[];
            obj.W(idx(2))=[];
            obj.V(pos([2 4]))=[];
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end

        end
        %% reduction
        function obj=updateE(obj)
            % dirty=~obj.E(2,:);
            dirty=ones(size(obj.E(2,:)));
            % 変更箇所のみを更新する改善がある
            for i=1:length(dirty)
                if ~dirty(i),continue;end
                % ["chord","notdirty","dilation","posdif","dilinv"]
                n=obj.E(1,i);
                pos=obj.getVidx([n -n]);
                obj.E(3,i)
                posdif=abs(pos(2)-pos(1));
                obj.E(3,i)=posdif==1;
                obj.E(4,i)=posdif;
                obj.E(5,i)=posdif==-1;
                obj.E(2,dirty)=1;
            end
        end
        function reduct1(obj)
            % dilationを掃き出す
            n=obj.E(1,:);
            N=length(n);
            pos=obj.getEidx([n -n]);
            pd=pos(1:N)-pos(N+1:end);
            % iD=find(abs(pd)==1);
            % if isempty(iD)
            %     return
            % end
            isdil=arrayfun(@isadj,pos(1:N),pos(N+1:end));
            for i=1:N
                if ~isdil(i)
                    continue
                end

            end
            function ret=isadj(arg1,arg2)
                [minpos,maxpos]=bounds(arg1,arg2);
                arr=obj.V(minpos:maxpos);
                ret=sum(abs(diff(sign(arr)))==2)==1;
            end
        end
        function ret=reduct2(obj)
            % 転倒を解消する
            idx=1;
            idx1=0;
            e1=false;
            cnt=0;
            while idx<=length(obj.V)&&cnt<50
                if idx<=0
                    idx=1;
                end
                switch sign(obj.V(idx))
                    case 1
                        if e1
                            arr=idx1:idx-1;
                            for k=arr
                                obj.swap(obj.V([k k+1]));
                            end
                            idx=idx-2;
                        end
                    case -1
                        if ~e1
                            idx1=idx;
                            e1=true;
                        else
                            idx1=idx;
                        end
                    otherwise
                        e1=false;
                end
                idx=idx+1;
                cnt=cnt+1;
            end
        end
        function [gp,gn]=getGroup(obj)
            n=obj.E(1,:);
            P=[0,find(diff(sign(obj.V)))].';
            posP=obj.getVidx(n);
            gp=sum(posP>P);
            posN=obj.getVidx(-n);
            gn=sum(posN>P);
        end
        function reduct3(obj)
            % addできる部分を簡約する
            N=length(obj.E(1,:));
            [gp,gn]=obj.getGroup();
            m=gp-gp.'==0&gn-gn.'==0;
            m=triu(m,1);
            d=[];
            cidx=zeros(0,2);
            for i=find(m).'
                [u,v]=ind2sub([N N],i);
                if ismember(v,d)
                    continue
                end
                cidx=[cidx;obj.E(1,[u,v])];
                d=[d,v];
            end
            for j=1:size(cidx,1)
                obj.add(cidx(j,1:2),true);
            end
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
        end
        function reduct4(obj,dilIdx)
            % dilationを掃き出す dilIdxは最後尾のdilationのEidx
            [gp,gn]=obj.getGroup();
            e1=gp-gn==[1;-1];
            e3=obj.getEidx(dilIdx);
            e2=setdiff(find(e1(1,:)|e1(2,:)),e3);
            e4=zeros(0,2);
            for i=1:length(e2)
                j=obj.E(1,e2(i));
                pos=obj.getVidx([j -j]);
                if e1(1,e2(i))
                    arr=obj.V(pos(2)+1:pos(1)-1);
                else
                    arr=obj.V(pos(1)+1:pos(2)-1);
                end
                if pos(1)>pos(2)
                    arr=fliplr(arr);
                end
                for k=arr(arr>0)
                    obj.swap([j,k]);
                end
                for k=fliplr(arr(arr<0))
                    obj.swap([-j,k]);
                end
                if e1(1,e2(i))
                    obj.swap([j -j]);
                end
                idx=obj.getTidx(j);
                pos2=obj.getVidx([-j dilIdx(idx)]);
                arr=obj.V(pos2(1)+1:pos2(2)-1);
                for k=arr
                    obj.dilation(j,k);
                end
                e4(end+1,1:2)=[dilIdx(idx) j];
            end
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
            for i=size(e4,1):-1:1
                obj.composeDilation(e4(i,1:2));
            end

        end
        function reduct5(obj)
            % reduct4による簡約化を行う
            e1=length(obj.V);
            obj.put0(e1,e1);
            e2=obj.V(end-1);
            obj.reduct4(e2);
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
        end
        function reduct6(obj)
            % 転倒を1つ解消する
            e1=false;
            for idx=1:length(obj.V)
                switch sign(obj.V(idx))
                    case 1
                        if e1
                            obj.swap(obj.V([idx idx-1]));
                            return
                        end
                    case -1
                        e1=true;
                end
            end
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
        end
        function [obj,isComplete]=reduct7(obj)
            % 各階の転倒を1つ解消する
            e1=false;
            e2=false;
            isComplete=false(1,0);
            for idx=1:length(obj.V)
                switch sign(obj.V(idx))
                    case 1
                        if e1&~e2
                            obj.swap(obj.V([idx idx-1]));
                            e2=true;
                        end
                    case -1
                        e1=true;
                    case 0
                        isComplete(end+1)=~e2;
                        e1=false;
                        e2=false;
                end
            end
            isComplete(end+1)=~e2;
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
        end

        function reduct8(obj)
            % reduct4(dilation掃き出し)による簡約化を行う
            pos0=[find(obj.V==0)-1 length(obj.V)];
            for i=1:length(pos0)
                pos0(i)=pos0(i)+2*(i-1);
                p=pos0(i);
                obj.put0(p,p);
            end
            e2=obj.V(pos0+1);
            obj.reduct4(e2);
            pos1=[-2,obj.getVidx(e2)];
            V=obj.V;
            for i=1:length(e2)
                arr=V(pos1(i+1)-1:-1:pos1(i)+3);
                for k=arr
                    if k>0
                        break
                    else
                        obj.dilation(e2(i),k);
                    end
                end
            end
            if TopologyConfig.H.URupdateDisp
                disp(obj)
            end
        end
        function ret=isBasic(obj)
            e1=false;
            ret=false;
            for idx=1:length(obj.V)
                switch sign(obj.V(idx))
                    case 1
                        if e1, return; end
                    case -1
                        e1=true;
                    case 0
                        e1=false;
                end
            end
            n=obj.E(1,:);
            T=obj.getTidx([n;-n].');
            if numel(unique(T,"rows"))<numel(T)
                return
            end
            ret=true;
        end
        function simplifySim(obj)
            obj.W=simplify(obj.W);
            obj.C=simplify(obj.C);
        end
        function simplify1(obj)
            % 1階の簡約化
            cnt=0;
            while numel(obj.V)>2&&cnt<50
                cnt=cnt+1;
                if TopologyConfig.H.URupdateDisp
                    disp("iter="+cnt)
                end
                % addできる部分を簡約する
                obj.reduct3();
                % dilationを掃き出す
                obj.reduct5();
                % 転倒を1つ解消する
                obj.reduct6();
                % obj.validateData
            end
        end
        function simplify2(obj)
            % N階の簡約化
            cnt=0;
            while cnt<50
                cnt=cnt+1;
                if mod(cnt,5)==0
                    obj.simplifySim;
                end
                if TopologyConfig.H.URupdateDisp
                    disp("iter="+cnt)
                end
                % addできる部分を簡約する
                obj.reduct3();
                % dilationを掃き出す
                obj.reduct8();
                % 転倒を1つ解消する
                [~,isComplete]=obj.reduct7();
                if all(isComplete)
                    break
                end
                % obj.validateData
            end
            obj.reduct3();
            obj.simplifySim;
        end
        %% invariants
        function ret=trace(obj)
            obj.simplify2;
            ret=obj.trace2;
        end
        function ret=trace1(obj)
            assert(isscalar(numel(obj.W)),"cannot calculate trace")
            e=simplify(-obj.C/obj.W);
            vs=symvar(e);
            vs=vs(contains(string(vs),"eps"));
            for v=vs
                e=limit(e,v,0);
            end
            ret=e;
        end
        function M=getMatrixDisplay(obj)
            assert(obj.isBasic,"cannot determine trace\n must be simplified")
            n=obj.E(1,:);
            T=obj.getTidx([n;-n].');
            M=zeros(obj.rank,'sym');
            for i=1:size(T,1)
                M(T(i,1),T(i,2))=obj.W(i);
            end

        end
        function ret=trace2(obj)
            M=obj.getMatrixDisplay;
            % 符号はテキトーなので変更の可能性あり
            detM=det(-M);
            e=simplify(obj.C/detM);
            vs=symvar(e);
            vs=vs(contains(string(vs),"eps"));
            for v=vs
                e=limit(e,v,0);
            end
            ret=e;
        end



        %% Construction Validation
        function obj = URDiagram()
            % Initialize the URDiagram object
        end
        function obj=setDataFromCVW(obj,C, V, W)
            % Set the vertex and edge data
            obj.C=C;
            obj.V = V;
            E=unique(V(V>0));

            obj.E = E;
            obj.E(2:5,:)=0;
            obj.W=W;
            obj.section=[0 find(V)];
        end
        function obj=setC(obj,c)
            obj.C=c;
        end
        function obj=setDataFromCVS(obj,C, V, S,eps)
            % Set the vertex and edge data
            obj.setDataFromCVW(C,V, S);
            if  nargin==4
                eps=str2sym("eps"+obj.E(1,:));
            end
            obj.W=obj.W+eps;
        end
        function obj=setDataFromV(obj, V)
            % Set the vertex
            obj.setDataFromCVW(1,V,nan*unique(V(V>0)));
        end
        function obj=setDataFromOData(obj,str)
            % 入力: O-data文字列 'v_{1}^{+r}v_{2}^{-l}⊗v_{3}^{+r}'
            % 出力: obj.V, obj.W を復元
            V = []; W = [];
            pattern = 'v_\{(\d+)\}\^\{([+-][rl])\}|⊗';
            tokens = regexp(str, pattern, 'tokens');

            for i = 1:length(tokens)
                token = tokens{i};
                if isempty(token)
                    V(end+1) = 0;
                else
                    idx = str2double(token{1});
                    sign = token{2}(1);
                    dir = token{2}(2);
                    if sign == '+'
                        V(end+1) = idx;
                    else
                        V(end+1) = -idx;
                    end
                    if dir == 'r'
                        W(idx) = +1;
                    else
                        W(idx) = -1;
                    end
                end
            end
            V=fliplr(V);
            obj.setDataFromCVS(1,V,W);
        end


        function validateData(obj)
            % Validate the vertex and edge data
            V=obj.V;
            Vabs=unique(abs(V));
            % if ~isequal([Vabs Vabs],[Vplus Vminus])
            if sum(V>0)~=length(Vabs)||sum(V<0)~=length(Vabs)

                Vplus=V(V>0);
                Vminus=-V(V<0);
                assert(numel(Vplus)==numel(unique(V(V>0))));
                assert(numel(Vminus)==numel(unique(-V(V<0))));
                Vnonpair=[setxor(Vplus,Vabs) -setxor(Vminus,Vabs)];
                if ~isempty(Vnonpair)
                    error('The following vertices are not paired: %s', mat2str(Vnonpair));
                end
            end
            assert(isequal(Vabs,sort(obj.E(1,:))),"error")

            % if isempty(obj.V) || isempty(obj.E)
            %     error('Vertex or edge data is empty.');
            % end
            % Add more validation checks as needed
        end
        function checkData(obj,C,V,W)
            try
                assert(all(obj.C==C),"C")
                assert(all(obj.V==V),"V")
                assert(all(obj.W==W),"W")
            catch ME
                disp(obj.C)
                disp(obj.V)
                disp(obj.W)
                rethrow(ME)
            end
        end


        %% Display method
        function ret=get.tbl(obj)
            ret=cell2table(num2cell(obj.E.'),Var=URDiagram.tblname);

            W=obj.W.';
            if isempty(symvar(sum(W)))&&~any(isnan(W))
                [n,d]=rat(W,1e-4);
                W=sym(n)./sym(d);
            end
            ret.weight=W;
        end
        function ret=get.rank(obj)
            ret=sum(obj.V==0)+1;
        end
        function ret=getOdata(obj)
            W=obj.lim0(obj.W);
            assert(all(abs(W)==1))
            ret='';
            for i=1:length(obj.V)
                V=obj.V(end+1-i);
                if V==0
                    ret(end+1)='⊗';
                else
                    eidx=obj.getEidx(abs(V));
                    if V>0, str2='+';
                    else,   str2='-'; end
                    if W(eidx)>0,   str2(2)='r';
                    else,           str2(2)='l'; end
                    str=sprintf('v_{%s}^{%s}', ...
                        string(abs(V)),str2);
                    ret=[ret str];
                end
            end
        end
        function ret=getWExpression(obj)
            W=obj.lim0(obj.W);
            assert(all(abs(W)==1))
            ret="";
            for i=1:length(obj.V)
                V=obj.V(i);
                if V==0
                    ret(end+1)='⊗';
                else
                    eidx=obj.getEidx(abs(V));
                    str3=string(abs(V));
                    if V>0, str3=str3+'2';
                    else,   str3=str3+'1'; end
                    
                    if W(eidx)>0,   str2="";
                    else,           str2="\overline"; end
                    
                    str=str2+"{\WO}_{"+str3+"}^{\epsilon_"+abs(V)+"}";
                    ret=ret+str;
                end
            end
        end

        function disp(obj)
            feval("disp"+ TopologyConfig.H.URDdisp,obj);
        end
        function disp0(obj)
            % Default display
            builtin('disp',obj);
        end
        function disp1(obj)
            % Short text display
            disp(obj.V);
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
        function plot(obj)
            figure
            L=length(obj.V)+1;
            quiver(L,0,-L,0,AutoScaleFactor=1,LineWidth=2,Color="b");
            hold on
            W=obj.W;
        end
        function ret=unicodeArt(obj,replace)
            V=obj.V;
            edge=obj.E(1,:);
            hei=length(edge);
            wid=length(V);
            d=dictionary([edge 0],[1:hei nan]);
            hpos=d(abs(V));
            str=char(zeros(hei+1,wid));
            str(end,:)="│";
            str(end,find(V==0))="⊗";
            for i=1:hei
                e1=false;
                for j=1:wid
                    if hpos(j)==i
                        if e1
                            e1=false;
                            str(end-i,j)="╮";
                        else
                            e1=true;
                            str(end-i,j)="╭";
                        end
                    else
                        if str(end-i+1,j)=="│"
                            str(end-i,j)="│";
                        elseif e1
                            str(end-i,j)="─";
                        else
                            str(end-i,j)=" ";
                        end
                    end
                end
                %  ╭─╮
                %  │╭│╮
                % ╭─│╯│
                % │╰╯ │
                % ╰───╯
            end
            for j=1:wid
                if V(j)>0
                    str(end,j)="+";
                elseif V(j)<0
                    str(end,j)="-";
                end
            end
            % var=symvar(obj.W);
            % var=var(contains(string(var),"eps"));
            if nargin==2
                W0=obj.W;
                C0=obj.C;
            else
                % W0=subs(obj.W,var,zeros(size(var)));
                % C0=subs(obj.C,var,zeros(size(var)));
                W0=obj.lim0(obj.W);
                C0=obj.lim0(obj.C);
                % W0(isnan(W0))=Inf;
            end

            weight=[fliplr(W0),C0].';
            diagram=categorical("･"+string(str));
            ret=table(weight,diagram, ...
                VariableNames=["weight","diagram"], ...
                RowNames=[string(edge) "C"]);
        end

        % Define other methods here

    end
    methods(Static)
        function ret=lim0(arg)
            vs=symvar(arg);
            vs=vs(contains(string(vs),"eps"));
            ret=arg;
            for v=vs
                ret=limit(ret,v,0);
            end
            ret(isnan(ret))=Inf;
        end

    end
end