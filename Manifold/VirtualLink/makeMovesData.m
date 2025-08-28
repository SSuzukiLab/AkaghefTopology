% This code generates move data for virtual links
% For the graph representation of virtual links, see Muramatsu-Suzuki-Taguchi(2024).
% gc:gauss code, gcFirst: matching prefix of gauss code, ori: orientation, 
% NV: number of vertices, param: parameter
% NOT SUPPORTED: PS move, CP move, moves with weight, H-move, weight move
cd(fullfile(TopologyConfig.ProjectPath,'Manifold','VirtualLink'))
%% MP move
T_MP_L=table(["A";"B";"C";"D"], ...
   [{[-1],[1,2],[-2]};
   {[1],[-1,-2],[2]};
   {[-1,2],[1],[-2]};
   {[2,-1],[1],[-2]}], ...
   [[1,1];[-1,-1];[-1,1];[1,-1]],VariableNames=["label","gc","ori"]);
T_MP_L=repmat(T_MP_L,4,1);
T_MP_L.num=repelem(1:4,1,4)';
T_MP_L.rev=repelem([0,0,0;0,0,1;0,1,0;0,1,1],4,1);
T_MP_L.param=T_MP_L.label+T_MP_L.num;
tmp=rowfun(@VirtualLink.calcReverseStrand,T_MP_L, ...
    InputVariables=["gc","ori","rev"],OutputVariableNames=["gc","ori"]);
T_MP_L.gc=tmp.gc;
T_MP_L.ori=tmp.ori;
T_MP_L.weight=cell(16,1); %format:[v0,v1,weight,isInternal]
for ii=1:16
    tmp=[];
    for jj=1:3
        cgc=T_MP_L.gc{ii,jj};
        tmp=[tmp;0,cgc(1),0,2];
        if length(cgc)==2
            tmp=[tmp;cgc(1:2),0,0];
        end
        tmp=[tmp;cgc(end),0,0,1];
    end
    T_MP_L.weight{ii}=tmp;
end
T_MP_L.weight(ismember(T_MP_L.param,["A4";"B4"]))={nan};
T_MP_L=modify(T_MP_L);
%%
T_MP_R=table(["A";"B";"C";"D"], ...
   [{[1,-3,2],[3],[-1,-2]};
   {[-1,3,-2],[-3],[1,2]};
   {[2,-3],[1,3],[-1,-2]};
   {[-3,1],[3,2],[-1,-2]}], ...
   [[1,-1,1];[-1,1,-1];[1,1,-1];[-1,-1,1]],VariableNames=["label","gc","ori"]);
T_MP_R=repmat(T_MP_R,4,1);
T_MP_R.num=repelem(1:4,1,4)';
T_MP_R.rev=repelem([0,0,0;0,0,1;0,1,0;0,1,1],4,1);
T_MP_R.param=T_MP_R.label+T_MP_R.num;
tmp=rowfun(@VirtualLink.calcReverseStrand,T_MP_R, ...
InputVariables=["gc","ori","rev"],OutputVariableNames=["gc","ori"]);
T_MP_R.gc=tmp.gc;
T_MP_R.ori=tmp.ori;
T_MP_R.weight=T_MP_L.weight;
MPWeightData=dictionary(["A2";"B2";"C3";"D3"],  ...
                    {[-2,1];[2,1];[1,1];[2,-1]});
for ii=1:16
    arr=[];
    for jj=1:3
        cgc=T_MP_R.gc{ii,jj};
        arr=[arr;0,cgc(1),0,2];
        for kk=1:length(cgc)-1
            arr=[arr;cgc(kk:kk+1),0,0];
        end
        arr=[arr;cgc(end),0,0,1];
    end
    param=T_MP_R.param{ii};
    if ismember(param,MPWeightData.keys)
        arr(arr(:,1)==MPWeightData{param}(1),3)=-MPWeightData{param}(2);
        arr(arr(:,2)==MPWeightData{param}(1),3)=MPWeightData{param}(2);
    end
    T_MP_R.weight{ii}=arr;
end
T_MP_R.weight(ismember(T_MP_R.param,["A4";"B4"]))={nan};
T_MP_R=modify(T_MP_R);
save("moveData/MPmoveData","T_MP_R","T_MP_L")
%% CP move
LHS=VirtualLink;
LHS.setData(ostring="v_{1}^{-r}v_{2}^{+r}v_{2}^{-r}v_{3}^{+r}v_{3}^{-r}v_{1}^{+r}");
LHS.convert("mirrorMfd");
LHS.isCut=true;
RHS=VirtualLink;
RHS.setData(ostring="v_{1}^{-l}v_{2}^{-l}v_{5}^{+r}v_{3}^{-l}v_{2}^{+l}v_{3}^{+l}v_{4}^{+r}v_{4}^{-r}v_{1}^{+l}v_{5}^{-r}");
RHS.convert("mirrorMfd");
RHS.isCut=true;
save("moveData/CPmoveObj.mat","LHS","RHS")
%% Pure Sliding move
T_PS_R=combinations({[-1,-2],[-2,-1]},{[1,-1],[-1,1]});
T_PS_R=table([repmat({[1,2]},4,1),T_PS_R.Var1],[[1;1;1;1],[-1;-1;-2;-2]], ...
    string(1:4)',vertcat(T_PS_R.Var2{:}),repmat(2,4,1), ...
    VariableNames=["gc","gcFirst","param","ori","NV"]);
T_PS_L=table({[],[]},nan(1,2),zeros(1,0),0,VariableNames=["gc","gcFirst","ori","NV"]);
T_PS_R=modify(T_PS_R);
save("moveData/PSmoveData","T_PS_R","T_PS_L")
%% Bumping MP move
T_BMP_L=table("",{[-1,-2],1,2},[1,-1] ...
    ,VariableNames=["param","gc","ori"]);
T_BMP_L=modify(T_BMP_L);
T_BMP_R=table("",{-1,-2,-3,[1,3,2]},[1,1,-1] ...
    ,VariableNames=["param","gc","ori"]);
T_BMP_R=modify(T_BMP_R);
save("moveData/BMPmoveData","T_BMP_L","T_BMP_R")

%% Bumping 02 move
T_B02_L=table("",{[-1,-2],1,2},[1,-1] ...
    ,VariableNames=["param","gc","ori"]);
T_B02_L=modify(T_B02_L);
save("moveData/B02moveData","T_B02_L");
%%

function T=modify(T)
    T.gcFirst=cellfun(@(x) x(1),T.gc,"ErrorHandler",@(~,~)nan);
    T.NV=max(cellfun(@max,T.gc),[],2);
    if ~ismember("weight", T.Properties.VariableNames)
        T.weight(:) = {nan}; % Initialize weight column if not present
    else
        for i=1:height(T)
            wdata=T.weight{i};
            if any(isnan(wdata)), continue; end
            jmiddle=~wdata(:,4);
            wdata=[wdata(~jmiddle,:);wdata(jmiddle,:)];
            T.weight{i}=wdata;
        end
    end
end