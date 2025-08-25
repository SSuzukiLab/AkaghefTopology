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
T_MP_L.detail=T_MP_L.label+T_MP_L.num;
tmp=rowfun(@VirtualLink.calcReverseStrand,T_MP_L, ...
    InputVariables=["gc","ori","rev"],OutputVariableNames=["gc","ori"]);
T_MP_L.gc=tmp.gc;
T_MP_L.ori=tmp.ori;
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
T_MP_R.detail=T_MP_R.label+T_MP_R.num;
tmp=rowfun(@VirtualLink.calcReverseStrand,T_MP_R, ...
InputVariables=["gc","ori","rev"],OutputVariableNames=["gc","ori"]);
T_MP_R.gc=tmp.gc;
T_MP_R.ori=tmp.ori;
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
%%
T_PS_R=combinations({[-1,-2],[-2,-1]},{[1,-1],[-1,1]});
T_PS_R=table([repmat({[1,2]},4,1),T_PS_R.Var1],[[1;1;1;1],[-1;-1;-2;-2]], ...
    string(1:4)',vertcat(T_PS_R.Var2{:}),repmat(2,4,1), ...
    VariableNames=["gc","gcFirst","detail","ori","NV"]);
T_PS_L=table({[],[]},nan(1,2),zeros(1,0),0,VariableNames=["gc","gcFirst","ori","NV"]);
save("moveData/PSmoveData","T_PS_R","T_PS_L")
modify(T_PS_L)
%% Bumping MP move
T_BMP_L=table("",{[-1,-2],1,2},[1,-1] ...
    ,VariableNames=["detail","gc","ori"]);
T_BMP_L=modify(T_BMP_L);
T_BMP_R=table("",{-1,-2,-3,[1,3,2]},[1,1,-1] ...
    ,VariableNames=["detail","gc","ori"]);
T_BMP_R=modify(T_BMP_R);
save("moveData/BMPmoveData","T_BMP_L","T_BMP_R")

%% Bumping 02 move
T_B02_L=table("",{[-1,-2],1,2},[1,-1] ...
    ,VariableNames=["detail","gc","ori"]);
T_B02_L=modify(T_B02_L);
save("moveData/B02moveData","T_B02_L");
%%

function T=modify(T)
    T.gcFirst=cellfun(@(x) x(1),T.gc,"ErrorHandler",@(~,~)nan);
    T.NV=max(cellfun(@max,T.gc),[],2);
end