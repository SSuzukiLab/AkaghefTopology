T_MP_L=table(["A";"B";"C";"D"], ...
   [{[-1],[1,2],[-2]};
   {[1],[-1,-2],[2]};
   {[-1,2],[1],[-2]};
   {[2,-1],[1],[-2]}], ...
   [[-1,-1];[-1,-1];[-1,1];[1,-1]],VariableNames=["label","gc","ori"]);
T_MP_L=repmat(T_MP_L,4,1);
T_MP_L.num=repelem(1:4,1,4)';
T_MP_L.rev=repelem([0,0,0;0,0,1;0,1,0;0,1,1],4,1);
T_MP_L.name=T_MP_L.label+T_MP_L.num;
tmp=rowfun(@VirtualLink.calcReverseStrand,T_MP_L, ...
    InputVariables=["gc","ori","rev"],OutputVariableNames=["gc","ori"]);
T_MP_L.gc=tmp.gc;
T_MP_L.ori=tmp.ori;
T_MP_L.gcFirst=cellfun(@(x) x(1),T_MP_L.gc);

T_MP_R=table(["A";"B";"C";"D"], ...
   [{[1,-3,2],[3],[-1,-2]};
   {[-1,3,-2],[-3],[1,2]};
   {[2,-3],[1,3],[-1,-2]};
   {[-3,1],[3,2],[-1,-2]}], ...
   [[1,-1,1];[-1,1,-1];[1,1,-1];[-1,-1,1]],VariableNames=["label","gc","ori"]);
T_MP_R=repmat(T_MP_R,4,1);
T_MP_R.num=repelem(1:4,1,4)';
T_MP_R.rev=repelem([0,0,0;0,0,1;0,1,0;0,1,1],4,1);
T_MP_R.name=T_MP_R.label+T_MP_R.num;
tmp=rowfun(@VirtualLink.calcReverseStrand,T_MP_R, ...
InputVariables=["gc","ori","rev"],OutputVariableNames=["gc","ori"]);
T_MP_R.gc=tmp.gc;
T_MP_R.ori=tmp.ori;
T_MP_R.gcFirst=cellfun(@(x) x(1),T_MP_R.gc);

save("moveData/MPmoveData","T_MP_R","T_MP_L")

%%
vl.setData(ostring="v_{1}^{-r}v_{2}^{+r}v_{2}^{-r}v_{3}^{+r}v_{3}^{-r}v_{1}^{+r}");
vl.convert("mirrorMfd");
vl.disp1
%%

vl.setData(ostring="v_{1}^{-l}v_{2}^{-l}v_{5}^{+r}v_{3}^{-l}v_{2}^{+l}v_{3}^{+l}v_{4}^{+r}v_{4}^{-r}v_{1}^{+l}v_{5}^{-r}");
vl.convert("mirrorMfd");
vl.getGaussCode