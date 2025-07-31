%[text] # テスト
%[text] ## 初期値
a=1;
vl=VirtualLink;
%[text] ##  テスト
%% test1
a=a+1;
assert(a==2); %[output:5a768cf0]
%[text] test2は連続的に実行すると失敗するが，独立に実施すると通る．
%% test2
a=a+2;
assert(a==3);
%%
%[text] ## Disk of spine
%% abalone disk
vl.setData(GaussCode={[1,-1]},orientation=1)
tbl=vl.getDiskTable;
weight=-(tbl.dx\tbl.cp)'
%%

%% negative abalone disk 
vl.setData(GaussCode={[1,-1]},orientation=-1)
tbl=vl.getDiskTable;
weight=-(tbl.dx\tbl.cp)'
%%

%% abalone knot compl
vl.setData(GaussCode={[1,-1]},orientation=1)
vl.convertKnotCompl;
tbl=vl.getDiskTable;
weight=-(tbl.dx\tbl.cp)' %[output:42e04107] %[output:9c2754d9]
tbl.dx*weight'+tbl.cp %[output:3c754ea0]
dx=tbl.dx;
dx=dx(:,[1,3,5,7]) %[output:46beeff3]
dx*[0;1;1;-1]+tbl.cp %[output:7a636727]
%%

vl.setData(GaussCode={[1,-1]},orientation=1)
vl.convertKnotCompl;
[gc,ori]=vl.getGaussCode;
gc{1}(end+1)=-5;
gc{3}=[5];
ori(end+1)=1;
vl.setData(GaussCode=gc,orientation=ori)
% tbl=vl.getDiskTable;
vl.calcWeight %[output:5045b4a0] %[output:2626302a]



%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:5a768cf0]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('assert')\" style=\"font-weight:bold\">assert<\/a>\nアサーションに失敗しました。"}}
%---
%[output:42e04107]
%   data: {"dataType":"warning","outputData":{"text":"警告: ランクが欠落しています、rank = 4, tol =  3.076740e-15。"}}
%---
%[output:9c2754d9]
%   data: {"dataType":"matrix","outputData":{"columns":8,"name":"weight","rows":1,"type":"double","value":[["0.4286","0","0.4286","-1.0000","0.4286","0","0","0"]]}}
%---
%[output:3c754ea0]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":6,"type":"double","value":[["0.2857"],["-0.0000"],["0.2857"],["0.5714"],["0.5714"],["0.5714"]]}}
%---
%[output:46beeff3]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"dx","rows":6,"type":"double","value":[["1","1","1","1"],["0","0","0","0"],["1","1","1","1"],["-1","0","0","-1"],["0","-1","0","0"],["0","0","-1","0"]]}}
%---
%[output:7a636727]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"ans","rows":6,"type":"double","value":[["0"],["1"],["0"],["2"],["0"],["0"]]}}
%---
%[output:5045b4a0]
%   data: {"dataType":"warning","outputData":{"text":"警告: unknown type conversion:py.sage.modules.vector_rational_dense.Vector_rational_dense"}}
%---
%[output:2626302a]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('py.list\/double')\" style=\"font-weight:bold\">py.list\/double<\/a>\n位置 1 の Python リスト\/タプル要素の型は 'double' でなければなりません。\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('SageWrapper.toDouble', '\/Users\/nisimoriyuuya\/Desktop\/topology\/SW\/SageWrapper.m', 66)\" style=\"font-weight:bold\">SageWrapper.toDouble<\/a> (<a href=\"matlab: opentoline('\/Users\/nisimoriyuuya\/Desktop\/topology\/SW\/SageWrapper.m',66,0)\">行 66<\/a>)\n                        ret=double(arg);\n                        ^^^^^^^^^^^^^^^^\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('SageWrapper.toDouble', '\/Users\/nisimoriyuuya\/Desktop\/topology\/SW\/SageWrapper.m', 70)\" style=\"font-weight:bold\">SageWrapper.toDouble<\/a> (<a href=\"matlab: opentoline('\/Users\/nisimoriyuuya\/Desktop\/topology\/SW\/SageWrapper.m',70,0)\">行 70<\/a>)\n                    ret=SageWrapper.toDouble(list(arg));\n                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\nエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('VirtualLink\/calcWeight', '\/Users\/nisimoriyuuya\/Desktop\/topology\/Manifold\/VirtualLink\/VirtualLink.m', 1002)\" style=\"font-weight:bold\">VirtualLink\/calcWeight<\/a> (<a href=\"matlab: opentoline('\/Users\/nisimoriyuuya\/Desktop\/topology\/Manifold\/VirtualLink\/VirtualLink.m',1002,0)\">行 1002<\/a>)\n            w=SW.toDouble(C{1})\n            ^^^^^^^^^^^^^^^^^^^"}}
%---
