%[text] # テスト
%[text] ## 初期値
a=1
%[text] ##  テスト
%%\ 
a=a+1;
assert(a==2); %[output:5a768cf0]
%[text] test2は連続的に実行すると失敗するが，独立に実施すると通る．
%% test2
a=a+2;
assert(a==3);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:5a768cf0]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"次を使用中のエラー: <a href=\"matlab:matlab.lang.internal.introspective.errorDocCallback('assert')\" style=\"font-weight:bold\">assert<\/a>\nアサーションに失敗しました。"}}
%---
