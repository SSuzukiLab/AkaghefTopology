function BEZ = bezier(TT, CPS)
%BEZIER ベジエ曲線の座標．
% 
%	BEZIER(TT, CPS) 
%	CPSを制御点とするN次ベジエ曲線のパラメータTにおける座標を返す．
%
%	パラメータTTは，0以上1以下の実数配列で与える必要がある．
%	制御点列CPSは，1列で各制御点を表し，3行以上の実数配列で与える必要がある．
%
%	例 :
%	平面の3次ベジエ曲線のt=1/2の座標．
%		bezier(1/2,[[0;0],[1;1],[2;0]])
%
%	3次元空間の5次ベジエ曲線のt=[0,1/2,1]の座標．
%		bezier([0:1/2:1],[[0;0;0],[1;1;1],[2;0;3],[3;2;1],[4;1;0],[5;1;1]])
%
%	参考 BINOMIAL

% --
%	Title : BEZIER() 
%	Author : Sach1o : hTTp://sach1o.blog80.fc2.com/
%	Created : 2006/12/06
% //--

error(nargchk(2, 2, nargin));
if (~isreal(TT) | any(TT<0) | any(1<TT)) error('第1引数は0以上1以下の実数配列である必要があります。'); end;
if (~isreal(CPS) | length(CPS(1,:))<2) error('第2引数は少なくとも2行以上の実数配列である必要があります。'); end;
if size(CPS,1)==2, BEZ=CPS(1,:)+
BEZ = CPS*bernstein(length(CPS(1,:))-1, TT);
