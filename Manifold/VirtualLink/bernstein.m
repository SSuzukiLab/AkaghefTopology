function BE = bernstein(NN, TT)
% BERNSTEIN バーンスタイン基底関数．
% 
%	BERNSTEIN(N, T) 
%	パラメータTに対するN次バーンスタイン基底関数の各項が配列として返ります．
%
%	例 :
%	t=0.5の3次バーンスタイン関数の各項．
%		bernstein(3, 0.5)
%
%	t=[0, 1/4, 1/2, 3/4, 1]の4次バーンスタイン関数の各項．
%		bernstein(4, [0:1/4:1])
%
%	参考 BEZIER PASCAL ROT90 DIAG 

% --
%	Title : BERNSTEIN()
%	Author : Sach1o : hTTp://sach1o.blog80.fc2.com/
%	Created : 2006/12/06
% //--

error(nargchk(2, 2, nargin));
if NN<1 | int32(NN)~=NN | ~isscalar(NN) error('第1引数は2以上の整数である必要があります．'); end;
if ~isreal(TT) | any(TT<0) | any(1<TT) error('第2引数は0以上1以下の実数配列である必要があります。'); end;

onest = ones(1, length(TT));
II = [0:NN]'*onest;
BI = diag(rot90(pascal(NN+1)))*onest;
TT = ones(NN+1,1)*TT;

BE = (1 - TT).^(NN - II).*BI.*(TT.^II);
