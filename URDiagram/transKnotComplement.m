function ret = transKnotComplement(arg)
    % transKnotComplement 
    % 入力arg:数値列　gauss code
    % 出力ret:数値列　o-graph gauss code
    arg=reshape(arg,1,[]);
    pos=4*arg-[2;3]*(arg>0)+[3;1]*(arg<0);
    neg=4*arg-[0;1]*(arg>0)+[2;0]*(arg<0);
    ret=[pos(1:end) 0 neg(end:-1:1)];
    % N=numel(arg)/2;
    % ret=
end

