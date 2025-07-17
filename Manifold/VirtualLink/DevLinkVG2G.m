% unknot
gc={[]};
sgn=[];
% abalone
gc={[1],[-1]};
sgn=[1];
% trefoil
gc={[1],[-1]};
sgn=[1];
% Hopf link
gc={[1],[-1]};
sgn=[1];

if isempty(sgn)
    assert(any(cellfun(@numel,gc))==0)
    vgc=gc;
    vsgn=sgn;
    return
end
regV={};
regE={};
vgc={};
vsgn=sgn;
usedV=[];
cntE=1;
faceG=graph;
for si=1:length(gc) %[output:group:576e22c6]
    cgc=gc{si};
    if isempty(cgc)
        vgc{si}=[];
        continue
    end
    if ismember(cgc(1),usedV)
        [fidx,idx]=findC(cgc(1))
    else
        fidx=length(regV)+1;
        regV{fidx}=[-1,0,1]*sgn(cgc(1)); %[output:8f1e3b06]
        regE{fidx}=[cntE,0,0,cntE];
    end
    for vi=2:length(cgc)

    end






end %[output:group:576e22c6]
function [idxC,idx]=findC(key,arr)
    idxC=find(cellfun(@(C)ismember(key,C),arr));
    idx=find(key==arr(idx));
end


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:8f1e3b06]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"配列インデックスは正の整数または logical 値でなければなりません。"}}
%---
