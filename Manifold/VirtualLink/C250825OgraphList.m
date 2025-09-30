N=2;
Tgc=table;
for ii=1:N
    Tgc=[Tgc;allGaussCode(ii)];
end
%%
Tgc.isClosed(:)=false;
for ii=1:height(Tgc)
    vl.setData(RG=Tgc.GaussCode(ii),ori=Tgc.Orientation{ii});
    Tgc.isClosed(ii)=vl.isClosed;
end

function tbl = allGaussCode(n)
    % Generate all Gauss codes with orientation for n crossings (virtual knots)
    % Reference: VirtualLink.m
    % @(n)factorial(2*n-1)/factorial(n-1)*2^n;(爆発)
    %   combination: 2 ,24 ,480 ,13440 ,483840 ,21288960
    h=factorial(2*n-1);
    GaussCode =[ones(h,1),perms([-n:-1,2:n])]; % Initialize the Gauss code with numbers from 1 to n
    GaussCode=mat2cell(GaussCode,ones(1,h),2*n);
    idx=cellfun(@(arr)isequal(arr(arr>0),1:n),GaussCode);
    GaussCode = GaussCode(idx); % Filter Gauss codes based on the index
    tmp=repmat({[-1,1]},1,n);
    tmp=combinations(tmp{:});
    Orientation=mat2cell(tmp{:,:},ones(1,2^n),n);
    tbl=combinations(GaussCode,Orientation);
    tbl.NV(:)=n;
    % tbl=table(GaussCode,Orientation);

end
