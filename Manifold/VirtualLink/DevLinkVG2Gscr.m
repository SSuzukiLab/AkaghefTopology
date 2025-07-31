arg=2;
switch arg
    case 1 % unknot
        vgc={[]};
        vsgn=[];
    case 2 % abalone
        vgc={[1,-1]};
        vsgn=[1];
    case 3 % trefoil
        vgc={[1,-2,3,-1,2,-3]};
        vsgn=[1,1,1];
    case 4 % Hopf link
        vgc={[1,-2],[-1,2]};
        vsgn=[1,1];
    case 5 % virtual Hopf link
        vgc={[1],[-1]};
        vsgn=[1];
    case 6 % mirror trefoil
        vgc={[1,-2,3,-1,2,-3]};
        vsgn=[-1,-1,-1];
    case 9 % virtual trefoil
        vgc={[-1,2,1,-2]};
        vsgn=[1,1];
    case 7 % Koda's example
        vgc={[-1,-2,-8,3,9,1,-4,6,5,4,-6,-5,-7,8,-3,7,2,-9]};
        vsgn=[-1,1,1,1,-1,1,1,1,1];
        % [[-6 -4], [-3 -8], [-5 -6 -4], [-7 -3 -8], [-5 -6], [-2 -7 -8], [-9 -2 -8 -3], [-9 -3 -7 -5 -4 -1], [-2 -1 -4 -6 -5 -7], [-9 -1 -2], [-1 -9]]
        % [[-10 7], [-14 3], [-11 -7 -9], [-15 -3 -13], [-8 11], [-16 13 -2], [-17 2 14 4], [-4 15 -12 9 -6 -18], [-1 6 10 8 12 16], [5 1 17], [-5 18]]
        %  1     1     2     2     3     3     4     4     5     6     6     7     7     8     8     8     9    10
        %  3     9     4     7     5     8     6     8     9     7     9     8    10     9     9    11    10    11
        %  7    10     3    14    11     9    13    15     8     2    16     4    17     6    12    18     1     5
    case 8 % Koda's example (virtual)
        vgc={[-1,-2,3,1,-4,5,4,-5,-3,2]};
        vsgn=[-1,1,1,1,-1];
    case 10
        vgc={[3,2],[1,-3],[-1,-2]};
        vsgn=[-1,1,1];

        
end

[gc,sgn]=vg2g(vgc,vsgn);
vl.setData(Ga=gc,ori=sgn)
%%
vl.plot


function [gc,sgn]=vg2g(vgc,vsgn)
    
    
    if isempty(vsgn)
        assert(~any(cellfun(@numel,vgc)))
        gc=vgc;
        sgn=vsgn;
        return
    end
    % store information about the regions of the virtual graph
    % using vertex/edge/face (regV/regE/regF)
    % regF is a graph object representing the adjacency of regions
    % and its Weights are the edge IDs connecting the regions
    regV={};
    regE={};
    gc=vgc;
    ge=cell(size(gc));
    sgn=vsgn;
    usedV=0;
    cntE=0;
    cntV=max(cellfun(@max,vgc))+1;
    regF=graph;
    for si=1:length(vgc)
        % loop through each strand
        cgc=[vgc{si} 0];
        if isscalar(cgc)
            vgc{si}=[];
            continue
        end
        cntE=cntE+1;
        ge{si}(end+1)=cntE;
        v1=abs(cgc(1));
        if ismember(v1,usedV)
            % start strand sweeping from a used vertex
            [f2,vidx2]=findC(-v1,regV);
            regV{f2}=[cntV+1,regV{f2}([vidx2:end,1:vidx2])];
            regE{f2}=[-cntE,regE{f2}(1:vidx2-1),regE{f2}(vidx2:end),cntE];
            regF=regF.addedge(f2,[f2 f2],[0 cntE]);
            [f1,vidx]=findC(v1,regV);
            regV{f1}=[regV{f1}(1:vidx-1),-v1,0,-v1,regV{f1}(vidx+1:end)];
            regE{f1}=[regE{f1}(1:vidx-1),0,0,regE{f1}(vidx:end)];
            % regF=regF.addedge(f1,f1,0);
            f1=f2;
        else
            % start strand sweeping from a new vertex
            f1=length(regV)+1;
            regV{f1}=[0,-1,0,1]*vsgn(v1)*cgc(1)+[cntV+1,0,0,0];%?
            regE{f1}=[-cntE,0,0,cntE];
            regF=regF.addedge(f1,[f1 f1],[0 cntE]);
            usedV(end+1)=v1;
        end
        for vi=2:length(cgc)
            % loop through each vertex in the strand
            fd
            cntE=cntE+1;
            ge{si}(end+1)=cntE;
            v1=abs(cgc(vi));
            if ismember(v1,usedV)
                % pattern 1,2,3,4,6
                [f2,vidx2]=findC(v1,regV);
                fnew=findC(-v1,regV);
                % Find the shortest path from f1 to f2 (and fnew) in the graph
                path = [shortestpath(regF, f1, f2) fnew];
                len=length(path);
                if len==1
                    % pattern 1: connect to another component
                    % f1⋯>f2->fnew
                    [fnew,vidx3]=findC(-v1,regV);
                    regV{f2}=[regV{f2}(1:vidx2-1),-v1,regV{f1}(2:end),-v1,regV{f2}(vidx2+1:end)];
                    regV{fnew}=[cntV+1,-v1,regV{fnew}([vidx3+1:end,1:vidx3-1]),-v1];
                    regV(f1)=[];
                    regE{f2}=[regE{f2}(1:vidx2-1),regE{f1},regE{f2}(vidx2:end)];
                    regE{fnew}=[-cntE,regE{fnew}([vidx3:end,1:vidx3-1]),cntE];
                    regE(f1)=[];
                    % contract the faces f1<->f2
                    ns=regF.neighbors(f1);
                    for n=ns(ns~=f1)'
                        bnd=regF.Edges.Weight(regF.findedge(f1,n));
                        regF=regF.addedge(f2,n*ones(size(bnd)),bnd);
                    end
                    bnd=regF.Edges.Weight(regF.findedge(f1,f1));
                    regF=regF.addedge(f2,f2*ones(size(bnd)),bnd);
                    regF=regF.rmnode(f1);
                    regF=regF.addedge(fnew,fnew,cntE);
                    f1=fnew;
                elseif path(1)==path(2)&&vi<length(cgc)
                    % pattern 2,3: connect to vertex
                    vidx=vidx2;
                    [f2,vidx2]=findC(-v1,regV);
                    fnew=length(regV)+1;
                    if vidx2<vidx
                        % pattern 2:
                        regV{fnew}=regV{f1}([vidx2,vidx+1:end]);
                        regV{f1}=regV{f1}([1,vidx2:vidx-1,vidx2,2:vidx2]);
                        regE{fnew}=regE{f1}(vidx:end);
                        regE{f1}=[-cntE regE{f1}([vidx2:vidx-1,1:vidx2-1]) cntE];
                    elseif vidx<vidx2
                        % pattern 3:
                        regV{fnew}=regV{f1}([vidx2,2:vidx-1]);
                        regV{f1}=regV{f1}([1,vidx2:end,vidx2,vidx+1:vidx2]);
                        regE{fnew}=regE{f1}(1:vidx-1);
                        regE{f1}=[-cntE regE{f1}([vidx2:end,vidx:vidx2-1]) cntE];
                    end
                    newbnd=abs(regE{fnew});
                    regF=updateFace(regF,newbnd,f1,f2,fnew,cntE);
                    f1=f2;
                else
                    % pattern 4,6
                    ge{si}(end)=[];
                    cntE=cntE-1;
                    for fj=1:len-1
                        f1=path(fj);
                        f2=path(fj+1);
                        if fj<len-1
                            % add virtual crossing if we have to cut through edges
                            cntE=cntE+1;
                            bnds=regF.Edges.Weight(findedge(regF,f1,f2));
                            bnd=bnds(1); %ここが選択の余地あり
                            vidx=find(bnd==abs(regE{f1}),1)+1;
                            eo=regE{f1}(vidx-1)>0;
                            es=2*eo-1;
                            regV{f1}=[regV{f1}(1:vidx-1),cntV,regV{f1}(vidx:end)];
                            regE{f1}=[regE{f1}(1:vidx-1-~eo),es*cntE,regE{f1}(vidx-1+eo:end)];
                            vidx2=find(bnd==abs(regE{f2}),1)+1;
                            regV{f2}=[regV{f2}(1:vidx2-1),-cntV,regV{f2}(vidx2:end)];
                            regE{f2}=[regE{f2}(1:vidx2-1-eo),-es*cntE,regE{f2}(vidx2-1+~eo:end)];
                            regF=regF.addedge(f1,f2,cntE);
                            vs=cntV;
                            [~,gcidx]=ismember(cgc(vi),[gc{si} 0]);
                            gc{si}=[gc{si}(1:gcidx-1),-es*vs,gc{si}(gcidx:end)];
                            [sidx,geidx]=findC(bnd,ge);
                            ge{sidx}=[ge{sidx}(1:geidx), cntE, ge{sidx}(geidx+1:end)];
                            gc{sidx}=[gc{sidx}(1:geidx), es*vs, gc{sidx}(geidx+1:end)];
                            sgn(cntV)=0;
                            v1=cntV;
                            cntV=cntV+1;
                        elseif vi==length(cgc)
                            % connect ends
                            break
                        else
                            v1=abs(cgc(vi));
                            [~,vidx]=ismember(v1,regV{f1});
                            [f2,vidx2]=findC(-v1,regV);
                        end
                        cntE=cntE+1;
                        fnew=length(regV)+1;
                        % pattern 4: cut through edges and add virtual crossing
                        regV{fnew}=[-v1 regV{f1}(vidx+1:end)];
                        regV{f1}=[-v1 regV{f1}(2:vidx-1)];
                        regV{f2}=[cntV+1 regV{f2}([vidx2:end,1:vidx2])];
                        regE{fnew}=regE{f1}(vidx:end);
                        regE{f1}=regE{f1}(1:vidx-1);
                        regE{f2}=[-cntE regE{f2}([vidx2:end,1:vidx2-1]) cntE];
                        regF=updateFace(regF,abs(regE{fnew}),f1,f2,fnew,cntE);
                        f1=f2;
                        ge{si}(end+1)=cntE;
                    end
                end
            else
                % pattern 5: add new crossing to temporary graph
                usedV(end+1)=v1;
                vs=vsgn(v1)*cgc(vi);
                regV{f1}=[cntV+1 -vs regV{f1}(2:end) vs];
                regE{f1}=[-cntE regE{f1} cntE];
                regF=regF.addedge(f1,f1,cntE);
            end
            
        end
        % pattern 6: connect ends of the strand
        vidx=find(regV{f1}==0,1);
        fnew=length(regV)+1;
        regV{fnew}=regV{f1}(vidx+1:end);
        regE{fnew}=regE{f1}(vidx+1:end);
        regV{f1}=regV{f1}(2:vidx-1);
        regE{f1}=regE{f1}([2:vidx-2,1]);
        regF=updateFace(regF,abs(regE{fnew}),f1,f2,fnew,cntE);
        regF=regF.rmedge(f2,f2);
        fd
        % regF.plot
        dispC(ge)
        dispC(gc)
    end
    
    function [idxC,idx]=findC(key,arr)
        % find the index of key in arr={[...],...,[...]}]}
        idxC=find(cellfun(@(C)ismember(key,C),arr));
        assert(isscalar(idxC))
        idx=find(key==arr{idxC});
    end
    function regF=updateFace(regF,bnd,f1,f2,fnew,cntE)
        % update the face in regF with the new boundary bnd
        % bnd=unique(bnd);
        [~,notSelfAdj]=ismember(bnd,bnd);
        bnd2=bnd;
        bnd2(notSelfAdj)=[];
        bnd=setdiff(bnd,bnd2);
        [~,eidx_]=ismember(bnd,regF.Edges.Weight);
        assert(all(eidx_~=0))
        BT=regF.Edges.EndNodes(eidx_,:);
        bidx=(BT(:,1)==f1).'+1;
        fbnd=arrayfun(@(x,y)BT(x,y),1:length(eidx_),bidx);
        regF=regF.rmedge(eidx_);
        regF=regF.addedge(fnew,fbnd,bnd);
        regF=regF.addedge(f2,f2,cntE);
        if ~isempty(bnd2)
            [~,eidx_]=ismember(bnd2,regF.Edges.Weight);
            regF=regF.rmedge(eidx_);
            regF=regF.addedge(fnew,fnew*ones(size(bnd2)),bnd2);
        end
        % disp(regF.Edges)
    end
    function fd()
        % for debugging
        try
            fprintf("f1=%d\n",f1)
            fprintf("si=%d, vi=%d, v_gc=%d, s_gc=%d\n",si,vi,cgc(vi),vsgn(abs(cgc(vi))))
        catch
        end
        dispC(ge)
        dispC(regV)
        dispC(regE)
        disp(regF.Edges{:,:}')
        
    end
    function dispC(arg)
        disp("["+join(cellfun(@mat2str,arg,Un=0),", ")+"]")
    end
end