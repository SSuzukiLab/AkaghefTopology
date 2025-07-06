classdef TopologyConfig<handle& matlab.mixin.SetGet
    %TOPOLOGYCONFIG このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties(Constant)
        H TopologyConfig=TopologyConfig() 
        ProjectPath=fileparts(mfilename("fullpath"))
    end
    properties
        URDdisp=2 % 0: default display, 1: short text, 2: long text
        URupdateDisp=1
    end
    
    methods

    end
end

