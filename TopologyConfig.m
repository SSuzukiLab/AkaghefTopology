classdef TopologyConfig < dynamicprops & matlab.mixin.SetGet
    %TOPOLOGYCONFIG Configuration manager class

    properties(Constant)
        % Singleton instance
        H TopologyConfig = TopologyConfig();
        % Path of this class file
        ProjectPath = fileparts(mfilename('fullpath'));
    end

    properties
        % Display mode for URD: 0=default, 1=short, 2=long
        URDdisp = 2;
        % Display mode for updates
        URupdateDisp = false;
        VL=struct( ...
            disp=1 ...
            )
        dic dictionary
        VL_vg2g_disp=false
        VL_move_method {mustBeMember(VL_move_method,["old","new"])}="new"
    end

    methods
        function obj = TopologyConfig()
            persistent flag
            if ~isempty(flag)
                obj=obj.H;
            end
            flag=1;
        end
    end
end
