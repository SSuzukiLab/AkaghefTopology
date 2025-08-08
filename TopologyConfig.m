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
        URupdateDisp = 1;
        VL=struct( ...
            disp=1 ...
            )
        dic dictionary
    end

    methods
        function obj = TopologyConfig()
            persistent flag
            if ~isempty(flag)
                error(['TopologyConfig is a singleton. ' ...
                    'Use TopologyConfig.H to access the instance.']);
            end
            flag=1;
        end
    end
end
