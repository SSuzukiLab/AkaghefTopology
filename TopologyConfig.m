classdef TopologyConfig < dynamicprops & matlab.mixin.SetGet
    %TOPOLOGYCONFIG Configuration manager class
    %   Inherits dynamicprops to manage XML settings as properties

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

    methods(Access=private)
        function obj = TopologyConfig()
            % This class must be singleton. get instance by accessing
            % the property .H
        end
    end
end
