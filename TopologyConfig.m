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
    end

    methods
        % function obj = TopologyConfig()
        %     % Constructor: load existing XML or create default settings
        % 
        % 
        %     obj.update()
        % end
        % function ret=struct(obj)
        %     % Convert properties to a struct
        %     ret = struct();
        %     props = properties(obj);
        %     for i = 1:length(props)
        %         ret.(props{i}) = obj.(props{i});
        %     end
        % end
        % function update(obj)
        % 
        %     config= yamlread("config.yaml");
        %     props=string(fieldnames(config))';
        %     for prop=props
        %         if ~obj.isprop(prop)
        %             addprop(prop);
        %         end
        %         obj.(prop)=config.(prop);
        %     end
        % end
        % 
        % function saveConfig(obj)
        %     % Write current Config struct to XML file
        %     yamlwrite("config.yaml", obj.struct());
        % end
    end
end
