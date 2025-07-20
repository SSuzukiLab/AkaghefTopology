classdef SmallObject
    properties
        value
    end

    methods
        function obj = SmallObject(v)
            if nargin > 0
                obj.value = v;
            end
        end

        function g = get.value(obj)
            disp("wow")
            g=obj.value;
        end
        function obj = set.value(obj,arg)
            % disp(obj.value)
            disp("a")
           obj.value = arg;
        end
    end
end