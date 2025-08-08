classdef TestClassVirtualLink < matlab.unittest.TestCase
    properties (TestParameter)
    vgc={ ...
        {[]}, ... % unknot
        {[1,-1]}, ... % abalone
        {[1,-2,3,-1,2,-3]}, ... % trefoil
        {[1,-2],[-1,2]}, ... % Hopf link
        {[1],[-1]}, ... % virtual Hopf link
        {[1,-2,3,-1,2,-3]}, ... % mirror trefoil
        {[-1,2,1,-2]}, ... % virtual trefoil
        {[-1,-2,-8,3,9,1,-4,6,5,4,-6,-5,-7,8,-3,7,2,-9]}, ... % Koda's example
        {[-1,-2,3,1,-4,5,4,-5,-3,2]}  % Koda's example (virtual)
    };
    vsgn={ ...
        [], ... % unknot
        [1], ... % abalone
        [1,1,1], ... % trefoil
        [1,1], ... % Hopf link
        [1], ... % virtual Hopf link
        [-1,-1,-1], ... % mirror trefoil
        [1,1], ... % virtual trefoil
        [-1,1,1,1,-1,1,1,1,1], ... % Koda's example
        [-1,1,1,1,-1]  % Koda's example (virtual)
    };
    
    end

    methods (Test, ParameterCombination = 'sequential')
        function setGaussCode(testCase, vgc, vsgn)
            obj=VirtualLink();
            obj.setGaussCode(vgc,vsgn);
            testCase.verifyEqual(obj.GaussCode,vgc);
            testCase.verifyEqual(obj.orientation,vsgn);
        end
    end
    methods (Test, ParameterCombination ='pairwise')
        % function testSum(testCase, a, b)
        %     actual = mySumFunction(a, b);  % Replace with your actual function
        %     expected = a + b;
        %     testCase.verifyEqual(actual, expected);
        % end
    end
end
