vl=VirtualLink
close all
vl.setData(Gauss={[1,-1,2,-3],[-2,3]},ori=[1,1,0])
vl.calcWeight(true)
vl.movable("MP",v=1:2)
% vl.plot
vl.move("H",p=-2,v=2)
% vl.plot
%%

vl.move("MP",p="C3",v=1:2)
% vl.plot
%%
vl.plot