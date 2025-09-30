vl=VirtualLink;
close all
vl.setData(odata={{[1,-1]},[1]});
vl.calcWeight(true);
vl.plot
%%
vl.movable("02",p="dir",v=[1,-1])
vl.move("02",p="dir",v=[1,-1])
vl.plot