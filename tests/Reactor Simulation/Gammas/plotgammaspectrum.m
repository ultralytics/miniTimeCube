clc; clear; close all
load gamma.mat
MWth=20;
plot(x,ybackground,'r',x,yreactor,'b',x,ybackground+yreactor*MWth,'g','linewidth',2); axis tight; grid on

sum(sort*MWth)
sum(ybackground)