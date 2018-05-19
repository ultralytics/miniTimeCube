clc; clear all; close all
tic
ns = 60*60*1; %number of seconds
y=zeros(ns,1);
rate = 500; % Hz
ov=ones(rate,1);

for i=1:ns
   t=rand(rate,1);
   t=t(:,ov);
   dt = abs(t-t');
   dt(dt==0)=1;
   mdt=min(dt);
   y(i) = sum(mdt<12E-6);
end
sy = sum(y)
toc