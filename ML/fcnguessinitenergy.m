% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

 function e0hat = fcnguessinitenergy(k,zN,p)
%p = sqrt(zN)'*input.cube.all.xyz/sum(sqrt(zN)); %position guess
k.reflections=0;

fall = fcnsolidanglevectorized(k,p);
e0hat = sum(zN(:))./(sum(fall(:))*(k.yield + 1));
