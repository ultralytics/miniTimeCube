% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function a = fcnwindows2channels(a)
%pa = 1536x32768x5 window format (for 5 pedestals)
%pb = 786432x64x5 channel format

ns=32768; %ns=4096; %number of samples
n=size(a,3);
if n>1
    a=permute(reshape(permute(a,[2 1 3]),[ns 1536 n]),[2 1 3]);
else
    a=reshape(a',[ns 1536])';
end
