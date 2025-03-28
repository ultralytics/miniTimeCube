% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function a = fcnchannels2windows(a)
%pa = 1536x32768x5 window format (for 5 pedestals)
%pb = 786432x64x5 channel format

nr=786432; %ns=98304; %number of rows

n=size(a,3);
if n>1
    a=permute(reshape(permute(a,[2 1 3]),[64 nr n]),[2 1 3]);
else
    a=reshape(a',[64 nr])';
end
