% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [] = fcnsavefile(input, flags, handles)

%CREATE NAME FOR MAT FILE
FileName = [input.cube.prettyvolume ' ' input.cube.prettyname sprintf(' %.0fpc %.0fpix %.0fqe %.0fre %s',input.cube.coverageFraction*100,input.cube.pixels, ...
    input.cube.QEmean*100, input.Material(1).mu(4)*100, input.Material(1).name) ''];

[filename,pathname] = uiputfile('*.mat','Save this Detector as:',[cd '/SAVED/' FileName '.mat']); %#ok<MCCD>
if filename(1)==0
    return
end

save([pathname filename],'input','flags');



