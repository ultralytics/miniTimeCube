% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [pdf, strip, pid] = fcnanalogvoltage(input,flags,PE,t,ystrip,i,addnoiseflag)
pmt = input.cube.pmt;

if addnoiseflag
    amplitude = PE.amplitude(i);
else
    [~,mi]=max(pmt.amplitude.pdf); amplitude = pmt.amplitude.x(mi);
end


if flags.status.lappd
    pid=[];  strip=[];
    nt = numel(t);
    ni = numel(i);  if ni==0; pdf = zeros(ni,nt); return; end %number of pmts/lappds
    ny = numel(ystrip); %number of positions down the strip (2 means just the ends)
    np = input.cube.pixels;
    ns = pmt.strips;
    speed = 190; %mm/ns down strip
    
    dt = (t - 3*pmt.spr) - PE.t(i); %ns
    yhit = PE.xpixel(i,2)/speed;
    ystrip = ystrip/speed; %2
    dist = abs(ystrip-yhit);
    x = dist - reshape(dt,[ni 1 nt]); %[ni 2 (1)] - [ni 1 nt] = [ni 2 nt]
    
    %ALONG TEMPLATE
    pdf = interp1c(-pmt.waveform.x,pmt.waveform.y,x).*amplitude;
    %12 x   2 x 256; %12 lappds,   2 positions, 256 time
    %12 x 100 x   1; %12 lappds, 100 positions, 1 time
    
    %ACROSS TEMPLATE (normal, 0.6 strip-gaps 1sigma)
    ag = pmt.pdfacross; %(mm) sigma
    acrosspdf =  fcnpdfnorm(pmt.stripx,PE.xpixel(i,3),ag);%/fcnpdfnorm(0,0,ag); %normalizes pdf max=1
    %strip = fcnindex1c(input.cube.pmt.stripx,PE.xpixel(i,3));

    na = 5; %number of brightest strips for acrosspdf
    [I,J] = sort(acrosspdf,2,'descend');  acrosspdf = I(:,1:na); 
    acrosspdf=acrosspdf./sum(acrosspdf,2); %normalizes 26-strip integral to 1
    ustrip=J(:,1:na);
    ustrip = ustrip + (PE.pixel(i)-1)*ns;
        
    pdf = zeros(1,na) + reshape(pdf,[ni 1 nt ny]);  pdf = reshape(  pdf, [na*ni ny nt]);
    
    ustrip = ustrip(:); %unique strip
    [~, I, J] = fcnunique(ustrip);  pid = ceil(ustrip(I)/ns);
    [strip,~] = fcnind2sub([ns np*ns*2],ustrip);
    
    pdf = fcnaccumrows(pdf,J,acrosspdf);  strip=strip(I);
    return
end

dt = t - PE.t(i); %ns
%pdf = interp1c(pmt.waveform.x,pmt.waveform.y,dt);  %pdf = pdf.*amplitude;  %vary amplitude
F = griddedInterpolant(pmt.waveform.x,pmt.waveform.y,'linear','none');  pdf=reshape(F(dt(:)),size(dt));  pdf(isnan(pdf))=0;
pdf = pdf.*amplitude;





