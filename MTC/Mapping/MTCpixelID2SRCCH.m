% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [S, R, C, CH, PMT, PMTR, PMTC, As, Ad, Si, CHs] = MTCpixelID2SRCCH(i)
load('MTCmapping.mat');  %X = [SRCCH PMT PMTR PMTC ASICs ASICd Si CHs FACE];
%SRCCH(SRCCH(:,1)==0,   1)=102;
%SRCCH(SRCCH(:,1)==255, 1)=13;

%SRCCH(SRCCH(:,1)==255, 1)=65; %5/2015
%SRCCH(SRCCH(:,1)==6,   1)=64;

% SRCCH(SRCCH(:,1)==255 | SRCCH(:,1)==13 | SRCCH(:,1)==65, 1)=3; %11/2015
% SRCCH(SRCCH(:,1)==101,                  1)=5;
% SRCCH(SRCCH(:,1)==102 | SRCCH(:,1)== 0, 1)=9;
% SRCCH(SRCCH(:,1)==64                  , 1)=6;

X = sortrows(X,1); %#ok<NODEF>

if nargout==1
    S = X(i,2:13);
else
    S       = X(i,2);   %SCROD              [1 2 4 6 7 8 10 11 12 255 101 0]
    R       = X(i,3);   %ASIC row           (0-3)
    C       = X(i,4);   %ASIC col           (0-3)
    CH      = X(i,5);   %ASIC channel       (0-7)
    
    if nargout>4
        PMT     = X(i,6);   %PMT                (1-24) 
        PMTR    = X(i,7);   %PMT row            (1-8)
        PMTC    = X(i,8);   %PMT col            (1-8)
        As      = X(i,9);   %ASIC on SCROD      (1-16)
        Ad      = X(i,10);  %ASIC on detector   (1-192)
        Si      = X(i,11);  %SCROD indices      (1-12)
        CHs     = X(i,12);  %CHannel on SCROD   (1-128)
        Face    = X(i,13);  %Face               (1-6)
    end
end
