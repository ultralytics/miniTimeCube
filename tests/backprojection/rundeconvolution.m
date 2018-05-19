clc
close all

load I.mat
fig
[nr,nc]=size(I);
%I0 = zeros(150,150);
%I0(11:10+nr,11:10+nc) = I;
%I = I0;
%[nr,nc]=size(I);

imshow(I,[]); title('Raw Image'); hold on; plot(84,48,'+g')

load PSF.mat
PSF = PSF(34:96,34:96); 

%PSF = I(31:93,1:53);
%PSF = fspecial('gaussian',121,5); %must be odd!
%PSF = PSF.^.0001;
PSF = PSF./sum(sum(PSF));
set(figure,'Position',[340 50 300 300]);
imshow(PSF,[]); title('PSF'); hold on; plot(64,64,'+g')

I = edgetaper(I,PSF);
%filter = 4;
%switch filter
%    case 1
        INITPSF=PSF;
        [J P]= deconvblind(I,INITPSF,30);
        set(figure,'Position',[10 350 300 300]); imshow(J,[]); title('Deconvolved Blind'); hold on; plot(84,48,'+g')
%    case 2
        J = deconvlucy(I, PSF);
        set(figure,'Position',[340 350 300 300]); imshow(J,[]); title('Deconvolved Lucy'); hold on; plot(84,48,'+g')
%    case 3
        J = deconvreg(I,PSF);
        set(figure,'Position',[670 350 300 300]); imshow(J,[]); title('Deconvolved Regular'); hold on; plot(84,48,'+g')
%    case 4
        NSR = .01;
        J = deconvwnr(I,PSF, NSR);
        set(figure,'Position',[1000 350 300 300]); imshow(J,[]); title('Deconvolved Weiner'); hold on; plot(84,48,'+g')
%end

% %JOCHER MANUAL DECONVOLUTION
% 
% [nrp,ncp]=size(PSF); %number of rows and columns in PSF
% cr = ceil(nrp/2); %center row
% cc = ceil(ncp/2); %center col
% sr = cr - 1; %number of side rows
% sc = cc - 1; %number of side rows
% 
% 
% for r = 1:nc
%     r
%     pr1 = max(cr - r + 1, 1);
%     pr2 = min(nr - r + cr, nrp);
%     ir1 = max(r - sr, 1); %image row 1
%     ir2 = min(r + sr, nr);
%     
%     for c = 1:nc
%         clc
%         
%         pc1 = max(cc - c + 1, 1);
%         pc2 = min(nc - c + cc, ncp);
%         ic1 = max(c - sc, 1);
%         ic2 = min(c + sc, nc);
% 
%         
%         I1 = I(ir1:ir2,ic1:ic2);
%         P1 = PSF(pr1:pr2,pc1:pc2);
%         
%         T = I1.*P1;
%         sumT = sum(sum(T));
%         
%         I1 = I1 - T;
% 
%         I(ir1:ir2,ic1:ic2) = I1;
%         I(r,c) = I(r,c) + sumT;
% 
%     end
% end
% set(figure,'Position',[1540 50 500 500]);
% imshow(I,[]); title('Manual Deconvolved Image'); hold on; plot(84,48,'+g')