% Ultralytics 🚀 AGPL-3.0 License - https://ultralytics.com/license

function [b, i] = fcnpruninglist(a)
if nargin==0; a=(1:1536)'; end

[b, i] = fcnprune(a); %b are the survivors and i are survivor indices
end


function [survivors, i] = fcnprune(a)
b = true(1536,1);

b(pruneasics) = false;
i = b(a);
survivors = a(i,:);
end


function i = pruneasics
S=MTCpixelID2SRCCH(1:1536);
i=false(1536,1);

% %ASICS
i = i | (S(:,2)==1 & S(:,3)==3)...
     | (S(:,1)==10 & S(:,2)==0 & S(:,3)==3)...  %PMT 20 bad asic in corner
     | (S(:,1)==3 & S(:,2)==0)... %SCROD 3 Carrier 0 bad 2017.5.17
     | (S(:,1)==3 & S(:,2)==1); %SCROD 3 Carrier 0 bad 2017.5.17
 
%PIXELS
i([48 56 63 64, 1217 1218 1225, 141, 1472])=true; %detaching pixels at corner of PMT 10 and PMT 21, bad pixel PMT 20, PMT 1

%PMTS
%i = i | any(S(:,5)==[3 4 5 6 8 11 12 13 14 15 16 17 18 23 24],2);  %PMT 6, 16 added April 2017
i = i | any(S(:,5)==[3 4 5 6 8 10 11 12 13 14 15 16 17 18 22 23 24],2);  %PMT 22 added April 28 2017, 10 dead
%i = i | any(S(:,5)==[3 4 5 8 10 11 12 13 14 15 16 17 18 22 23 24],2);  %NEED PMT 23, 24 as NTC SCROD
end