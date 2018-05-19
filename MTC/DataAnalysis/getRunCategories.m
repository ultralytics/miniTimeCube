function [run,category] = getRunCategories()
X = GetGoogleSpreadsheet('1LEZYzBG8Ae_jQcYmxaTPxrfccSzXTMRNtyozElfS98k');
X=X(3:end,[1 13]);  
X=cellfun(@str2double,X);

i=X(:,1);
n=max(i);

run=(1:n)';
category=zeros(n,1);
category(i)=X(:,2);

