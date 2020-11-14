a = [0,1,2,3;4,5,6,7;8,9,10,11];

a1 = a(:).'; %矩阵列向量 按行排出来
a2 = a(:); %矩阵列向量 按列排出来
a3 = a.'; %矩阵转置
%a4 = a(:,);
b1 = fliplr(a);    % For row-vector a.
b2 = flipud(a);    % For column-vector a.
b3 = wrev(a);      % For any vector a.
b4 = a(end:-1:1); % This is the implementation of function wrev.
