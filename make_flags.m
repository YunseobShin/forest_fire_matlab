A = csvread('sample2.csv');

sfs = A(:,[1 5]);
dfd = A(:,[2 6]);

ns = vertcat(sfs, dfd);
mx = max(ns(:,1));
ns(:,2) = ns(:,2) * mx;
isSpam = ns(:,1) - ns(:,2);

len = max(unique(A(:,[1 2])));
SA = sparse(A(:,1), A(:,2), 1, len, len);

spams = find(unique(isSpam) > 0);

flag2 = ones(1, size(SA, 1));

flag2(spams) = 0;

save('SA.mat', 'SA');
save('flag2.mat', 'flag2');


% sfs = horzcat(s, fs);
% dfd = horzcat(d, fd);
