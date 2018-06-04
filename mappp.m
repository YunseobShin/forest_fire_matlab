%% mapping_page
[s,d,hs,hd,fs,fd]=textread('ffdata.txt','%d %d %d %d %d %d');
global T;
T = [s d];
global host; global flag; global M; global F;
M = [s d hs hd];
assert(size(T,1)==size(M,1));
F = [fs fd];
assert(size(F,1)==size(M,1));
m = size(A,1);
host = NaN(m,1);
flag = NaN(m,1);

%multiple workers
parpool(8, 'IdleTimeout', Inf);
disp('mapping');
mapping_page(A, T, M, F);
delete(gcp('nocreate'));
