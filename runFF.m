tic;
% Sample 3
% pre-processing the data
% Node : 955,789
% Edge : 4,150,368

[s,d,hs,hd,fs,fd]=textread('ffdata.txt','%d %d %d %d %d %d');
s=s+1;
d=d+1;
hs=hs+1;
hd=hd+1;
node = vertcat(s,d);
Nnode =length(unique(node));
A = sparse(s,d,ones(length(s),1),Nnode,Nnode);
A = double(A>0); % unweighted
fprintf('[1] *** No. of nodes: %d\n',size(A,1));
fprintf('[2] *** No. of edges: %d\n\n',nnz(A));

% 95만개 중에 2만9천개 정도 밖에 in-link와 out-link 둘 다 가지고 있음
[Dgoal0 Dgoal1 Hgoal] = checkDH(s, d, 0);
fprintf('Dgoal0: %.4f, Dgoal1: %.4f, Hgoal: %.4f\n', Dgoal0, Dgoal1, Hgoal);
% Forest Fire
fprintf('\n===Forest Fire Sampling starts.===\n');
fid = fopen('output.txt', 'w');
fclose(fid);
fid = fopen('DHscores.txt', 'w');
fprintf(fid, 'D0      D1      H\n');
fprintf(fid, '%5.4f %5.4f %5.4f\n', Dgoal0, Dgoal1, Hgoal);
fclose(fid);
fid = fopen('./samples/num_of_spams.txt', 'w');

fprintf(fid, '%4s %6s %8s %5s %6s %6s %8s %7s %5s %5s\n', 't', '#nodes', '#edges', 'density', 'normal', 'spam', 'n-n', 'n-s', 's-s', 'rate');
get_nms(s, d, 0);
fclose(fid);

% [B C] = largestcomponent(A);
% flcc = (size(B,1)*size(B, 2))/(size(A, 1)*size(A,2));
% [C1 C2 C] = clustCoeff(A);
sv = svds(A, 1);

fid = fopen('./samples/T_scores.txt', 'w');
fprintf(fid, '=================\n');
fprintf(fid, 'round: 0\n');
% fprintf(fid, 'T3: %.4f\n', flcc);
fprintf(fid, 'T4: %.4f\n', sv);
% fprintf(fid, 'T5: %.4f\n', C1);
fprintf(fid, '=================\n');
fclose(fid);

pf = 0.4;
pb = 0.35;
max_t = 100;
wff = 0.15;
wfb = ((pf-1) * wff)/(pf-1+2*wff);
wbf = wff*(pb/pf);
wbb = ((pb-1) * wbf)/(pb-1+2*wbf);
% avgp = pf/(1-pf);
% avgwf = (pf+wff)/(1-(pf+wff));
% avgwb = (pf-wfb)/(1-(pf-wfb));
% fprintf('%.3f %.3f\n', avgwf - avgp, avgp - avgwb);
rate = 4.19;
len = size(A, 1);
new_gene = 5000;

k = 1;
seed_q = [];
nol = 0;
while k <= max_t
    seeds = unique(randi(len, new_gene, 1));
    adv_FF_rate(A, A, rate, pf, pb, 0, seeds, max_t, k, 0, wff, wfb, wbf, wbb, len, new_gene, seed_q, nol);
    load round.mat;
    load seed_q.mat;
    load nol.mat;
end

% k = 1;
% seed_q = [];
% while k <= max_t
%     seeds = unique(randi(len, new_gene, 1));
%     adv_FF(A, A, Dgoal0, Dgoal1, Hgoal, pf, pb, 0, seeds, max_t, k, 0, 0, 0, wff, wfb, wbf, wbb, len, new_gene, seed_q);
%     load round.mat;
%     load seed_q.mat;
% end

toc;
