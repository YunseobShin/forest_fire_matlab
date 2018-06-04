function Sample = uniq_sample(max_t, len)
    % fprintf('saving %dth round sample...\n', max_t);
    [s,d]=textread('output.txt','%d %d');
    nodes = vertcat(s, d);
    Nnode = length(unique(nodes));
    Sample = sparse(s, d, 1, len, len);
    Sample = double(Sample>0);

    % temp = Sample + Sample';
    % sums1 = sum(temp, 2);
    % disnode = find(sums1==0);
    % Sample(disnode, :) = [];
    % Sample(:, disnode) = [];

    % output = sprintf('./samples/FF_sample');
    % save(output, 'Sample');

    % B = largestcomponent(Sample);
    % flcc = (size(B,1)*size(B, 2))/(size(Sample, 1)*size(Sample,2));
    % [C1 C2 C] = clustCoeff(Sample);
    sv = svds(Sample, 1);

    fid = fopen('./samples/T_scores.txt', 'a+');
    fprintf(fid, '=================\n');
    fprintf(fid, 'round: %d\n', max_t);
    % fprintf(fid, 'T3: %.4f\n', flcc);
    fprintf(fid, 'T4: %.4f\n', sv);
    % fprintf(fid, 'T5: %.4f\n', C1);
    fprintf(fid, '=================\n');


end
