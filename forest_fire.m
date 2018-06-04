function forest_fire(A, pf, pb, prev_length, seeds, max_t, k, extime)
    if k > max_t
        return;
    end
    tic;
    fprintf('===========%dth round===========\n', k);
    ct = fix(clock);
    t = cputime;
    fprintf('current time: %d:%d:%d\n', ct(4), ct(5), ct(6));
    fprintf('# of seeds: %d\n', length(seeds));
    if k > 1
        fprintf('Expected time: %.2f seconds.\n', extime * prev_length);
    end
    len = size(A,1);
    next_seeds=[];
    fid = fopen('output.txt', 'a+');
    % for i = 1:length(seeds)
    %     for_nexts = find(randAf(seeds(i),:));
    %     back_nexts = find(randAb(:,seeds(i)));
    %     next_seeds = horzcat(next_seeds, for_nexts);
    %     next_seeds = horzcat(next_seeds, back_nexts');
    %     to_add = ones(1,length(for_nexts)) * seeds(i);
    %     fprintf(fid, '%8d %8d\n', [to_add; for_nexts]);
    %     to_add = ones(1,length(back_nexts)) * seeds(i);
    %     fprintf(fid, '%8d %8d\n', [back_nexts'; to_add]);
    %     % fprintf('for_nexts: %d, back_nexts: %d\n', length(for_nexts), length(back_nexts));
    %     % A(seeds(i),:) = 0;
    %     % A(:,seeds(i)) = 0;
    % end
    rand_row = (rand(1, len) <= pf);
    rand_col = (rand(len, 1) <= pb);
    RF = A;
    RB = A;
    lins = linspace(1,len, len);
    RF(~ismember(lins, seeds), :) = 0;
    RB(:, ~ismember(lins, seeds)) = 0;

    RF(seeds,:) = RF(seeds, :) .* rand_row;
    RB(:, seeds) = RB(:, seeds) .* rand_col;

    [cur_seed1 for_nexts] = find(RF);
    [back_nexts cur_seed2] = find(RB);
    A(seeds,:) = 0;
    A(:,seeds) = 0;

    fprintf(fid, '%8d %8d\n', [cur_seed1'; for_nexts']);
    fprintf(fid, '%8d %8d\n', [back_nexts'; cur_seed2']);

    next_seeds = unique([for_nexts; back_nexts]);

    % for_nexts = find(A(seeds, :) .* (rand(1,len) <= pf));
    % for_nexts = reshape(for_nexts, 1, size(for_nexts, 1)*size(for_nexts));
    % back_nexts = find(A(seeds,:) .* (rand(len, 1) <= pb));
    % back_nexts = reshape(back_nexts, 1, size(back_nexts, 1)*size(back_nexts));
    % next_seeds = horzcat(next_seeds, for_nexts);
    % next_seeds = horzcat(next_seeds, back_nexts');
    % randAf(seeds,:) = 0;
    % randAb(:,seeds) = 0;
    fclose(fid);

    Sample = uniq_sample(k);
    [s, d] = find(Sample);

    [Dt0 Dt1 Ht] = checkDH(s, d);
    fid = fopen('DHscores.txt', 'a+');
    fprintf(fid, '%5.4f %5.4f %5.4f\n', Dt0, Dt1, Ht);
    fclose(fid);

    fprintf('\n');
    toc;
    elapsed = (cputime - t)/length(seeds);
    forest_fire(A, pf, pb, length(seeds), next_seeds, max_t, k+1, elapsed);

end
