function adv_FF(oA, A, Dg0, Dg1, Hg, pf, pb, prev_length, seeds, max_t, k, d0, d1, h, wff, wfb, wbf, wbb, len, new_gene, seed_q)
    if k > max_t
        k = max_t;
        save('round', 'k');
        return;
    end
    tic;
    seed_q = vertcat(seed_q, seeds);
    save('seed_q', 'seed_q');
    load flag.mat
    flag(find(flag)) = 2;
    flag(find(flag==0)) = 1;    % normal
    flag(find(flag==2)) = 0;    % spam
    spams = find(flag==0);
    norms = find(flag);

    % load flag2.mat
    % spams = find(flag2==0);
    % norms = find(flag2);

    fprintf('===========%dth round===========\n', k);
    ct = fix(clock);
    t = cputime;
    fprintf('current time: %d:%d:%d\n', ct(4), ct(5), ct(6));
    fprintf('# of seeds: %d\n', length(seeds));

    len = size(A,1);
    next_seeds=[];
    fid = fopen('output.txt', 'a');

    rand_row = zeros(1, len);
    rand_col = zeros(len, 1);

    spam_seeds = intersect(seeds, spams);
    norm_seeds = intersect(seeds, norms);

    RF = A;
    RB = A;
    lins = linspace(1,len, len);
    RF(~ismember(lins, seeds), :) = 0;
    RB(:, ~ismember(lins, seeds)) = 0;

    if k > 1
        if d0 > 0.1  % Dg >> D(t-1)
            xss = min( geornd(1 - (pf + weight), length(spam_seeds), 1), sum(RF(spam_seeds, spams), 2));
            yss = min( geornd(1 - (pb + weight), 1, length(spam_seeds)), sum(RB(spams, spam_seeds), 1));
            % RF(spam_seeds, spams) = RF(spam_seeds, spams) .* (rand(1, length(spams)) <= pf + weight);
            % RB(spams, spam_seeds) = RB(spams, spam_seeds) .* (rand(length(spams), 1) <= pb + weight);
        elseif d0 < -0.1
            xss = min( geornd(1 - (pf - weight), length(spam_seeds), 1), sum(RF(spam_seeds, spams), 2));
            yss = min( geornd(1 - (pb - weight), 1, length(spam_seeds)), sum(RB(spams, spam_seeds), 1));

            % RF(spam_seeds, spams) = RF(spam_seeds, spams) .* (rand(1, length(spams)) <= pf - weight);
            % RB(spams, spam_seeds) = RB(spams, spam_seeds) .* (rand(length(spams), 1) <= pb - weight);
        else
            xss = min( geornd(1 - pf, length(spam_seeds), 1), sum(RF(spam_seeds, spams), 2));
            yss = min( geornd(1 - pb, 1, length(spam_seeds)), sum(RB(spams, spam_seeds), 1));
            % RF(spam_seeds, spams) = RF(spam_seeds, spams) .* (rand(1, length(spams)) <= pf);
            % RB(spams, spam_seeds) = RB(spams, spam_seeds) .* (rand(length(spams), 1) <= pb);
        end
        for i = 1:length(spam_seeds)
            if xss & length(find(RF(spam_seeds(i), spams)))
                randsel_xss = randsample(find(RF(spam_seeds(i), spams)), xss(i));
                RF(spam_seeds(i), spams) = 0;
                RF(spam_seeds(i), randsel_xss) = 1;
            else
                RF(spam_seeds(i), spams) = 0;
            end

            if yss & length(find(RB(spams, spam_seeds(i))))
                randsel_yss = randsample(find(RB(spams, spam_seeds(i))), yss(i));
                RB(spams, spam_seeds(i)) = 0;
                RB(randsel_yss, spam_seeds(i)) = 1;
            else
                RB(spams, spam_seeds(i)) = 0;
            end
        end
        if d1 > 0.1
            xnn = min( geornd(1 - (pf + weight), length(norm_seeds), 1), sum(RF(norm_seeds, norms), 2));
            ynn = min( geornd(1 - (pf + weight), 1, length(norm_seeds)), sum(RB(norms, norm_seeds), 1));
            % RF(norm_seeds, norms) = RF(norm_seeds, norms) .* (rand(1, length(norms)) <= pf + weight);
            % RB(norms, norm_seeds) = RB(norms, norm_seeds) .* (rand(length(norms), 1) <= pb + weight);
        elseif d1 < -0.1
            xnn = min( geornd(1 - (pf - weight), length(norm_seeds), 1), sum(RF(norm_seeds, norms), 2));
            ynn = min( geornd(1 - (pf - weight), 1, length(norm_seeds)), sum(RB(norms, norm_seeds), 1));
            % RF(norm_seeds, norms) = RF(norm_seeds, norms) .* (rand(1, length(norms)) <= pf - weight);
            % RB(norms, norm_seeds) = RB(norms, norm_seeds) .* (rand(length(norms), 1) <= pb - weight);
        else
            xnn = min( geornd(1 - pf, length(norm_seeds), 1), sum(RF(norm_seeds, norms), 2));
            ynn = min( geornd(1 - pf, 1, length(norm_seeds)), sum(RB(norms, norm_seeds), 1));
            % RF(norm_seeds, norms) = RF(norm_seeds, norms) .* (rand(1, length(norms)) <= pf);
            % RB(norms, norm_seeds) = RB(norms, norm_seeds) .* (rand(length(norms), 1) <= pb);
        end
        for i = 1:length(norm_seeds)
            if xnn & length(find(RF(norm_seeds(i), norms)))
                randsel_xnn = randsample(find(RF(norm_seeds(i), norms)), xnn(i));
                RF(norm_seeds(i), norms) = 0;
                RF(norm_seeds(i), randsel_xnn) = 1;
            else
                RF(norm_seeds(i), norms) = 0;
            end

            if ynn & length(find(RB(norms, norm_seeds(i))))
                randsel_ynn = randsample(find(RB(norms, norm_seeds(i))), ynn(i));
                RB(norms, norm_seeds(i)) = 0;
                RB(randsel_ynn, norm_seeds(i)) = 1;
            else
                RB(norms, norm_seeds(i)) = 0;
            end
        end
        if h > 0.1   % Hg >> H(t-1)
            xns = min( geornd(1 - (pf + weight), length(norm_seeds), 1), sum(RF(norm_seeds, spams), 2));
            xsn = min( geornd(1 - (pf + weight), length(spam_seeds), 1), sum(RF(spam_seeds, norms), 2));
            yns = min( geornd(1 - (pf + weight), 1, length(spam_seeds)), sum(RB(norms, spam_seeds), 1));
            ysn = min( geornd(1 - (pf + weight), 1, length(norm_seeds)), sum(RB(spams, norm_seeds), 1));
            % RF(norm_seeds, spams) = RF(norm_seeds, spams) .* (rand(1, length(spams)) <= pf + weight);
            % RF(spam_seeds, norms) = RF(spam_seeds, norms) .* (rand(1, length(norms)) <= pf + weight);
            % RB(spams, norm_seeds) = RB(spams, norm_seeds) .* (rand(length(spams), 1) <= pb + weight);
            % RB(norms, spam_seeds) = RB(norms, spam_seeds) .* (rand(length(norms), 1) <= pb + weight);
        elseif h < -0.1
            xns = min( geornd(1 - (pf - weight), length(norm_seeds), 1), sum(RF(norm_seeds, spams), 2));
            xsn = min( geornd(1 - (pf - weight), length(spam_seeds), 1), sum(RF(spam_seeds, norms), 2));
            yns = min( geornd(1 - (pf - weight), 1, length(spam_seeds)), sum(RB(norms, spam_seeds), 1));
            ysn = min( geornd(1 - (pf - weight), 1, length(norm_seeds)), sum(RB(spams, norm_seeds), 1));
            % RF(norm_seeds, spams) = RF(norm_seeds, spams) .* (rand(1, length(spams)) <= pf - weight);
            % RF(spam_seeds, norms) = RF(spam_seeds, norms) .* (rand(1, length(norms)) <= pf - weight);
            % RB(spams, norm_seeds) = RB(spams, norm_seeds) .* (rand(length(spams), 1) <= pb - weight);
            % RB(norms, spam_seeds) = RB(norms, spam_seeds) .* (rand(length(norms), 1) <= pb - weight);
        else
            xns = min( geornd(1 - pf, length(norm_seeds), 1), sum(RF(norm_seeds, spams), 2));
            xsn = min( geornd(1 - pf, length(spam_seeds), 1), sum(RF(spam_seeds, norms), 2));
            yns = min( geornd(1 - pf, 1, length(spam_seeds)), sum(RB(norms, spam_seeds), 1));
            ysn = min( geornd(1 - pf, 1, length(norm_seeds)), sum(RB(spams, norm_seeds), 1));
            % RF(norm_seeds, spams) = RF(norm_seeds, spams) .* (rand(1, length(spams)) <= pf);
            % RF(spam_seeds, norms) = RF(spam_seeds, norms) .* (rand(1, length(norms)) <= pf);
            % RB(spams, norm_seeds) = RB(spams, norm_seeds) .* (rand(length(spams), 1) <= pb);
            % RB(norms, spam_seeds) = RB(norms, spam_seeds) .* (rand(length(norms), 1) <= pb);
        end
        for i = 1:length(norm_seeds)
            if xns & length(find(RF(norm_seeds(i), spams)))
                randsel_xns = randsample(find(RF(norm_seeds(i), spams)), xns(i));
                RF(norm_seeds(i), spams) = 0;
                RF(norm_seeds(i), randsel_xns) = 1;
            else
                RF(norm_seeds(i), spams) = 0;
            end
            if ysn & length(find(RB(spams, norm_seeds(i))))
                randsel_ysn = randsample(find(RB(spams, norm_seeds(i))), ysn(i));
                RB(spams, norm_seeds(i)) = 0;
                RB(randsel_ysn, norm_seeds(i)) = 1;
            else
                RB(spams, norm_seeds(i)) = 0;
            end
        end
        for i = 1:length(spam_seeds)
            if xsn & length(find(RF(spam_seeds(i), norms)))
                randsel_xsn = randsample(find(RF(spam_seeds(i), norms)), xsn(i));
                RF(spam_seeds(i), norms) = 0;
                RF(spam_seeds(i), randsel_xsn) = 1;
            else
                RF(spam_seeds(i), norms) = 0;
            end
            if yns & length(find(RB(norms, spam_seeds(i))))
                randsel_yns = randsample(find(RB(norms, spam_seeds(i))), yns(i));
                RB(norms, spam_seeds(i)) = 0;
                RB(randsel_yns, spam_seeds(i)) = 1;
            else
                RB(norms, spam_seeds(i)) = 0;
            end
        end
    else
        x = min(geornd(1-pf, length(seeds), 1), sum(RF(seeds,:), 2));
        y = min(geornd(1-pb, 1, length(seeds)), sum(RB(:,seeds), 1));
        % 링크 있는 것 중에 x개 골라야함
        for i = 1:length(seeds)
            if x(i) > 0
                randsel_x = randsample(find(RF(seeds(i),:)), x(i));
                RF(seeds(i),:) = 0;
                RF(seeds(i), randsel_x) = 1;
            else
                RF(seeds(i),:) = 0;
            end
            if y(i) > 0
                randsel_y = randsample(find(RB(:,seeds(i))), y(i));
                RB(:,seeds(i)) = 0;
                RB(randsel_y,seeds(i)) = 1;
            else
                RB(:,seeds(i)) = 0;
            end
        end
        % RF(seeds,:) = RF(seeds, :) .* (rand(1, len) <= pf);
        % RB(:, seeds) = RB(:, seeds) .* (rand(len, 1) <= pb);
    end

    % RF(seeds,:) = RF(seeds, :) .* (rand(1, len) <= pf);
    % RB(:, seeds) = RB(:, seeds) .* (rand(len, 1) <= pb);


    [cur_seed1 for_nexts] = find(RF);
    [back_nexts cur_seed2] = find(RB);
    A(seeds,:) = 0;
    A(:,seeds) = 0;

    fprintf(fid, '%8d %8d\n', [cur_seed1'; for_nexts']);
    fprintf(fid, '%8d %8d\n', [back_nexts'; cur_seed2']);

    for_nexts = unique(for_nexts);
    back_nexts = unique(back_nexts);
    next_seeds = [for_nexts; back_nexts];

    [os, od]=textread('output.txt','%d %d');
    sampled_nodes = unique(vertcat(os, od));
    sampled_nodes = unique(vertcat(sampled_nodes, seed_q));
    R = oA;
    R(~ismember(lins, sampled_nodes), :) = 0;
    R(:, ~ismember(lins, sampled_nodes)) = 0;
    [ix iy] = find(R);
    fprintf(fid,'%8d %8d\n', [ix'; iy']);
    fclose(fid);

    Sample = uniq_sample(k, len);
    [s, d] = find(Sample);
    [Dt0 Dt1 Ht r] = checkDH(s, d, k);
    [n1, n0, m11, m10, m00 r] = get_nms(s, d, k);

    d0 = Dg0 - Dt0;
    d1 = Dg1 - Dt1;
    h = Hg - Ht;
    fprintf('d0:%.2f d1:%.2f h:%.2f\n', d0, d1, h);

    if mod(k, max_t/20) == 0
        fprintf('saving %dth round sample...\n', max_t);
        Sample = uniq_sample(k, len);
        [s, d] = find(Sample);
        fid = fopen('./samples/num_of_spams.txt', 'a+');
        fprintf(fid, '%5d %6d %8d %7.2f %7d %6d %8d %7d %5d %5.2f\n', k, n1+n0, m11+m10+m00, (m11+m10+m00)/(n0+n1), n1, n0, m11, m10, m00, r);
        fclose(fid);
        gid = fopen('DHscores.txt', 'a');
        fprintf(gid, '%5.4f %5.4f %5.4f\n', Dt0, Dt1, Ht);
        fclose(gid);
    end

    toc;
    if length(next_seeds) == 0
        fprintf('---fire dies---\n');
        k = k + 1;
        save('round', 'k');
        return;
    end

    fprintf('\n');
    adv_FF(oA, A, Dg0, Dg1, Hg, pf, pb, length(seeds), next_seeds, max_t, k+1, d0, d1, h, wff, wfb, wbf, wbb, len, new_gene, seed_q);

end
