function adv_FF_rate(oA, A, rate, pf, pb, prev_length, seeds, max_t, k, r, wff, wfb, wbf, wbb, len, new_gene, seed_q, nol)
    % seeds
    fid = fopen('output.txt', 'a');
    if k > max_t
        k = max_t+1;
        save('round', 'k');
        return;
    end
    if pf - wfb < 0
        k = max_t+1;
        save('round', 'k');
        disp('pf - wfb < 0');
        return;
    end
    if pf + wff > 1
        k = max_t+1;
        save('round', 'k');
        disp('pf + wff > 1');
        return;
    end
    if pb - wbb < 0
        k = max_t+1;
        save('round', 'k');
        disp('pb - wbb < 0');
        return;
    end
    if pb + wbf > 1
        k = max_t+1;
        save('round', 'k');
        disp('pb + wbf < 1');
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
    fprintf('current time: %d:%d:%d\n', ct(4), ct(5), ct(6));
    fprintf('# of seeds: %d\n', length(seeds));

    len = size(A,1);
    next_seeds=[];

    RF = A;
    RB = A;
    lins = linspace(1,len, len);
    RF(~ismember(lins, seeds), :) = 0;
    RB(:, ~ismember(lins, seeds)) = 0;

    if k > 1
        if r > rate + 0.1
            x_spams = min(geornd(1 - (pf - wfb), length(seeds), 1), sum(RF(seeds,spams), 2));
            y_spams = min(geornd(1 - (pb - wbb), 1, length(seeds)), sum(RB(spams,seeds), 1));
            x_norms = min(geornd(1 - (pf + wff), length(seeds), 1), sum(RF(seeds,norms), 2));
            y_norms = min(geornd(1 - (pb + wbf), 1, length(seeds)), sum(RB(norms,seeds), 1));
            for i = 1:length(seeds)
                if x_spams(i) & length(find(RF(seeds(i), spams)))
                    randsel_xs = randsample(find(RF(seeds(i), spams)), x_spams(i));
                    RF(seeds(i), spams) = 0;
                    RF(seeds(i), randsel_xs) = 1;
                else
                    RF(seeds(i), spams) = 0;
                end
                if y_spams(i) & length(find(RB(spams, seeds(i))))
                    randsel_ys = randsample(find(RB(spams, seeds(i))), y_spams(i));
                    RB(spams, seeds(i)) = 0;
                    RB(randsel_ys, seeds(i)) = 1;
                else
                    RB(spams, seeds(i)) = 0;
                end
                if x_norms(i) & length(find(RF(seeds(i), norms)))
                    randsel_xn = randsample(find(RF(seeds(i), norms)), x_norms(i));
                    % fprintf('xn: %d\n', length(randsel_xn));
                    RF(seeds(i), norms) = 0;
                    RF(seeds(i), randsel_xn) = 1;
                    % randsel_xn
                else
                    RF(seeds(i), norms) = 0;
                end
                if y_norms(i) & length(find(RB(norms, seeds(i))))
                    randsel_yn = randsample(find(RB(norms, seeds(i))), y_norms(i));
                    % fprintf('yn: %d\n', length(randsel_yn));
                    RB(norms, seeds(i)) = 0;
                    RB(randsel_yn, seeds(i)) = 1;
                    % randsel_yn
                else
                    RB(norms, seeds(i)) = 0;
                end
            end
            % RF(seeds, spams) = RF(seeds, spams) .* (rand(1, length(spams)) <= pf);
            % RB(spams, seeds) = RB(spams, seeds) .* (rand(length(spams), 1) <= pb);
            % RF(seeds, norms) = RF(seeds, norms) .* (rand(1, length(norms)) <= pf + weight);
            % RB(norms, seeds) = RB(norms, seeds) .* (rand(length(norms), 1) <= pb + weight);
        elseif r < rate - 0.1
            x_spams = min(geornd(1 - (pf + wff), length(seeds), 1), sum(RF(seeds,spams), 2));
            y_spams = min(geornd(1 - (pb + wbf), 1, length(seeds)), sum(RB(spams,seeds), 1));
            x_norms = min(geornd(1 - (pf - wfb), length(seeds), 1), sum(RF(seeds,norms), 2));
            y_norms = min(geornd(1 - (pb - wbb), 1, length(seeds)), sum(RB(norms,seeds), 1));
            for i = 1:length(seeds)
                if x_spams(i) & length(find(RF(seeds(i), spams)))
                    randsel_xs = randsample(find(RF(seeds(i), spams)), x_spams(i));
                    % fprintf('xs: %d\n', length(randsel_xs));
                    RF(seeds(i), spams) = 0;
                    RF(seeds(i), randsel_xs) = 1;
                else
                    RF(seeds(i), spams) = 0;
                end
                if y_spams(i) & length(find(RB(spams, seeds(i))))
                    randsel_ys = randsample(find(RB(spams, seeds(i))), y_spams(i));
                    % fprintf('ys: %d\n', length(randsel_ys));
                    RB(spams, seeds(i)) = 0;
                    RB(randsel_ys, seeds(i)) = 1;
                else
                    RB(spams, seeds(i)) = 0;
                end
                if x_norms(i) & length(find(RF(seeds(i), norms)))
                    randsel_xn = randsample(find(RF(seeds(i), norms)), x_norms(i));
                    % fprintf('xn: %d\n', length(randsel_xn));
                    RF(seeds(i), norms) = 0;
                    RF(seeds(i), randsel_xn) = 1;
                    % randsel_xn
                else
                    RF(seeds(i), norms) = 0;
                end
                if y_norms(i) & length(find(RB(norms, seeds(i))))
                    randsel_yn = randsample(find(RB(norms, seeds(i))), y_norms(i));
                    % fprintf('yn: %d\n', length(randsel_yn));
                    RB(norms, seeds(i)) = 0;
                    RB(randsel_yn, seeds(i)) = 1;
                    % randsel_yn
                else
                    RB(norms, seeds(i)) = 0;
                end
            end
            % RF(seeds, spams) = RF(seeds, spams) .* (rand(1, length(spams)) <= pf + weight);
            % RB(spams, seeds) = RB(spams, seeds) .* (rand(length(spams), 1) <= pb + weight);
            % RF(seeds, norms) = RF(seeds, norms) .* (rand(1, length(norms)) <= pf - weight);
            % RB(norms, seeds) = RB(norms, seeds) .* (rand(length(norms), 1) <= pb - weight);
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
            % RF(seeds, spams) = RF(seeds, spams) .* (rand(1, length(spams)) <= pf);
            % RB(spams, seeds) = RB(spams, seeds) .* (rand(length(spams), 1) <= pb);
            % RF(seeds, norms) = RF(seeds, norms) .* (rand(1, length(norms)) <= pf);
            % RB(norms, seeds) = RB(norms, seeds) .* (rand(length(norms), 1) <= pb);
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

    % x = min(geornd(1-pf, length(seeds), 1), sum(RF(seeds,:), 2));
    % y = min(geornd(1-pb, 1, length(seeds)), sum(RB(:,seeds), 1));
    % for i = 1:length(seeds)
    %     if x(i) > 0
    %         randsel_x = randsample(find(RF(seeds(i),:)), x(i));
    %         RF(seeds(i),:) = 0;
    %         RF(seeds(i), randsel_x) = 1;
    %     else
    %         RF(seeds(i),:) = 0;
    %     end
    %     if y(i) > 0
    %         randsel_y = randsample(find(RB(:,seeds(i))), y(i));
    %         RB(:,seeds(i)) = 0;
    %         RB(randsel_y,seeds(i)) = 1;
    %     else
    %         RB(:,seeds(i)) = 0;
    %     end
    % end

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
    next_seeds = unique([for_nexts; back_nexts]);

    sampled_nodes = unique(seed_q);
    R = oA;
    R(~ismember(lins, sampled_nodes), :) = 0;
    R(:, ~ismember(lins, sampled_nodes)) = 0;
    [ix iy] = find(R);
    nol = nol + length([cur_seed2; cur_seed1]);

    fprintf('induced: %d FF:%d\n', length(ix), nol);
    fprintf(fid, '%8d %8d\n', [ix'; iy']);
    fclose(fid);
    Sample = uniq_sample(k, len);
    [s, d] = find(Sample);
    fid = fopen('output.txt', 'w');
    fprintf(fid, '%8d %8d\n', [s'; d']);
    fclose(fid);
    [n1, n0, m11, m10, m00 r] = get_nms(s, d, k);

    if mod(k, max_t/20) == 0
        fprintf('saving %dth round sample...\n', max_t);
        Sample = uniq_sample(k, len);
        [s, d] = find(Sample);
        [Dt0 Dt1 Ht r] = checkDH(s, d, k);
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
        save('nol', 'nol');
        return;
    end
    fprintf('\n');
    adv_FF_rate(oA, A, rate, pf, pb, length(seeds), next_seeds, max_t, k+1, r, wff, wfb, wbf, wbb, len, new_gene, seed_q, nol);

end
