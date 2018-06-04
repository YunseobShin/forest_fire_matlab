function [n1, n0, m11, m10, m00 r] = get_nms(s, d, k)
    nodes = vertcat(s, d);
    uq_nodes = unique(nodes);

    load flag.mat
    flag(find(flag)) = 2;
    flag(find(flag==0)) = 1;
    flag(find(flag==2)) = 0;

    n1 = length(find(flag(uq_nodes)));
    n0 = length(find(flag(uq_nodes)==0));
    fs = flag(s);
    fd = flag(d);

    % load flag2.mat
    % n1 = length(find(flag2(uq_nodes)));
    % n0 = length(find(flag2(uq_nodes)==0));
    % fs = flag2(s);
    % fd = flag2(d);

    spams = fs+fd;
    % fid = fopen('./samples/num_of_spams.txt', 'a+');

    m00 = length(find(spams == 0));
    m10 = length(find(spams == 1));
    m11 = length(find(spams == 2));
    r = n0/(n1+n0)*100;

    % fprintf(fid, '%5d %6d %8d %7.2f %7d %6d %8d %7d %5d %5.2f\n', k, n1+n0, m11+m10+m00, (m11+m10+m00)/(n0+n1), n1, n0, m11, m10, m00, r);
    % fclose(fid);
end
