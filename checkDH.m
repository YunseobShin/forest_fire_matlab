function [D0 D1 H r] = checkDH(s, d, k)
    node = vertcat(s,d);
    N =length(unique(node));
    M = length(s);
    p = (2*M)/(N*(N-1));
    [n1, n0, m11, m10, m00 r] = get_nms(s, d, k);
    fprintf('normal: %d, spam: %d\n', n1, n0);

    if n1 > 1
        e_m11 = nchoosek(n1, 2) * p;
    else
        e_m11 = n1;
    end
    if n0 > 1
        e_m00 = nchoosek(n0, 2) * p;
    else
        e_m00 = n0;
    end
    if n1*n0 > 0
        e_m10 = n1 * n0 * p;
    else
        e_m10 = 0;
    end

    if e_m00 > 0
        D0 = m00/e_m00;
    else
        D0=0;
    end
    if e_m11 > 0
        D1 = m11/e_m11;
    else
        D1 = 0;
    end
    if e_m10 > 0
        H = m10/e_m10;
    else
        H = 0;
    end

end
