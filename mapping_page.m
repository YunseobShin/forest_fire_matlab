function [] = mapping_page(A, T, M, F)

global host;
global flag;
m = size(A,1);
a = cell(m,1);
b = cell(m,1);
parfor node=1:m
    %fprintf('[mapping_page] %3.3f%%.\n',node/m*100);
    [ti,tj] = find(T==node,1);
    if tj==1 % source
        a{node} = [node M(ti,3)];
        if F(ti,1) == '-'
            b{node} = [node 0]; %normal
        else
            b{node} = [node 1]; %spam
        end
    elseif tj==2 % destination
        a{node}= [node M(ti,4)];
        if F(ti,2) == '-'
            b{node} = [node 0];
        else
            b{node} = [node 1];
        end
    else
        error('invalid index');
    end
end

DV = cell2mat(a);
host(DV(:,1)) = DV(:,2);
NS = cell2mat(b);
flag(NS(:,1)) = NS(:,2);

fprintf('\n---Docs Level---\n');
fprintf('[1] *** No. of normal docs: %d (%3.2f%%)\n',nnz(flag==0),nnz(flag==0)/m*100);
fprintf('[2] *** No. of spam docs: %d (%3.2f%%)\n',nnz(flag==1),nnz(flag==1)/m*100);
