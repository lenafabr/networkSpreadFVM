%% test weighted random selection in fortran

data = dlmread('./testout.txt');

for c = 1:10
    counts(c) = nnz(data(:,2)==c);
end
%%
list = 1:10;
plot(list,counts,'.-',list,list.^2/sum(list.^2)*size(data,1))

%% check sampling without replacement

data = dlmread('./testout.txt');

any(data(:,2)==data(:,3))

for c = 1:10
    counts(c) = nnz(data(:,2)==c);
end