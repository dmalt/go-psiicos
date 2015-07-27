%% drawlcurve

N1 = [];
N2 = [];
cd ./lcurve
out = dir('./Output_*');
numfiles = length(out);
mydata = cell(1, numfiles);
for k=1:numfiles
	mydata{k} = load(out(k).name);
	N1 = [N1,mydata{k}.N1];
	N2 = [N2,mydata{k}.N2];
end
figure; plot(N1,N2, 'c.');