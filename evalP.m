
% Code to evaluate the outputs of the submitted code. The value of the objective function
% will be written to a file. 

clear all;

%%%%%%%%%%%% Read Data Files %%%%%%%%%%%%%%%%%%%%%%

gainFile ='gains.csv';
Gains = csvread(gainFile);

paramFile = 'params.csv';
Params = csvread(paramFile);

powFile = 'pow.csv';
P = csvread(powFile);

n=Params(1); N=Params(2); theta=Params(3); c=Params(4);

% Read the partitions in to a cell named 'Part'.
partFile = 'partition.csv';
fid = fopen(partFile);
tline = fgetl(fid);
k = str2num(tline);

for i=1:k
	tline = fgetl(fid);
	Part{i} = str2num(tline);
end
fclose(fid);

%%%%%%%%%%%% Sanity checks %%%%%%%%%%%%%%%%%%%
errFile = 'err.txt';
errfid = fopen(errFile,'w');

if ~isempty(P(P<0))
	fprintf(errfid,'Negative power value\n');
	fclose(errfid);
	return;
end

Spart = [];
for i=1:k
	Spart = [Spart Part{i}];
end

if ~isequal(sort(Spart),1:n)
	fprintf(errfid,'Invalid Partitioning\n');
	fclose(errfid);
	return;
end

%%% Construct the G matrix
for i=1:k
	G{i} = -theta*Gains(Part{i}, Part{i})';
	l = length(Part{i});
	for j=1:l
		G{i}(j,j) = -1*G{i}(j,j)/theta;
	end
end

%%% Get the powers of each partition
for i=1:k
	Pow{i}=P(Part{i});
end

%%% Verify that the thresholds are met (eq.1)
for i=1:k
	Thr = round(1e3*G{i}*Pow{i}')/1e3;
	Thr = (Thr>= theta*N);
	if ~isempty(Thr(Thr==0))
		fprintf(errfid,'Threshold not met for partition %u \n',i);
		fclose(errfid);
		return;
	end
end

%%%%%%%%% Compute Objective Function %%%%%%%%%%%%%
objFile = 'obj.csv';
objF = c/k - sum(P);
dlmwrite(objFile,objF,'-append');
fclose(errfid);

