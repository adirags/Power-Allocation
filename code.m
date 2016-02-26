clear
clc
g = csvread('gains.csv');
params = csvread('params.csv');
n = params(1);        %number of users
N = params(2);      %Noise Power
theta = params(3);  %SNR threshold for decoding
C = params(4);      %Constant
G = zeros(n,n);
for i=1:n
    for j=1:n
        if(i==j)
            G(i,j) = g(i,j);
        else
            G(i,j) = -1*theta*g(j,i);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Need to maximize the objective function given by                  %
%           (C/k)-sum(P), where k is the number of partitions             %
%       subject to the constraint that:                                   %
%           P_S_i >= inv(G_S_i)*theta*N*1                                 %
%*************************************************************************%
% We need to perform joint power allocation and link scheduling. This can %
% be done in the following ways:                                          %
%       1. Perform link scheduling, and then find optimal power for each  %
%          set of active links.                                           %       
%       2. Perform joint link scheduling and power allocation.            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Attempting method wherein both link scheduling and power allocation     %
% is done alternatively.                                                  %
% The algorithm proceeds as follows :                                     %
%       1. Initialize A which is the set of active links as {}, and       %
%          initialize last_set = 1                                        %
%       2. Compute Z, which is the sum of the gains for each pair of      %
%          Tx-Rx. Z = { g(i,i)          , if i=j                          %
%                     { g(i,j)+g(j,i)   , if i!=j                         %
%       3. Consider the min(sum(Z)) which would give us the wireless node %
%          that has the least effect on all other nodes. Add this node to %
%          temp_A{last_set}                                               %
%       4. Compute Power Allocation for this subset. If all powers are    %
%          positive, add this node to A{last_set} and make                %
%          temp_A{last_set = A{last_set}, else last_set = last_set+1,     % 
%          and add this node to A{last_set}.                              %
%       5. Set the entire column of Z corresponding to the node above to  %
%          inf.                                                           %
%       6. Repeat from Step 3, till all nodes assigned to a set.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = zeros(n,n);
for i=1:n
    for j=1:n
        if(i==j)
            Z(i,j) = g(i,j);
        else
            Z(i,j) = g(i,j) + g(j,i);
        end
    end
end
A = {};
temp_A = A;
last_set = 1;
B_factor = inf;     %Factor which tells us how the activity of one node affects the others.
sum_Z = sum(Z);
num_allocated_nodes = 0;
loc = 0;
for i=1:n
    if(B_factor>sum_Z(i))
        B_factor = sum_Z(i);
        loc = i;
    end
end
temp_A{last_set} = loc;
S = loc; %use [loc] if error
G1 = G(loc,loc);
P1 = inv(G1)*theta*N;
P{last_set} = P1; %use [P1] if error
num_allocated_nodes = num_allocated_nodes + 1;
Z(:,loc) = inf;
A{last_set} = temp_A{last_set};
temp_loc = loc;
while(num_allocated_nodes<=n)
    B_factor = inf;
    sum_Z = sum(Z);
    poor_choice_flag = 0;
    for i=1:n
        if(sum_Z(i) ~= inf)
            if(B_factor>sum_Z(i))
                B_factor = sum_Z(i);
                loc = i;
            end
        end
    end
    if(loc ~= temp_loc)
        temp_A{last_set} = [temp_A{last_set},loc];
        S = [S,loc];
    end
    G1 = G(S,S);
    [~,b] = size(S);
    if(cond(G1)<10000)
        P1 = inv(G1)*theta*N*ones(b,1);
    else
        poor_choice_flag = 1;
    end
    for i=1:b
        if(P1(i)<=0)
            poor_choice_flag = 1;
        end
    end
    if(poor_choice_flag == 0)
        A{last_set} = temp_A{last_set};
        P{last_set} = P1;
        Z(:,loc) = inf;
        num_allocated_nodes = num_allocated_nodes + 1;
        temp_loc = loc;
    else
        last_set = last_set + 1;
        temp_A{last_set} = loc;
        S = loc; %use [loc] if error
        temp_loc = loc;
    end
end
[~,b] = size(A);
k = b;          %Number of Time Divided sets
P_sum = 0;
for i=1:3
    P1 = P{i};
    [a,~] = size(P1);
    for j=1:a
        P_sum = P_sum + P1(j);
    end
end
objFunc_method1 = C/k - P_sum;
P_method1 = P;
A_method1 = A;
k_method1 = k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Method :                                                         %
% Find the link or link combination that has the lowest value of Z.       %
% Select it and add the corresponding nodes to temp_A, Find power         %
% allocation for this set. If it is feasible, add to A, else create new   %
% partition.                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n
    for j=1:n
        if(i==j)
            Z(i,j) = g(i,j);
        else
            Z(i,j) = g(i,j) + g(j,i);
        end
    end
end

temp_A = {};
A2 = {};
P = {};
last_set = 1;
S = [];
min = inf;
for i=1:n
    for j=1:n
        if(min>Z(i,j))
            min = Z(i,j);
            loc_i = i;
            loc_j = j;
        end
    end
end
if(loc_i~=loc_j)
    S = [loc_i, loc_j];
    temp_A{last_set} = [loc_i, loc_j];
else
    S = loc_i;
    temp_A{last_set} = loc_i;
end
G1 = G(S,S);
[a,b] = size(S);
P1 = inv(G1)*theta*N*ones(b,1);
P2{last_set} = P1;
A2{last_set} = temp_A{last_set};
num_allocated_nodes = b;
Z(loc_i,:) = inf;
Z(:,loc_j) = inf;
Z(loc_j,loc_i) = inf;

while(num_allocated_nodes <= n)
    min = inf;
    poor_choice_flag = 0;
    for i=1:n
        for j=1:n
            if(min>Z(i,j))
                min = Z(i,j);
                loc_i = i;
                loc_j = j;
            end
        end
    end
    [~,b] = size(S);
    initial_size = b;
    a = 0;
    for i=1:b
        if(loc_i~=S(i)) %Checking if newly identified node to be added isn't already there
            a = a+1;
        end
    end
    b1 = 0;
    a1 = 0;
    if(last_set>1)
        [~,b1] = size(A2{last_set-1});
        S1 = A2{last_set-1};
        for i=1:b1
            if(loc_i~=S1(i))
                a1 = a1+1;
            end
        end
    end
    if(a==b)&&(a1==b1)
        S = [S,loc_i];
    end
    a = 0;
    for i=1:b
        if(loc_j~=S(i))
            a = a+1;
        end
    end
    b1 = 0;
    a1 = 0;
    if(last_set>1)
        [~,b1] = size(A2{last_set-1});
        S1 = A2{last_set-1};
        for j=1:b1
            if(loc_j~=S1(j))
                a1 = a1+1;
            end
        end
    end
    if(a==b)&&(a1==b1)
        if(loc_i~=loc_j)
            S = [S,loc_j];
        end
    end
    G1 = G(S,S);
    [a,b] = size(S);
    temp_A{last_set} = S;
    if(cond(G1)<100000)
        P1 = inv(G1)*theta*N*ones(b,1);
    else
        poor_choice_flag = 1;
    end
    [a,b] = size(P1);
    for i=1:b
        if(P1(i)<=0)
            poor_choice_flag = 1;
        end
    end
    if(poor_choice_flag == 0)
        A2{last_set} = temp_A{last_set};
        P2{last_set} = P1;
        num_allocated_nodes = num_allocated_nodes + b;
        Z(loc_i,:) = inf;
        Z(:,loc_j) = inf;
        Z(loc_j,loc_i) = inf;
    else
        last_set = last_set + 1;
        tempS = S;
        S = [];
        if(loc_i==loc_j)
            temp_A{last_set} = loc_i;
            S = loc_i;
        else
            a1 = 0;
            a2 = 0;
            b = initial_size;
            for i=1:b
                if(loc_i ~= tempS(i))
                    a1 = a1+1;
                end
                if(loc_j ~= tempS(i))
                    a2 = a2+1;
                end
            end
            if (a1 == b)&&(a2==b)
                S = [loc_i,loc_j];
                temp_A{last_set} = [loc_i,loc_j];              
            elseif (a1==b)
                S = [S,loc_i];
                temp_A{last_set} = [loc_i];
            elseif (a2 == b)
                S = [S,loc_j];
                temp_A{last_set} = [loc_j];
            end
        end
    end
end
[~,b] = size(A2);   %b gives us the number of TDM slots
k = b;
P_sum = 0;
for i=1:b
    P1 = P2{i};
    [x,y] = size(P1);
    for j=1:x
        P_sum = P_sum + P1(j);
    end
end
objFunc_method2 = C/k - P_sum;
A_method2 = A2;
P_method2 = P2;
k_method2 = k;

% Choosing optimal case

if(objFunc_method1>objFunc_method2)
%     Method 1 gave a better result
    A_final = A_method1;
    P_final = P_method1;
    k_final = k_method1;
    objFunc_final = objFunc_method1;
else
    A_final = A_method2;
    P_final = P_method2;
    k_final = k_method2;
    objFunc_final = objFunc_method2;
end
display(objFunc_method1);
display(objFunc_method2);
display('Optimal value of objective function');
display(objFunc_final);
[~,b] = size(P_final);
P_write = zeros(1,n);
for i=1:b
    a = A_final{i};
    p = P_final{i};
    [~,y] = size(a);
    for j=1:y
        P_write(a(j)) = p(j);
    end
end
csvwrite('pow.csv',P_write);
csvwrite('partition.csv',k_final);
for i=1:b
    a = A_final{i};
    dlmwrite('partition.csv',a,'-append');
end
