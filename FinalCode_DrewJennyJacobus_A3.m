%The purpose of this code is to act as an analyzer for an inputted data
%set file to determine the validity of the proposed truss design

clearvars %prevents potential variable interference
clc %clears command windows for optimal viewing of results

%% Step 1: Read in the data file (must be .mat format) and localize variables
    % This data file will contain the following matrices
        %Connection Matrix
        %Sx matrix of reaction forces in x-direction
        %Sy matrix of reaction forces in y-direction
        %Joint location vectors X and Y
        %Vector of applied external loads L
load('TrussDesignVALIDATION_DrewJennyJacobus_A3.mat', '-mat', 'C', 'Sx', 'Sy', 'X', 'Y', 'L') %change filename as needed to load in correct parameters
C = C;
Sx = Sx;
Sy = Sy;
X = X;
Y = Y;
L = L;

%% Step 2: Construct equilibrium equation
[joints,members] = size(C);

Ax = zeros(joints, members);
Ay = Ax;

S = [Sx; Sy];

for i = 1:joints
    index = find(C(i,:));
    for j = 1:length(index)
        second_index = find(C(:,index(j)));
        if second_index(1) ~= i
            Ax(i,index(j)) = (X(second_index(1)) - X(second_index(2))) / distance(X(second_index(1)),X(second_index(2)),Y(second_index(1)),Y(second_index(2))) ;
        else
            Ax(i,index(j)) = (X(second_index(2)) - X(second_index(1))) / distance(X(second_index(1)),X(second_index(2)),Y(second_index(1)),Y(second_index(2)));
        end
    end
end

for i = 1:joints
    index = find(C(i,:));
    for j = 1:length(index)
        second_index = find(C(:,index(j)));
        if second_index(1) ~= i
            Ay(i,index(j)) = (Y(second_index(1)) - Y(second_index(2))) / distance(X(second_index(1)),X(second_index(2)),Y(second_index(1)),Y(second_index(2))) ;
        else
            Ay(i,index(j)) = (Y(second_index(2)) - Y(second_index(1))) / distance(X(second_index(1)),X(second_index(2)),Y(second_index(1)),Y(second_index(2)));
        end
    end
end

A = [Ax ; Ay];
A = [A S];

%% Step 3: Use linear algebra to solve the unknown forces
invA = A^-1;
T = invA * L;

%% Step 4: Determine maximum theoretical load 
T2 = zeros(length(T)-3,1);
for i = 1:length(T)-3
    T2(i) = T(i);
end

L_vert = zeros(length(L)/2,1);
load_index = find(L);
vert_load_index = load_index - (length(L)/2);
L_vert(vert_load_index) = L(load_index);

R = zeros(length(T2),1);
for i = 1:length(T2)
    R(i) = T2(i) / L_vert(vert_load_index);
end
T_live = R * L_vert(vert_load_index);

w_unit_length = .0843;
T_dead = zeros(joints,members);

for i = 1:members
    index = find(C(:,i));
    dist = distance(X(index(1)),X(index(2)),Y(index(1)),Y(index(2)));
    member_dead_load = w_unit_length * dist * 0.5;
    T_dead(index(1),i) = member_dead_load;
    T_dead(index(2),i) = member_dead_load;
end

T_dead = sum(T_dead);
T_dead = T_dead';
T_dead = -T_dead;
T_total = T_live + T_dead;

P_crit = zeros(members,1);
fit_coeff = 3127.25;
member_lengths = zeros(members,1);

for i = 1:members
    index = find(C(:,i));
    dist = distance(X(index(1)),X(index(2)),Y(index(1)),Y(index(2)));
    member_lengths(i) = dist;
    P_crit(i) = - (fit_coeff / dist^2);
end
crit_total = -P_crit;
[crit_val,crit_index] = min(crit_total);

W_failure = (crit_val - T_dead(crit_index)) / R(crit_index);


%% Step 5: Determine cost of truss and load/cost ratio in oz/$
C1 = 10;
C2 = 1;

cost = C1*joints + C2*sum(member_lengths);

ratio = abs(W_failure) / cost;

%% Step 6: Display results in proper format
fprintf('\nEK301, Section A3, Group S(truss)ed Out, Drew A., Jenny D., Jacobus K., 4/2/2021\n')
fprintf('Load: %.2f oz\n', L(find(L)))
fprintf('Member forces in oz\n') %does not include dead load correction
for k = 1:length(T)-3
    if T(k) >= 0
        fprintf('m%d: %.3f (T)\n', k, abs(T(k)))
    else
        fprintf('m%d: %.3f (C)\n', k, abs(T(k)))
    end
end
fprintf('Reaction forces in oz: \n')
fprintf('Sx1: %.2f\n', T(length(T)-2))
fprintf('Sy1: %.2f\n', T(length(T)-1))
fprintf('Sy2: %.2f\n', T(length(T)))
fprintf('Cost of truss: $%.2f \n', cost)
fprintf('Theoretical max load/cost ratio in oz/$: %.4f\n', ratio)

%% Function to calculate distance between 2 points 
function d = distance(x1,x2,y1,y2)
d = sqrt((x2-x1)^2+(y2-y1)^2);
end