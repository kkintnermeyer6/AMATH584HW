clear all; close all; clc;
c = newline;


%% Question 1 ----------------------------------------------------
disp("Problem 1 Output")
for p = 1:3
   disp("Iteration number: "+string(p)+c)
   A = 100*randn(20,7);
   if p == 2
       A(:,end) = A(:,1);
   end
   disp("Condtion number of Matrix: ")
   cond(A)
   [Q,R] = factorToQR(A)
   [Q_matlab,R_matlab] = qr(A,0)
   [Q_other,R_other] = qrfactor(A)
%    Q_other = Q_other(:,1:10);
%    R_other = R_other(1:10,:);
    disp("Sanity Check to ensure Q*R matrices equal A Matrix")
    disp("Diff indicates A-Q*R")
    Diff = A-Q*R
    Diff_matlab = A-Q_matlab*R_matlab;
    Diff_other = A-Q_other*R_other;
    disp("Sanity Check to ensure Q Matrix is orthogonal")
    disp("Take Q'*Q to ensure this equals Identity")
    OrthogonalityCheck = Q'*Q
    disp(""+c+c)
end



% ----------------------------------------------------------------


%% Question 2 ----------------------------------------------------
disp(""+c+c+"Problem 2 output")
delta_x = 0.001;
x_min = 1.920;
x_max = 2.080;

f_LHS = @(x) (x-2)^9;
f_RHS = @(x) (x^9)-18*(x^8)+144*(x^7)-672*(x^6)+2016*(x^5)-4032*(x^4)+5376*(x^3)-4608*(x^2)+2304*(x)-512;


index = 0;
Values = zeros(3,161);

tic
for k = x_min:delta_x:x_max
    index = index+1;
    Values(1,index) = k;
    Values(2,index) = f_RHS(k);
    Values(3,index) = f_LHS(k);
    
end
toc


f1 = figure;
plot(Values(1,:),Values(2,:));
hold on;
plot(Values(1,:),Values(3,:))
title("plots of function evaluated with RHS and LHS expressions")
xlabel('x')
ylabel('f(x)')


% Vectorized calc -------------------------------
y = [x_min:delta_x:x_max]';
tic
soln_lhs = (y-2).^9;

soln_rhs = y.^9-18*y.^8+144*y.^7-672*y.^6+2016*y.^5-4032*y.^4+5376*y.^3-4608*y.^2+2304*y-512;
toc
y = y';
soln_rhs = soln_rhs';
soln_lhs = soln_lhs';


plot(y,soln_rhs)
plot(y,soln_lhs)

legend('RHS','LHS','RHS Vectorized','LHS Vectorized')

% ----------------------------------------------

% ----------------------------------------------------------------


%% Question 3 -----------------------------------------------------
disp(""+c+c+"Problem 3 output")
m = 5;
n = 3;
numTimes = 60;

CondNumBySize = zeros(3,numTimes);
% Part a -------------------
for k = 1:numTimes
   A = randn(m+10*(k-1),n+10*(k-1));
   CondNumBySize(1,k) = m+10*(k-1);
   CondNumBySize(2,k) = n+10*(k-1);
   CondNumBySize(3,k) = cond(A);
   
end
numCells = CondNumBySize(1,:).*CondNumBySize(2,:);
f3 = figure;
plot(numCells,CondNumBySize(3,:))
title("Condition Number vs Size of Matrix")
xlabel("Number of Cells (Rows x Columns)")
ylabel("Condition Number")
% --------------------------


% Part b -------------------
A = randn(m+100,n+100-2);
A = [A,A(:,1)];
disp("" + c + c)
disp("Condition number of A with last column identical to first column")
Condtion_A = cond(A)

% Append the column three more times to make matrix square
A_prime = [A,A(:,1),A(:,1),A(:,1)];
disp(c+ "Determinant of A with last 4 columns identical to first column")
Determinant_A = det(A_prime)
% --------------------------

numTimes2 = 10;
CondNumByNoise = zeros(2,numTimes2);
noise = randn(size(A,1),1);
factor = 10^-9;
% Part c -------------------
for k = 1:numTimes2
    A_Prime = A;
    A_Prime(:,end) = A_Prime(:,end) + factor*(k-1)*noise;
    CondNumByNoise(1,k) = factor*(k-1);
    CondNumByNoise(2,k) = cond(A_Prime);
    
end

f2 = figure;
semilogy(CondNumByNoise(1,:),CondNumByNoise(2,:))
title("Condition Number as a function of amount of noise added")
ylabel("Condition Number")
xlabel('$\epsilon$','Interpreter','latex','fontsize',14)

% --------------------------


% -----------------------------------------------------------------