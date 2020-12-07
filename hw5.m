%% Problem 1
clear; close all; clc
m = 10;
A = randn(m);
A = A + A.';
[V,D] = eig(A);

% Part b
num = 1000;
BiggestEigValPerIteration = zeros(1,num);
x = randn(m,1);
for j = 1:num
    x = (A*x);
    x = x/norm(x,2);
    BiggestEigValPerIteration(1,j) = x.'*A*x;
end
mag = norm(x,2);
x = x/mag;

f4 = figure;
plot(BiggestEigValPerIteration,'r.');hold on;
plot(BiggestEigValPerIteration(end,end)*ones(size(BiggestEigValPerIteration)),'b.');

LargestEigVal = (A*x)./x;
LargestEigval = sum(LargestEigVal)/m;
PrevEigval = LargestEigval;
PrevEigVec = x;

% Part b testing around to find second largest magnitude
% eigenvalue-eigenvector
% num = 1000;
% x = randn(m,1);
% A = A-PrevEigval*eye(m);
% for j = 1:num
%     x = (A*x);
%     x = x/norm(x,2);
% end
% mag = norm(x,2);
% x = x/mag;
% 
% LargestEigVal = (A*x)./x;
% LargestEigval = sum(LargestEigVal)/m;
% LargestEigval = LargestEigval + PrevEigval;

% Part c
AllEigVals = zeros(m,1);
AllEigVecs = zeros(size(A));
NumEigsLeft = m;
index = 1;
num = 20;
BiggestEigValPerIteration2 = zeros(m,num);

while  NumEigsLeft > 0
    
z_EigValPerIter = zeros(1,num);
x = randn(m,1);
x = x/norm(x,2);

% Subtract off orthogonal eigenvectors to speed up process
componentsInEigVecDirections = (x.'*AllEigVecs);
x = x - AllEigVecs*componentsInEigVecDirections.';
x = x/norm(x,2);

lambda = 0;
for j = 1:num
    lambda = x.'*A*x;
    w = (A-lambda*eye(m))\x;
    x = w/norm(w,2);
    componentsInEigVecDirections = (x.'*AllEigVecs);
    x = x - AllEigVecs*componentsInEigVecDirections.';
    x = x/norm(x,2);
    z_EigValPerIter(1,j) = x.'*A*x;
end
% LargestEigVal = (A*x)./x;
% LargestEigval = sum(LargestEigVal)/m;
LargestEigValAlt = x.'*A*x;

diff = AllEigVals - LargestEigValAlt*ones(size(AllEigVals));
bool = abs(diff) < 10^-8;
EigvalAlreadyFound = sum(bool);
if EigvalAlreadyFound < 1
    AllEigVals(index,1) = LargestEigValAlt;
    AllEigVecs(:,index) = x;
    BiggestEigValPerIteration2(index,:) = z_EigValPerIter;
    index = index+1;
    NumEigsLeft = NumEigsLeft - 1;
end
end

f6 = figure;
plot(BiggestEigValPerIteration2(1,:)); hold on;
for l = 2:size(BiggestEigValPerIteration2,1)
    plot(BiggestEigValPerIteration2(l,:));
end

[AllEigValsCheck,I] = sort(AllEigVals);
D_Check = diag(D);
V_Check = zeros(size(V));

for k = 1:length(I)
    V_Check(:,k) = AllEigVecs(:,I(k));
end

%% Problem 1 part d for non-symmetric matrix
clear; close all; clc
m = 10;
A = randn(m)+1i*randn(m); % have to make this a complex-valued guess in order to get somewhere, otherwise it will not converge
[V,D] = eig(A);

% Part b
num = 1000;
BiggestEigValPerIteration = zeros(1,num);
x = randn(m,1);
for j = 1:num
    x = (A*x);
    x = x/norm(x,2);
    BiggestEigValPerIteration(1,j) = x'*A*x;
end
x = x/norm(x,2);

f10 = figure;
plot(BiggestEigValPerIteration);

LargestEigval = x'*A*x;
PrevEigval = LargestEigval;
PrevEigVec = x;

% Part c
noiseScalingFactor = 0.1;
InitialGuesses = V + noiseScalingFactor*randn(size(V))+1i*noiseScalingFactor*randn(size(V));
AllEigVals = zeros(m,1);
AllEigVecs = zeros(size(A));
NumEigsLeft = m;
index = 1;
num = 20;
BiggestEigValPerIteration2 = zeros(m,num);

while  NumEigsLeft > 0
    
z_EigValPerIter = zeros(1,num);
% x = randn(m,1);
% x = x/norm(x,2);
x = InitialGuesses(:,NumEigsLeft);
x = x/norm(x,2);

% Subtract off orthogonal eigenvectors to speed up process
% componentsInEigVecDirections = (x.'*AllEigVecs);
% x = x - AllEigVecs*componentsInEigVecDirections.';
% x = x/norm(x,2);

lambda = 0;
for j = 1:num
    lambda = x'*A*x;
    w = (A-lambda*eye(m))\x;
    x = w/norm(w,2);
%     componentsInEigVecDirections = (x.'*AllEigVecs);
%     x = x - AllEigVecs*componentsInEigVecDirections.';
%     x = x/norm(x,2);
    z_EigValPerIter(1,j) = x'*A*x;
end
% LargestEigVal = (A*x)./x;
% LargestEigval = sum(LargestEigVal)/m;
LargestEigValAlt = x'*A*x;

diff = AllEigVals - LargestEigValAlt*ones(size(AllEigVals));
bool = abs(diff) < 10^-8;
EigvalAlreadyFound = sum(bool);
if EigvalAlreadyFound < 1
    AllEigVals(index,1) = LargestEigValAlt;
    AllEigVecs(:,index) = x;
    BiggestEigValPerIteration2(index,:) = z_EigValPerIter;
    index = index+1;
    NumEigsLeft = NumEigsLeft - 1;
end
end

f6 = figure;
plot(BiggestEigValPerIteration2(1,:),'x'); hold on;
for l = 2:size(BiggestEigValPerIteration2,1)
    plot(BiggestEigValPerIteration2(l,:),'x');
end

[~,I] = sort(abs(AllEigVals));
AllEigValsCheck = AllEigVals(I);
D_Check = diag(D);
[~,I2] = sort(abs(D_Check));
D_Check_Sorted = D_Check(I2);
V_Check = zeros(size(V));

for k = 1:length(I)
    V_Check(:,k) = AllEigVecs(:,I(k));
end


%% Problem 2
clear; close all; clc
z_pic = imresize(imread("C:\Users\Karl\Documents\Documents\AMATH584\hw2\yalefaces_cropped\CroppedYale\yaleB01\yaleB01_P00A+000E+00.pgm","pgm"),[120 80]);

load("CorrelationMatrix")
num = 1000;
x = randn(size(Corr_Matrix,2),1);
for j = 1:num
    x = (Corr_Matrix*x);
    x = x/norm(x,2);
end
mag = norm(x,2);
x = x/mag;
LargestEigVal = x.'*Corr_Matrix*x;

load("eigenDecomp")
load("SVD")

leadingOrderSVDValue = S(1,1);
leadingOrderSVDMode = U(:,1);

% Calculate the norm of the difference between the power iterated mode and
% SVD mode
diff = norm(x-leadingOrderSVDMode,2);

% Show image version of power iterated mode vs SVD mode.
f1 = figure;
leadingOrderSVDMode=reshape(leadingOrderSVDMode,size(z_pic,1),size(z_pic,2));
subplot(1,2,1)
pcolor(flipud(leadingOrderSVDMode))
shading interp
colormap gray
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title("Dominant SVD Mode")

x=reshape(x,size(z_pic,1),size(z_pic,2));
if x(1,1) > 0
    x = -x;
end
subplot(1,2,2)
pcolor(flipud(x))
shading interp
colormap gray
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title("Dominant Power-iterated Mode")

% Random Sampling
numSamples = 50;
Omega = randn(size(Corr_Matrix,2),numSamples);
Y = Corr_Matrix*Omega;
[Q,R] = qr(Y,0);
B = Q.'*Corr_Matrix;
[U_randomSampled,S_randomSampled,V_randomSampled] = svd(B,0);
U_approxFromRandSampling = Q*U_randomSampled;
U_short = U_approxFromRandSampling(:,1:10);

% S_short = S_randomSampled(:,1:100);
% S_Short = diag(S_short);
S_Short = diag(S_randomSampled);
S_large = diag(S);
S_Large = S_large(1:numSamples,:);
f2 = figure;
plot(S_Short,'r.'); hold on;
plot(S_Large,'bo');
legend("Singular Values From Random Sampling","Singular Values From Full SVD")

f3 = figure;
semilogy(S_Short,'r.'); hold on;
semilogy(S_Large,'bo');
legend("Singular Values From Random Sampling","Singular Values From Full SVD")

f5 = figure;
n_modesPlotted = 3;
for k = 1:n_modesPlotted
    temp=reshape(U(:,k),size(z_pic,1),size(z_pic,2));
    if temp(1,1) > 0
        temp = -temp;
    end
    subplot(2,n_modesPlotted,k);
    pcolor(flipud(temp))
    shading interp
    colormap gray
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    title("True Mode " + string(k))
    
    temp=reshape(U_approxFromRandSampling(:,k),size(z_pic,1),size(z_pic,2));
    if temp(1,1) > 0
        temp = -temp;
    end
    subplot(2,n_modesPlotted,k+n_modesPlotted);
    pcolor(flipud(temp))
    shading interp
    colormap gray
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    title("Random Sampled Approx. Mode " + string(k))    
    
end

