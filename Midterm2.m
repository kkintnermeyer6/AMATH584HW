clear all; close all; clc;
m = 8;
n = 4;


A = randn(m,n)

[P,L,U] = LUfactor(A);
[L2,U2,P2] = lu(A);

P
P2

L
L2

U
U2



error = P*A-L*U
error_matlab = P2*A-L2*U2