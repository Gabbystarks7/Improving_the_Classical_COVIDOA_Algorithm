%_______________________________________________________________________________________
%  The Coronavirus Disease Optimization Algorithm (COVIDOA) source codes demo version 1.0                  
%                                                                                       
%  Developed in MATLAB R2016a                                                     
%                                                                                       
%  Authors: Asmaa M. Khalid, Khalid M. Hosny, & Seyedali Mirjalili                      
%                                                                                       
%  E-Mail: asmaa.elhenawy@gmail.com  (Asmaa Khalid)                                               
%  Homepage: https://www.researchgate.net/profile/Asmaa-Khalid-3                      
%                                                                                      
% Main paper:   COVIDOA: a novel evolutionary optimization algorithm based on coronavirus disease replication lifecycle
% Reference: Khalid, A. M., Hosny, K. M., & Mirjalili, S. (2022). COVIDOA: a novel evolutionary optimization algorithm based on coronavirus disease replication lifecycle. Neural Computing and Applications, 1-28.?
%
%_______________________________________________________________________________________

clear all 
clc
MaxIt=500;
nPop=300;
minVal=10;
maxVal=50;
D=30;
CostFunction=@F3;
shifttingNo=1;
numOfSubprotiens=6;
MR=0.1;

[Best_FF,Best_P,Conv_curve]=COVIDOA(nPop,MaxIt,minVal,maxVal,D,CostFunction,MR,shifttingNo,numOfSubprotiens);

subplot(1,2,2);
semilogy(Conv_curve,'Color','r','LineWidth',2)
title('Convergence curve')
xlabel('Iterations');
ylabel('Best Cost');
axis tight
legend('Improved COVIDOA')