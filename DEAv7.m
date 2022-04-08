%% This code contains a MATLAB code which uses the stochastic input oriented data envelopment analysis model to evaluate the efficiency scores 
%% of different decision making units. Together with the deterministic variables, this code can handle two type of stochastic input variables, 
%% one is following normal distribution and the other is following Poisson distribution.   
clear;
%%%%%%%%%%%%%%%%%%%%%%Given
InputD=[6,15,10];% %Deterministic inputs
Output=[29,27,2;23,8,12];%
mean=[17,28,12];%mean normal
covM=[1.3,0.8,0.7;0.8,1.4,0.6;0.7,0.6,1.3];
meanP=[13,18,16];%mean Poisson
disp(meanP);
covMP=diag(meanP);
epsilon=0.3;
e=norminv(epsilon);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ID=fopen('resultDEA.doc','w');%the results will be written in this file
n=size(Output,2);%number of DMUs
r=size(Output,1);% number of outputs
m_D=size(InputD,1);%number of deterministic inputs
Eff=zeros(n,3);
%Lambdas are between 0 and 1
lb=zeros(n+1,1);
ub=ones(n+1,1);
%Sum of lambdas equal one
Aeq=[ones(1,n),0];
Beq=1;
for np=0:1 
if(np==0)
    fprintf(ID,'Normal distribution\n');
else
    fprintf(ID,'Poisson distribution\n');
end
for p=1:n
%Input linear constraint
X=zeros(n+1,1);% X(n+1)=theta and X(i)=lambda_i
A=[InputD,-InputD(p)];
%Output linear constraint
A(m_D+1:m_D+r,:)=-[Output,zeros(r,1)];
B(m_D+1:m_D+r,1)=-Output(:,p);
options=optimoptions('fmincon','Algorithm','sqp');
X=fmincon(@objfun,X,A,B,Aeq,Beq,lb,ub,@(X)confun(X,p,e,epsilon,mean,covM,meanP,covMP,np),options);
%disp(X);
Eff(p,:)=[p,X(n+1),X(p)];
end
for p=1:n
fprintf(ID,'p=%i, theta =%f,lambda=%f\n',Eff(p,:));
end
fprintf(ID,'\n');
end
fclose(ID);
function f =objfun(X)
n=length(X)-1;
f=X(n+1,1);
end
function [c,ceq]=confun(X,p,e,epsilon,mean,covM,meanP,covMP,np)
% Nonlinear inequality constraints
n=length(X)-1;
alpha(:,1)=X(1:n,1);
nu=abs(meanP*alpha-meanP(p)*alpha(p,1));%!!!
alpha(p,1)=alpha(p,1)-X(n+1,1);%\labda_p-\theta_p
c(1)=mean*alpha-e*sqrt(transpose(alpha)*covM*alpha);
mu=abs(-meanP(p)*alpha(p,1));%!!!
if(np==0)
    c(2)=meanP*alpha-e*sqrt(transpose(alpha)*covMP*alpha);
elseif(np==1)
    c(2)=marcumq(sqrt(2*mu),0)-marcumq(sqrt(2*mu),sqrt(2*nu))-epsilon;
end
%Nonlinear equality constraints
ceq=[];
end
