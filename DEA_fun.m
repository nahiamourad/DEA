%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The following can be used to read the data from Excel file, indicate
%%% the parameters for the choice of the DEA model and write the results in
%%% the same Excel file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% file_name='DEA_test.xlsx';%%% Wite the name of the file containing the data
% 
% %%Read the data from the given file
% sheet_name=sheetnames(file_name);
% num_sheet=max(sheet_name);
% data=cell(num_sheet);
% for k=1:num_sheet
%     data{k}=readmatrix(file_name, 'Sheet', sheet_name(k));
%     data{k} = data{k}(:,all(~isnan(data{k}))); % to remove nan - columns
% end
% 
% %% Determining the inputs and output
% Data=data{4}; %% the data to be used are in sheet number 4
% I_D=2:4;    %% the columns for deterministic input variables
% O_D=5:6;    %% the columns for deterministic output variables
% VRS=1;      %% indicates the returns orientation: VRS=0 for CRS model, VRS=1 for VRS model
% Ort_O=0;    %% indicates the orientation: Ort_O=0 for input orientaion, Ort_O=1 for output orientaion
% 
% epsilon=0.3;%% epsilon for chance constrained if exists
% I_S=[];     %% the columns for stochastic input variables
% O_S=[];     %% the columns for stochastic output variables
% covM_I=data{5}; %% Covariance matrix for stochastic input varaibles are on sheet number 5;
% % if there are more than one stochastic variable the covariance matrices should be placed
% % next to each other in the same excel sheet
% covM_O=[];  %%%% Covariance matrix for stochastic output varaible
% 
% MI=1;       %% MI=1 indicates that the Malmquist index is to be used
% I_D2=8:10;  %% if MI=1, I_D2 is the columns for deterministic input variables for the second period
% O_D2=11:12; %% if MI=1, O_D2 is the columns for deterministic output variables for the second period
% 
% %% Call the DEA function
% [Scores,tfpch,catch_up,frontier_shift]=DEA_fun(VRS,Ort_O,Data,I_D,O_D,I_S,O_S,covM_I,covM_O,epsilon,MI,I_D2,O_D2);
% 
% %% Write the results
% if(MI==1)
%     writematrix([catch_up,frontier_shift,tfpch],file_name,'Sheet','result','Range','A2');
% else
%     writematrix(Scores,file_name,'Sheet','result','Range','A2');
% end


function [Scores,tfpch,catch_up,frontier_shift]=DEA_fun(VRS,Ort_O,Data,I_D,O_D,I_S,O_S,covM_I,covM_O,epsilon,MI,I_D2,O_D2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% This function calculates the deterministic/stochastic DEA efficiency
%%%%%% scores for VRS/CRS with input/output orientations.
%%%%%% Moreover, it calculates the tfpch for Malmquisit index based on DEA
%%%%%% models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The variables are defined as follows
% Scores: efficiency scores
% tfpch: total factor productivity change in case Malmquist index is used
% catch_up: the catch up term assochiated to the Malmquisit Index DEA
% forntier_shift: the forntier shift term assochiated to the Malmquisit Index DEA
% VRS: parameter indicating the returns to scale (VRS=0 for CRS model,
% VRS=1 for VRS model)
% Ort_O: parameter indicating orientation (Ort_O=0 for input orientaion,
% Ort_O=1 for output orientaion)
% Data: Matrix containing the data to be used
% I_D: the number of columns in the Data matrix for the deterministic input
% O_D: the number of columns in the Data matrix for the deterministic output
% I_S: the number of columns in the Data matrix for the stochastic input
% represnted by their means
% O_S: the number of columns in the Data matrix for the stochastic output
% covM_I: Covariance matrix for stochastic input varaibles,
% if there are more than one stochastic variable the covariance matrices
% should be placed next to each other in one matrix
% covM_O: Covariance matrix for stochastic output varaible
% epsilon: epsilon for the chance constraint of the stochastic variable
% MI: parameter for the Malmquist index (MI=1 for MI-DEA model)
% I_D2: the columns for deterministic input variables for the second
% period, if MI=1
% O_D2: the columns for deterministic output variables for the second
% period, if MI=1


%%%% Define the matrix of deterministic variables
per=1;
Input(:,:,1)=Data(:,I_D)';
Output(:,:,1)=Data(:,O_D)';
if(MI==1)
    Input(:,:,2)=Data(:,I_D2)';
    Output(:,:,2)=Data(:,O_D2)';
    per=2;
end

%%% Define the mean of stochastic variables
if(~isempty(I_S)||~isempty(O_S))
    e=norminv(epsilon);
    mean_I=Data(:,I_S)';
    mean_O=Data(:,O_S)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(Output,2);%number of DMUs
m_O=size(Output,1);% number of deterministic outputs
m_I=size(Input,1);%number of deterministic inputs
Eff=zeros(N,3);

delta=1-kronecker_delta(1,Ort_O);

%%%% upper and lowerbounds on \lambda and \theta
lb=[zeros(N,1);Ort_O];
ub=[Inf(N,1);dirac_delta(1-Ort_O)+1];

%%%% Equality constrained
Aeq=VRS*[ones(1,N),0];
beq=VRS;

%%%% Objective function used as a vector for the linear model
objvec=[zeros(N,1);(-1)^Ort_O];%(-1)^O=1 for input oriented and (-1)^O=-1 for output oriented

Scores=zeros(N,per,per);
%%%%
for s=1:per
    for t=1:per
        A=[Input(:,:,s),zeros(m_I,1);
            -Output(:,:,s),zeros(m_O,1)];
        %%%% Calculating the effeiciency score for each DMU
        for p=1:N
            A(:,N+1)=zeros(m_I+m_O,1);
            A(:,N+1)=[-delta*Input(:,p,t);(1-delta)*Output(:,p,t)];
            b(:,1)=[(1-delta)*Input(:,p,t);-delta*Output(:,p,t)];

            if(isempty(I_S)&&isempty(O_S))
                %%%% Solve the linear model
                options = optimoptions('linprog','Algorithm','dual-simplex');
                [X,~,exitflag] = linprog(objvec,A,b,Aeq,beq,lb,ub,options); % X(N+1)=theta and X(i)=lambda_i
                if(exitflag==-2)%%% change the algorithm
                    options = optimoptions('linprog','Algorithm','interior-point-legacy');
                    [X,~] = linprog(objvec,A,b,Aeq,beq,lb,ub,options);
                end
            else
                if(MI==1)
                    disp('MI-DEA is not defined for stochastic variables')
                    return
                end
                %%%% Solve the nonlinear model
                X=zeros(N+1,1);
                options=optimoptions('fmincon','Algorithm','sqp');
                X=fmincon(@(X)objfun(X,Ort_O),X,A,b,Aeq,beq,lb,ub,@(X)confun(delta,X,p,e,I_S,mean_I,covM_I,O_S,mean_O,covM_O),options);
            end
            %%%% efficiency scores
            Eff(p,:)=[p,(1-Ort_O)*X(N+1)+Ort_O*1/X(N+1),X(p)]; %% number of DMU,efficiency scoe, and it's weight!
        end
        Scores(:,s,t)=Eff(:,2);
        clear A b
    end
end
Scores=squeeze(Scores);

if(MI==1)
    catch_up=Scores(:,2,2)./Scores(:,1,1);
    catch_up=squeeze(catch_up);
    frontier_shift=sqrt(Scores(:,1,1)./Scores(:,2,1).*Scores(:,1,2)./Scores(:,2,2));
    frontier_shift=squeeze(frontier_shift);
    tfpch=catch_up.*frontier_shift;
    %tfpch=squeeze(tfpch);
end
end

function delta=kronecker_delta(i,j)
if(i==j)
    delta=1;
else
    delta=0;
end
end

function delta=dirac_delta(t)
if(t==0)
    delta=inf;
else
    delta=0;
end
end

function f =objfun(X,d)
n=length(X)-1;
f=(-1)^d*X(n+1,1);
end

function [c,ceq]=confun(delta,X,p,e,I_S,mean_I,covM_I,O_S,mean_O,covM_O)
% Nonlinear inequality constraints
n=length(X)-1;
c_I=zeros(size(I_S));
c_O=zeros(size(O_S));
if(~isempty(I_S))
    Lambda(:,1)=X(1:n,1);
    Lambda(p,1)=Lambda(p,1)-(delta*X(n+1,1)+(1-delta));
    for i=1:max(size(I_S))
        c_I(i)=mean_I(i,:)*Lambda-e*sqrt(transpose(Lambda)*covM_I(:,n*(i-1)+1:n*(i-1)+n)*Lambda);
    end
end
if(~isempty(O_S))
    Lambda(:,1)=-X(1:n,1);
    Lambda(p,1)=Lambda(p,1)+(delta+(1-delta)*X(n+1,1));
    for i=1:max(size(O_S))
        c_O(i)=mean_O(i,:)*Lambda-e*sqrt(transpose(Lambda)*covM_O(:,n*(i-1)+1:n*(i-1)+n)*Lambda);
    end
end
%Nonlinear equality constraints
c=[c_I,c_O];
ceq=[];
end