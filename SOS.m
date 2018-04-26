% ----------------------------------------------------------------------- %
% Symbiotic Organisms Search(SOS) for unconstrained benchmark problems
% a simplified version, last revised: 2014.08.27
% ----------------------------------------------------------------------- %
% Files of the Matlab code used in the artile:
%
% Min-Yuan Cheng, Doddy Prayogo,
% Symbiotic Organisms Search: A new metaheuristic optimization algorithm,
% Computers & Structures 139 (2014), 98-112
% http://dx.doi.org/10.1016/j.compstruc.2014.03.007
%
% ----------------------------------------------------------------------- %
% Written by Doddy Prayogo at National Taiwan University of Science and
% Technology (NTUST)
% Email: doddyprayogo@ymail.com
% ----------------------------------------------------------------------- %

%% --- MAIN OPTIMIZER ---
function [bestOrganism bestFitness pp]=SOS(ecosize,funnum,xx,yy)
tic;
% Outputs: best organism/solution and best fitness
% Inputs: ecosystem/population size and # of benchmark problems
% Example: [A,B]=SOS (50,17), SOS will solve Sphere (F17) with 50 organisms
% (please see the "OBJECTIVE FUNCTIONS" and "SETUP" sub-functions)
%format compact
fprintf('-------------------------------------------------------------------------\n');
fprintf('  Symbiotic Organisms Search(SOS) for unconstrained benchmark problems\n');
fprintf('-------------------------------------------------------------------------\n\n');

% --- Counters, Parameters & Matrix Initialization
[globalMin lb ub n maxFE]=terminate(funnum);
fprintf(' Ecosystem Size: %d\t\tMaxFE: %d\t\tFunctionNumber: %d',ecosize,maxFE,funnum);
fprintf('\n\n');
fprintf('-------------------------------------------------------------------------\n\n');

FE=0;                           % Function of Evaluation Counter
eco=zeros(ecosize,n);           % Ecosystem Matrix
fitness =zeros(ecosize,1);      % Fitness Matrix

% --- Ecosystem Initialization
for i=1:ecosize
    % Initialize the organisms randomly in the ecosystem
    eco(i,:)=rand(1,n).*(ub-lb)+lb;
    % Evaluate the fitness of the new solution
    fitness(i,:)=fobj(eco(i,:),funnum,xx,yy);
    % Increase the number of function evaluation counter
    %FE=FE+1;
end

% --- Main Looping
while FE<maxFE
    
    for i=1:ecosize % Organisms' Looping
        
        % Update the best Organism
        [bestFitness,idx]=min(fitness); bestOrganism=eco(idx,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mutualism Phase
        % Choose organism j randomly other than organism i
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        
        % Determine Mutual Vector & Beneficial Factor
        mutualVector=mean([eco(i,:);eco(j,:)]);
        BF1=round(1+rand); BF2=round(1+rand);
        
        % Calculate new solution after Mutualism Phase
        ecoNew1=eco(i,:)+rand(1,n).*(bestOrganism-BF1.*mutualVector);
        ecoNew2=eco(j,:)+rand(1,n).*(bestOrganism-BF2.*mutualVector);
        ecoNew1=bound(ecoNew1,ub,lb);
        ecoNew2=bound(ecoNew2,ub,lb);
        
        % Evaluate the fitness of the new solution
        fitnessNew1=fobj(ecoNew1,funnum,xx,yy);
        fitnessNew2=fobj(ecoNew2,funnum,xx,yy);
        
        % Accept the new solution if the fitness is better
        if fitnessNew1<fitness(i)
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
        end
        if fitnessNew2<fitness(j)
            fitness(j)=fitnessNew2;
            eco(j,:)=ecoNew2;
        end
        
        % Increase the number of function evaluation counter
        FE=FE+2;
        
        % End of Mutualism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Commensialism Phase
        
        % Choose organism j randomly other than organism i
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        
        % Calculate new solution after Commensalism Phase
        ecoNew1=eco(i,:)+(rand(1,n)*2-1).*(bestOrganism-eco(j,:));
        ecoNew1=bound(ecoNew1,ub,lb);
        
        % Evaluate the fitness of the new solution
        fitnessNew1=fobj(ecoNew1,funnum,xx,yy);
        
        % Accept the new solution if the fitness is better
        if fitnessNew1<fitness(i)
            fitness(i)=fitnessNew1;
            eco(i,:)=ecoNew1;
        end
        
        % Increase the number of function evaluation counter
        FE=FE+1;
        
        % End of Commensalism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Parasitism Phase
        
        % Choose organism j randomly other than organism i
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        
        % Determine Parasite Vector & Calculate the fitness
        parasiteVector=eco(i,:);
        seed=randperm(n);
        pick=seed(1:ceil(rand*n));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
        fitnessParasite=fobj(parasiteVector,funnum,xx,yy);
        
        % Kill organism j and replace it with the parasite
        % if the fitness is lower than the parasite
        if fitnessParasite < fitness(j)
            fitness(j)=fitnessParasite;
            eco(j,:)=parasiteVector;
        end
        
        % Increase the number of function evaluation counter
        FE=FE+1;
        
        % End of Parasitism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end % End of Organisms' Looping
    
    % Checking the termination criteria
    if bestFitness<globalMin
        break
    end
    
end % End of Main Looping

% --- Update the best Organism
[bestFitness,idx]=min(fitness); bestOrganism=eco(idx,:);

% --- Display the result
%disp(['Funnum: ', num2str(funnum)])
disp(['Best Fitness: ', num2str(bestFitness)])
disp(['Best Organism: ', num2str(bestOrganism)])
toc;


%% --- Boundary Handling ---
function a=bound(a,ub,lb)
a(a>ub)=ub(a>ub); a(a<lb)=lb(a<lb);

%% --- Objective Functions ---
function y=fobj(x,funnum,xx,yy)
% Benchmark function
%F28 = PolinomFitting [];

if funnum==28
    error=0;
    P=xx;
    F=yy; %gerçek y de?eri
    for j=1:size(xx,1)
        S=...
            x(1)*P(j)^20+...
            x(2)*P(j)^19+...
            x(3)*P(j)^18+...
            x(4)*P(j)^17+...
            x(5)*P(j)^16+...
            x(6)*P(j)^15+...
            x(7)*P(j)^14+...
            x(8)*P(j)^13+...
            x(9)*P(j)^12+...
            x(10)*P(j)^11+...
            x(11)*P(j)^10+...
            x(12)*P(j)^9+...
            x(13)*P(j)^8+...
            x(14)*P(j)^7+...
            x(15)*P(j)^6+...
            x(16)*P(j)^5+...
            x(17)*P(j)^4+...
            x(18)*P(j)^3+...
            x(19)*P(j)^2+...
            x(20)*P(j)^1+...
            x(21);%hesaplanan y de?eri
        error=error+abs(F(j)-S); %gerçek ve hesaplanan aras?ndaki fark?n mutlak de?eri
    end
    
    y = polyval(x,xx);    
    subplot(2,1,1)
    plot(xx,yy);
    subplot(2,1,2)
    plot(xx,y);
    pause(0.1)
    
    y=error;
end

%% --- SETUP (boundary limit, variable numbers, stopping criterion) ---
function [globalMin Lb Ub nd maxFE]=terminate(funnum)
maxFE=5e+8; % maximum number of function evaluation
Tol = 1e-12;

if funnum==28
    Lb=-1e+30; %alt sýnýr
    Ub=1e+30; %üst sýnýr
    nd=21; %optimize edilecek deðiþken sayýsý
    Lb = ones(1,nd)*Lb;
    Ub = ones(1,nd)*Ub;
    globalMin=0;
    globalMin = globalMin+Tol;
end
