% 25 trials with parameters: 20 particles, 200,000 max evals
% 25 trials with parameters: 20 particles, 300,000 max evals
% 25 trials with parameters: 20 particles, 400,000 max evals
% 25 trials with parameters: 50 particles, 200,000 max evals
% 25 trials with parameters: 50 particles, 300,000 max evals
% 25 trials with parameters: 50 particles, 400,000 max evals
% 25 trials with parameters: 100 particles, 200,000 max evals
% 25 trials with parameters: 100 particles, 300,000 max evals
% 25 trials with parameters: 100 particles, 400,000 max evals

rng('shuffle');
it=1;
maxtrials=3;
maxIt=200000;
swarmSize=20;
while swarmSize<=100
    while maxIt<=400000
        while it<=maxtrials
            clc;
            clearvars -except it maxIt maxtrials swarmSize;
            close all;

            %% Problem Definition

            costFunction = @(x) schwefel222fcn(x);                     % Cost Function

            prob_Dimensions=30-0;                              % Number of Unknown (Decision) Variables- based on number of dimensions

            varSize=[1 prob_Dimensions];                       % Matrix Size of Decision Variables

            prob_lower_bound=-10;                             % Lower Bound of Decision Variables
            prob_upper_bound=10;                              % Upper Bound of Decision Variables


            %% Parameters
            
            w=inf;
            c=1.49618;
            CR=0.05;

            %% Initialize
            empty_particle.position=[];             %Current Particle
            empty_particle.velocity=[];
            empty_particle.cost=[];             
            empty_particle.best.position=[];        %Personal Best
            empty_particle.best.permutation=[];
            empty_particle.best.cost=[];

            particle=repmat(empty_particle,swarmSize,1);

            % Initialize global best
            globalBest.cost=inf;

            new_particle.position=[];
            new_particle.cost=[];
            G=zeros(swarmSize,1);
            V=repmat(new_particle,swarmSize,1);
            U=repmat(new_particle,swarmSize,1);
            PiNEW=repmat(new_particle,swarmSize,1);
            for i=1:swarmSize
                PiNEW(i).G=0;
                PiNEW(i).history=[];
            end

            bestCosts=repmat(new_particle,maxIt,1); %holds all best costs costs for every iteration
            bc=zeros(maxIt,1);
            %% t=0

            for i=1:swarmSize
                %random solution generation for t=0
                particle(i).position=deg2rad(unifrnd(prob_lower_bound, prob_upper_bound, varSize));

                %initialize velocity between [vmin=xmin vmax=xmax]
                particle(i).velocity=prob_lower_bound+rand(1,prob_Dimensions)*(prob_upper_bound-(prob_lower_bound));

                %evaluation
                particle(i).cost=costFunction(particle(i).position);

                %update personal best to current location and cost
                particle(i).best.position=particle(i).position;
                particle(i).best.cost=particle(i).cost;
                localBests(i)=particle(i).cost;     %Arumugam: used to store current localBests at each particle

                %update Global Best
                if particle(i).best.cost<globalBest.cost
                    globalBest=particle(i).best;
                end

                %fprintf('[%f,%f,%f] \t %f\n',particle(i).position,particle(i).cost);
            end

            %% Main Loop
            t=1;
            while t<=maxIt
                %fprintf('Iteration %d:\n',t);
                %% Exemplar Part
                r1=rand(1,prob_Dimensions);
                r2=rand();
                jrand=randi(varSize);
                for i=1:swarmSize
                    %artithmetical crossover
                    j=randperm(prob_Dimensions);
                    particle(i).best.permutation=particle(i).best.position(:,j);
                    V(i).position=r1.*particle(i).best.position+...
                        (1-r1).*particle(i).best.permutation;

                    %differential evolution crossover
                    for j=1: prob_Dimensions
                        if (r2 <= CR) || (j==jrand)
                            U(i).position(1,j)=V(i).position(1,j);
                        else
                            U(i).position(1,j)=particle(i).best.position(1,j);
                        end
                    end

                    %Evaluate U
                    U(i).cost=costFunction(U(i).position);

                    %competitive selection
                    if U(i).cost<=particle(i).best.cost
                        PiNEW(i).position=U(i).position;
                        PiNEW(i).cost=U(i).cost;
                    else
                        PiNEW(i).position=particle(i).best.position;
                        PiNEW(i).cost=particle(i).best.cost;
                    end
                %fprintf('P_%dNEW: [%f,%f,%f]\t %f\n',i,PiNEW(i).position,PiNEW(i).cost);    
                end

                %% Adjustment Part
                w=(1.1-(globalBest.cost/(mean(localBests(:))))); %Arumugam Method
                weights(t)=w;
                for i=1:swarmSize
                    for d=1:prob_Dimensions
                        particle(i).velocity(1,d)=w*particle(i).velocity(1,d)+...
                            c*rand()*(PiNEW(i).position(1,d)-particle(i).position(1,d));
                        particle(i).position(1,d)=particle(i).position(1,d)+...
                            particle(i).velocity(1,d);
                    end
                    %fprintf('[%f,%f,%f] \t %f\n',particle(i).position,particle(i).cost);

                    %evaluation of new position
                    particle(i).cost=costFunction(particle(i).position);

                    %update personal best
                    if particle(i).cost<particle(i).best.cost
                        particle(i).best.position=particle(i).position;
                        particle(i).best.cost=particle(i).cost;
                        localBests(i)=particle(i).cost;     %Arumugam: used to store current localBests at each particle

                        %update Global Best
                        if particle(i).best.cost<globalBest.cost
                            globalBest=particle(i).best;
                        end
                    end
                end

                %store best cost value
                bestCosts(t).position=globalBest.position;
                bestCosts(t).cost=globalBest.cost;
                bc(t)=globalBest.cost;
                %fprintf('Best Cost: [%f,%f,%f]\t %f\n',bestCosts(t).position,bestCosts(t).cost);

                %% check G=7 for PiNEW
                for i=1:swarmSize
                    if isequal(PiNEW(i).position,PiNEW(i).history)
                        PiNEW(i).G=PiNEW(i).G+1;
                        if PiNEW(i).G==7
                            %pause(5);
                            random_exemplars=randperm(swarmSize,floor(swarmSize*0.5));
                            best_rand_index=random_exemplars(1); %set first random exemplar as the lowest
                            for ii=2:swarmSize*0.5 %for every random exemplar from 2 to N*0.5
                                if PiNEW(random_exemplars(ii)).cost<PiNEW(best_rand_index).cost %check if current random exemplar is less than best rand exemplar
                                    best_rand_index=random_exemplars(ii); %change best rand exemplar if so
                                end
                            end
                            %check best fit, set to current PiNEW
                            PiNEW(i).position=PiNEW(best_rand_index).position;
                            PiNEW(i).cost=PiNEW(best_rand_index).cost;
                            PiNEW(i).G=0;
                        end
                    else
                        PiNEW(i).history=PiNEW(i).position;
                        PiNEW(i).G=0;
                    end
                end

            %% End Main Loop
                t=t+1;
                
                %saves current progress at every 1000th iteration
                if mod(t,30000)==0
                    num=it;
                    filename=sprintf('%s_%d_%d_%d.mat','MPSOCO_f2_schwefel222fcn_Arumugam',num,maxIt,swarmSize);
                    save(filename);
                end
            end

            %% Results

    %         figure;
    %         plot(bestCosts.cost,'LineWidth',2);
    %         semilogy(bc,'LineWidth',2);
    %         xlabel('iteration');
    %         ylabel('best cost');

            num=it;
            filename=sprintf('%s_%d_%d_%d.mat','MPSOCO_f2_schwefel222fcn_Arumugam',num,maxIt,swarmSize);
            save(filename);
            it=it+1;
        end
        maxIt=maxIt+100000;
        it=1;
    end
    if swarmSize==20
        swarmSize=swarmSize+30;
    elseif swarmSize==50
        swarmSize=swarmSize+50;
    else
        swarmSize=swarmSize+1;
    end
    maxIt=200000;
end

fprintf('***FINISHED RUNNING AT: %s***\n',datetime('now'));
%%Benchmark functions
function z= sphere(x)%[-100,100]
    z=sum(x.^2);
end

function scores = schwefel222fcn(x) %[-10,10]

    absx = abs(x);
    scores = sum(absx, 2) + prod(absx, 2);
end

function scores_ss = sumsquaresfcn(x) %[-10,10]
   
   [m, n] = size(x);
   x2 = x .^2;
   I = repmat(1:n, m, 1);
   scores_ss = sum( I .* x2, 2);
   
end

function scores_step= step(x) %[-100,100]
    scores_step=sum(floor(x+0.5).^2);
end

function scores_q = quarticfcn(x) %[-1.28,1.28]

    n = size(x, 2);
    
    scores_q = 0;
    for i = 1:n
        scores_q = scores_q + i *(x(:, i) .^ 4);
    end
     
    scores_q = scores_q + rand;
end

function scores = rosenbrockfcn(x) %[-10,10]
    scores = 0;
    n = size(x, 2);
    assert(n >= 1, 'Given input X cannot be empty');
    a = 1;
    b = 100;
    for i = 1 : (n-1)
        scores = scores + (b * ((x(:, i+1) - (x(:, i).^2)) .^ 2)) + ((a - x(:, i)) .^ 2);
    end
end