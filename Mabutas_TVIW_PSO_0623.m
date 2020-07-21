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
maxtrials=25;
maxIt=200000;
swarmSize=20;
while swarmSize<=100
	while maxIt<=400000
		while it<=maxtrials
			clc;
			clearvars -except it bestCosts maxIt maxtrials swarmSize;
			close all;

			%%Problem
			costFunction = @(x) sphere(x);                     % Cost Function
            prob_Dimensions=30-0;                              % Number of Unknown (Decision) Variables- based on number of dimensions
            varSize=[1 prob_Dimensions];                       % Matrix Size of Decision Variables
            prob_lower_bound=-100;                             % Lower Bound of Decision Variables
            prob_upper_bound=100;

            %%Parameters
            wmax=0.9;
            wmin=0.4;
            w=inf;
            c1=2.0;
            c2=2.0;

            %% Initialize
            empty_particle.position=[];             %Current Particle
            empty_particle.velocity=[];
            empty_particle.cost=[];             
            empty_particle.best.position=[];        %Personal Best
            empty_particle.best.permutation=[];
            empty_particle.best.cost=[];

            particle=repmat(empty_particle,swarmSize,1);
            new_particle.position=[];
            new_particle.cost=[];
            % Initialize global best
            globalBest.cost=inf;

            bestCosts=repmat(new_particle,maxIt,1); %holds all best costs costs for every iteration
            bc=zeros(maxIt,1);

            %%t=0
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

                %update Global Best
                if particle(i).best.cost<globalBest.cost
                    globalBest=particle(i).best;
                end

                %fprintf('[%f,%f,%f] \t %f\n',particle(i).position,particle(i).cost);
            end

            %% Main Loop
            t=1;
            while t<=maxIt
            	w=(wmax-wmin)*((maxIt-t)/maxIt)+wmin;
            	weights(t)=w;
            	for i=1:swarmSize
            		for d=1:prob_Dimensions
            			particle(i).velocity(1,d)=w*particle(i).velocity(1,d)+...
            				c1*rand().*(particle(i).best.position(1,d)-particle(i).position(1,d))+...
            				c2*rand().*(globalBest.position(1,d)-particle(i).position(1,d));

            			%Velocity clamping
            			if particle(i).velocity(1,d)>prob_upper_bound || particle(i).velocity(1,d)<prob_lower_bound
            				particle(i).velocity(1,d)=sign(particle(i).velocity(1,d)).*prob_upper_bound;
            			end
            			particle(i).position(1,d)=particle(i).position(1,d)+...
                            particle(i).velocity(1,d);
                        %Position Clamping
                        if particle(i).position(1,d)>prob_upper_bound ||particle(i).position(1,d)<prob_lower_bound
                        	particle(i).position(1,d)=sign(particle(i).position(1,d)).*prob_upper_bound;
                        end
            		end

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

            	t=t+1;
                
                if mod(t,1000)==0
                    num=it;
                    filename=sprintf('%s_%d_%d_%d.mat','TVIWPSO_f1_sphere',num,maxIt,swarmSize);
                    save(filename);
                end
            end
            
            figure;
            plot(bestCosts.cost,'LineWidth',2);
            semilogy(bc,'LineWidth',2);
            xlabel('iteration');
            ylabel('best cost');
            
            num=it;
            filename=sprintf('%s_%d_%d_%d.mat','TVIWPSO_f1_sphere',num,maxIt,swarmSize);
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
end

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