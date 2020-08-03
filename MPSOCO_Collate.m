%Goal: Collate costs and weight
%Include Semilog Graphing all setups while iterating (set gca color etc)
clear;
clc;
close all;
for i=10000:10000:100000
    allCosts=zeros(i,10);
    allWeights=zeros(i,10);
    aveCosts=zeros(i,1);
    aveWeights=zeros(i,1);
    maxIt=i;
    popsize=20;
    while popsize<=100
        for k=1:10
            filename='MPSOCO_f1_sphere_Arumugam';
            file=[filename '_' num2str(k) '_' num2str(maxIt) '_' num2str(popsize) '.mat'];
            fprintf(['Processing ' file '\n']);
            load(file);
            
            clearvars -except i j k maxIt popsize bestCosts weights bc allCosts allWeights aveCosts aveWeights
            
            allCosts(:,k)=bc;
            allWeights(:,k)=weights';
        end
            aveCosts=mean(allCosts,2);
            aveWeights=mean(allWeights,2);
            %Graph and hold on
            
            %break;
        if popsize==20
            popsize=popsize+30;
        elseif popsize==50
            popsize=popsize+50;
        elseif popsize==100
            popsize=popsize+1;
        end
    end
    clear;
    %break;
end


% for index=1:25
%     clearvars -except index allCosts allcosts;
%     filename='PSOCO_test_rosenbrockfcn';
%     file=[filename '____' num2str(index) '.mat'];
%     load(file);
%     for j=1:300000
%         %allCosts(j,index).cost=bestCosts(j).cost;
%         allcosts(j,index)=bestCosts(j).cost;
%         %allcosts(j,index)=bestCosts(j).cost;
%         %allCosts(j,index).position=bestCosts(j).position;
%     end
%     %allPos(:,i)=bestCosts(:).position;
%     
%     %mean_ac=(mean2(allCosts));
%     %std_ac=(std2(allCosts));
%     %min_ac=min(allCosts);
%     
%     %fprintf('Mean: %.10f\nSD: %.10f\n Best: %.10f\n',mean_ac,std_ac,min_ac);
% end
% %semilogy(allcosts);
% filename=sprintf('%s.mat','PSOCO_F6');
% save(filename);

