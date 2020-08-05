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
    collator=zeros(30,6);
    maxIt=i;
    popsize=20;
    while popsize<=100
        for k=1:10
            filename='MPSOCO_f1_sphere_Arumugam';
            file=[filename '_' num2str(k) '_' num2str(maxIt) '_' num2str(popsize) '.mat'];
            fprintf(['Processing ' file '\n']);
            load(file);
            
            clearvars -except i j k collator maxIt popsize bestCosts weights bc allCosts allWeights aveCosts aveWeights
            
            allCosts(:,k)=bc;
            allWeights(:,k)=weights';
            collatind=k;
            if popsize==50
                collatind=collatind+10;
            elseif popsize==100
                collatind=collatind+20;
            end
            
            collator(collatind,1)=popsize;
            collator(collatind,2)=k;
            [temp1,temp2]=min(bc(:));
            collator(collatind,3)=temp1;
            collator(collatind,4)=temp2;
            collator(collatind,5)=temp2;
            collator(collatind,6)=maxIt;
        end
            aveCosts=mean(allCosts,2);
            aveWeights=mean(allWeights,2);
            %Graph and hold on
            if popsize==20
                graph20=semilogy(aveCosts,'-r');
                hold on;
            elseif popsize==50
                graph50=semilogy(aveCosts,'-b');
                hold on;
            elseif popsize==100
                graph100=semilogy(aveCosts,'-m');
                hold off;
            end
            xlim([0 maxIt]);
            xlabel('iteration number');
            ylabel('cost value');
            title(['Arumugam Method Sphere Function Evaluation for ' num2str(maxIt) ' Iterations']);
            legend('20 Particles','50 Particles','100 Particles');
            set(gcf,'position',[350 100 700 500]);
            
            %break;
        if popsize==20
            popsize=popsize+30;
        elseif popsize==50
            popsize=popsize+50;
        elseif popsize==100
            popsize=popsize+1;
        end
    end
    %clear;
    break;
end


