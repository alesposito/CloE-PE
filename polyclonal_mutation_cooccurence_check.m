%% V7

% Initialize
clear all
close all
rng('shuffle')

% Simulation parameters
r0 = 1e-6; % mutational rate (1e-6)
dt = 1;    % time base
N  = 1000; % N^2 is total number of cells cells

tn = 1000000; % number of time points (100000)
rn = 2000;     % number of times the simulation is repeated (200)

kernel = [0 1 0; 1 0 1; 0 1 0]; % Omega = 4 kernel
%kernel = [1 1 1; 1 0 1; 1 1 1]; % Omega = 8 kernel
%kernel = [0 0 1 0 0; 0 1 1 1 0; 1 1 0 1 1; 0 1 1 1 0; 0 0 1 0 0]; % Omega = 12 kernel
%kernel = [1 1 1 1 1; 1 1 1 1 1; 1 1 0 1 1; 1 1 1 1 1; 1 1 1 1 1]; % Omega = 24 kernel

p0 = r0*dt;  
omega = sum(kernel(:))
%% Initialize GUI 
hf=figure
tic
if rn==1 %if single shot sim rn==1

    nx = 2; ny = 4;
    img_idx = [1 2 3 5 6 7];
    sTitles = {'SA','SB','2muts-1cell','n muts','SC','SD','1mut-2cells','nab/ncd'};
    
    for fi=1:nx*ny
        subplot(nx,ny,fi)
        if ismember(fi,img_idx)
            h{fi} = imagesc(zeros(N));
            axis off
            set(gca,'clim',[0 1])
        elseif fi==4        
            hp1=plot(0,0);
            hold on
            hp2=plot(0,0);
            hold off
            legend('ab','cd')
        elseif fi==8
            hp3=plot(0,0);
       
        end
        axis square
        title(sTitles{fi})

    end
else
    
    subplot(1,2,1)
    plot(0,0)
    ylabel('frequency')
    xlabel('time')
    axis square
    
    subplot(1,2,2)
    plot(0,0)
    ylabel('frequency')
    xlabel('number of ab muts at cd occurence')
    axis square
    
end

%%
% store, time of first co-mutations (ab, or cd) and the number of ab
% mutations when the first cd-mutation appears
memo_tab = [];
memo_tcb = [];
memo_nab = [];
%%
for ri=1:rn
    
    % initialize
    memo_ab =[];
    memo_cd =[];
    [SA SB SC SD] = deal(zeros(N));

    % cycle time
    for ti=1:tn
        
        % generate mutations with probability p0
        SA = (SA+(rand(N)<=p0))>=1;
        SB = (SB+(rand(N)<=p0))>=1;
        SC = (SC+(rand(N)<=p0))>=1;
        SD = (SD+(rand(N)<=p0))>=1;

        % identify cells with both mutations (reoccurence of mutations is
        % neglected, i.e. once the mutation occurred, the fact that within
        % the simulation might reoccurr is excluded by design
        SAB = SA & SB;

%       check if mutation co-occured in the first neighborhood (8-cells)
        SCD = SD&(ordfilt2(SC,omega,kernel));        

        memo_ab(ti)= nnz(SAB)/N^2;
        memo_cd(ti)= nnz(SCD)/N^2;

        if mod(tn,100)==0 & rn==1
            set(h{1},'cdata',SA);
               set(h{2},'cdata',SB);
               set(h{3},'cdata',SC);
               set(h{5},'cdata',SD);
               set(h{6},'cdata',SAB);
               set(h{7},'cdata',SCD);
             
               hp1.XData=(1:ti);
               hp1.YData=memo_ab;
               hp2.XData=(1:ti);
               hp2.YData=memo_cd;
               hp3.XData=(1:ti);
               hp3.YData=memo_ab./memo_cd;
              
               drawnow
        end   

        if rn>1

            if memo_ab(ti)*memo_cd(ti)>0
                break
            end
        end

    end
    
    memo_tab(ri) = min(find(memo_ab.*N^2>=1));
    memo_tcd(ri) = min(find(memo_cd.*N^2>=1));
    memo_ncd(ri)=N^2 * memo_cd(ti);

    if mod(ri,10)==0
    
        xt = (0:.25:8);
        xs = (0:10:100);
        h1 = hist(memo_tab/365,xt);
        h2 = hist(memo_tcd/365,xt);
        h3 = hist(memo_ncd,xs);

        subplot(1,2,1)
        plot(xt,[h1/sum(h1); h2/sum(h2)]')
        axis square

        subplot(1,2,2)
        boxplot(memo_ncd,'plotstyle','compact')
        axis square
        set(gca,'ylim',[0 100])

        omega*median(memo_tcd)/365
        omega*mad(memo_tcd,1)/365

        omega*median(memo_tab)/365
        omega*mad(memo_tab,1)/365

        median(memo_ncd)
        mad(memo_ncd,1)

        nnz(memo_ncd==1)/rn
        
        drawnow
        ri
         saveas(hf,['figdump_' num2str(omega) '_' num2str(ri) '.fig'])
         save(['wsdump_' num2str(omega) '_' num2str(ri) '.mat'])
    end
end

beep

toc


