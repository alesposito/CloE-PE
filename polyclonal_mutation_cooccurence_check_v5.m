% NOTE: this code was not optimized for speed
% Initialize
clear all
close all
rng('shuffle')

% Simulation parameters
r0 = 1e-6; % mutational rate
dt = 1;    % time base
N  = 1000; % N^2 is total number of cells cells

tn = 100000; % number of time points
rn = 2000;    % number of times the simulation is repeated

p0 = r0*dt;  

%% Initialize GUI 
hf=figure

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
        else        
            plot(0,0)        
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

%        check if mutation co-occured in the first neighborhood (8-cells)
        SCD = ...
            (SD & imfilter(SC,[0 0 0; 1 0 0; 0 0 0])) |...
            (SD & imfilter(SC,[1 0 0; 0 0 0; 0 0 0])) |...
            (SD & imfilter(SC,[0 1 0; 0 0 0; 0 0 0])) |...
            (SD & imfilter(SC,[0 0 1; 0 0 0; 0 0 0])) |...
            (SD & imfilter(SC,[0 0 0; 0 0 1; 0 0 0])) |...
            (SD & imfilter(SC,[0 0 0; 0 0 0; 0 0 1])) |...
            (SD & imfilter(SC,[0 0 0; 0 0 0; 0 1 0])) |...
            (SD & imfilter(SC,[0 0 0; 0 0 0; 1 0 0]));

%        check if mutation co-occured in the first neighborhood (12-cells)
%         SCD = ...
%             (SD & imfilter(SC,[0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0])) |...
%             (SD & imfilter(SC,[0 0 0 0 0; 0 0 0 0 0; 1 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0])) |...
%             (SD & imfilter(SC,[0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 1; 0 0 0 0 0; 0 0 0 0 0])) |...
%             (SD & imfilter(SC,[0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0])) |...                        
%             (SD & imfilter(SC,[0 0 0; 1 0 0; 0 0 0])) |...
%             (SD & imfilter(SC,[1 0 0; 0 0 0; 0 0 0])) |...
%             (SD & imfilter(SC,[0 1 0; 0 0 0; 0 0 0])) |...
%             (SD & imfilter(SC,[0 0 1; 0 0 0; 0 0 0])) |...
%             (SD & imfilter(SC,[0 0 0; 0 0 1; 0 0 0])) |...
%             (SD & imfilter(SC,[0 0 0; 0 0 0; 0 0 1])) |...
%             (SD & imfilter(SC,[0 0 0; 0 0 0; 0 1 0])) |...
%             (SD & imfilter(SC,[0 0 0; 0 0 0; 1 0 0]));
% 
%         
        % check if mutation co-occured in the first neighborhood (4-cells)
%         SCD = ...
%             (SD & imfilter(SC,[0 0 0; 1 0 0; 0 0 0])) |...            
%             (SD & imfilter(SC,[0 1 0; 0 0 0; 0 0 0])) |...            
%             (SD & imfilter(SC,[0 0 0; 0 0 1; 0 0 0])) |...            
%             (SD & imfilter(SC,[0 0 0; 0 0 0; 0 1 0]));
            

        
        memo_ab(ti)= nnz(SAB)/N^2;
        memo_cd(ti)= nnz(SCD)/N^2;

        if mod(tn,100)==0 & rn==1
            set(h{1},'cdata',SA);
               set(h{2},'cdata',SB);
               set(h{3},'cdata',SC);
               set(h{5},'cdata',SD);
               set(h{6},'cdata',SAB);
               set(h{7},'cdata',SCD);

               subplot(2,4,4)
               hold off
               plot(memo_ab)
               hold on
               plot(memo_cd)

               subplot(2,4,8)
               hold off
               plot(memo_cd./memo_ab)

               drawnow
        end   

        if rn>1

            if memo_ab(ti)*memo_cd(ti)>0
                break
            end
        end

    end
    
    
    if mod(ri,10)==0
        memo_tab(ri) = min(find(memo_ab.*N^2>=1));
        memo_tcd(ri) = min(find(memo_cd.*N^2>=1));
        memo_ncd(ri)=N^2 * memo_cd(ti);
    
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

        12*median(memo_tcd)/365
        12*mad(memo_tcd,1)/365

        12*median(memo_tab)/365
        12*mad(memo_tab,1)/365

        median(memo_ncd)
        mad(memo_ncd,1)

        nnz(memo_ncd==1)/rn
        
        drawnow
        ri
        saveas(hf,['figdump' num2str(ri) '.fig'])
        save(['wsdump' num2str(ri) '.mat'])
    end
end

beep




