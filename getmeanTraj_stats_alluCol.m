%% scripts to plot the mean trajectories without binning using nuc:cyto ratio
% Get means over time but separately for different colony sizes
% new figure for each new colony size ( data drawn from multiple .mat
% files)
fr_stim = 5;        %22(jan8data)  %38 nov12data              %16 feb16     16 (july7 data)   %(fr_stim = 13) july 26 data
delta_t = 17;%15       % 12 min       % 5 minutes               15min             17min              (17min)
p = fr_stim*delta_t/60;
trajmin = 30;%30
tpt =99;%100 81 99 83
strdir = '*_outPluriset.mat';%outFebsetBGan
ff = dir(strdir);%'*_60X_testparam_allT.mat'ff = dir('*60Xjan8_R*.mat');
C1 = {'c','b','r'};% for colorcoding colony size ( if look at means of all nc-colonies
signalbin = [];
notbinned = [];
lastT = [];
clear traces
totalcol = zeros(6,1);
q = 1;
ucol = 3;% up to which size colonies to look at
N = 40;
for nc = 1:ucol
    
    for k=1:length(ff)
        outfile = ff(k).name; %nms{k};
        load(outfile,'colonies');
        if ~exist('colonies','var');
            disp('does not contain colonies structure')
        end
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        traces = cell(1,numcol);
        nucgfp= cell(1,numcol);
        for j = 1:numcol
            if size(colonies(j).ncells_actual,1)>fr_stim  %&&             % new segmentation
                colSZ1 =colonies(j).ncells_actual(fr_stim) ;
                colSZ2 =colonies(j).ncells_actual(end-1) ;   % check what was the colony size at the last timepoint
                
                % new segmentation
                %colSZ =colonies(j).numOfCells(timecolSZ-1) ;                                                      % for old mat files analysis
                if colSZ1 == nc %&& (colSZ2 == nc)
                    
                    traces{j} = colonies(j).NucSmadRatio;                                                 % colonies(j).NucSmadRatio(:)
                    %traces{j} = colonies(j).NucOnlyData;   % to look at only
                    %nuclear GFP times the cell nucArea
                    %traces{j} = colonies(j).NucSmadRatioOld;                                                     % for old mat files analysis
                    totalcol(nc) = totalcol(nc)+1; % determine the number of colonies (not traces)
                    for h = 1:size(traces{j},2)
                        [r,~] = find(isfinite(traces{j}(:,h)));                  %
                        dat = zeros(tpt,1);
                        dat(r,1) = traces{j}(r,h);
                        s = mean(nonzeros(dat(end-N:(end-1))));                        % take the mean of the last chunk of N time points
                        if length(nonzeros(dat))>trajmin                               %FILTER OUT SHORT TRAJECTORIES
                            disp(['filter trajectories below' num2str(trajmin)]);
                            disp(['use' num2str(length(nonzeros(dat)))]);
                            figure(nc), plot(dat,'-*','color',C1{nc});hold on          % here plot the traces for specific colony size,no binning
                            notbinned{nc}(:,q+size(dat,2)-1) = dat;
                            errnotbinned{nc}(:,q+size(dat,2)-1) = std(nonzeros(dat(end-N:(end-1))));  % here store the sd for traces that meat condition
                            lastT{nc}(:,q+size(dat,2)-1) = s;
                            % disp(q+sz-1)
                        end
               
                    end
                    q = q+size(dat,2);
         
                end
            end
        end                       % new segmentation
    end
    
    figure(nc)
    ylim([0 2.5]);
    xlim([0 (tpt+10)]);
    ylabel('mean Nuc/Cyto smad4  ');
    xlabel('frames');
end
%end
% figure, plot(1:size(totalcol,1),totalcol,'r-*','markersize',18,'linewidth',3);
% xlabel('cells per colony','fontsize',20);
% ylabel('totla colonies','fontsize',20);
% title('colony size distribution','fontsize',20)
%% average those trajectories and get the other stats
C = {'c','b','r'};
tpt = tpt-1;
vect = (1:tpt)';
ucolmean = zeros(tpt,ucol);
ucolerr= zeros(tpt,ucol);
startpoint = 5;
lastmean = zeros(tpt,ucol);

for j =1:size(notbinned,2)% loop over colony sizes;
    for k=1:(size(notbinned{j},1)-1)
        for ii = 1:(size(notbinned{j},2))
            if (notbinned{j}(k,ii) > 2.1) || (notbinned{j}(k,ii) < 0.5);                        % remove outliers
                
                notbinned{j}(k,ii) = 0;
            end
        end
    end
end


% average over time points, each colony separately
for j =1:size(notbinned,2)% loop over colony sizes;         
for k=1:(size(notbinned{j},1)-1)
    ucolmean(k,j) = mean(nonzeros(notbinned{j}(k,:)));   % mean over nonzero values of signaling at each time point
    ucolerr(k,j) = std(nonzeros(notbinned{j}(k,:)));     % error over nonzero values of signaling at each time point
end
end
for k=1:ucol
    figure(4), plot(ucolmean(startpoint:end,k),'-*','color',C{k},'Markersize',16), hold on
    ylim([0.4 2]);title('Signaing','FontSize',21);xlabel('Frames')
    legend('1','2','3','Location','northwest');
    figure(5), plot(ucolerr(startpoint:end,k),'*','color',C{k},'Markersize',16), hold on
    ylim([0 0.5]);title('SD','FontSize',21);xlabel('Frames')
    legend('1','2','3');
    figure(6), plot(totalcol,'m-*','Markersize',16,'linewidth',3), hold on,
    ylim([0 130]);
    title('Colonies','FontSize',21);xlabel('Colony Size');
end
 % see if there is reason to bin
for k=1:ucol
   
    figure(7+k), histogram((lastT{k}(isfinite(lastT{k}))),'normalization','pdf'), legend([ num2str(k) '-cell colony' ]);xlim([0 2]);
    title('Last tpt distribution','FontSize',21);xlabel('Final signaling');
    
    
end

for k=1:ucol
    
    figure(10), plot(nonzeros((errnotbinned{k}(isfinite(errnotbinned{k})))),'*','color',C{k},'Markersize',16),
    title('SD during last time period','FontSize',21); ylabel('Final signaling');
    xlim([0 2]);
    ylim([0 2]);
    hold on
    legend('1','2','3');
end

    
    
    
    
    

 