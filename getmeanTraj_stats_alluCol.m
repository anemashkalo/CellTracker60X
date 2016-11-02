%% scripts to plot the mean trajectories without binning using nuc:cyto ratio
% Get means over time but separately for different colony sizes
% data drawn from multiple .mat
% 
fr_stim = 16;          % 22(jan8data)  % 38 nov12data              %16 feb16          16 (july7 data)    (fr_stim = 13) july 26 data
delta_t = 17;         % 12 min        % 5 minutes                  15min             17min              (17min)
p = fr_stim*delta_t/60;
trajmin = 30;
tpt =100;
strdir = '*_out.mat';
ff = dir(strdir);
C1 = {'c','b','r','m','k'};
colormap = prism;
signalbin = [];
notbinned = [];
lastT = [];
clear traces
totalcol = zeros(10,1);
q = 1;
ucol = 3;% up to which size colonies to look at
N = 30;
   
for k=1:length(ff)
    outfile = ff(k).name;
    load(outfile,'colonies');
    if ~exist('colonies','var');
        disp('does not contain colonies structure')
    end
    numcol = size(colonies,2); % how many colonies were grouped within the frame
    traces = cell(1,numcol);
    for j = 1:numcol
        fr_stim_col = find(colonies(j).onframes == fr_stim);
        if ~isempty(fr_stim_col)           % new segmentation
            colSZ1 =colonies(j).ncells_actual(fr_stim_col) ;
            colSZ2 =colonies(j).ncells_actual(end-1) ;          % check what was the colony size at the last timepoint
            
            %colSZ =colonies(j).numOfCells(timecolSZ-1) ;                                          % for old mat files analysis
            nc = colSZ1 ;  %&& (colSZ2 == nc)
            traces{j} = colonies(j).NucSmadRatio;
            %traces{j} = colonies(j).NucSmadRatioOld;                                            % for old mat files analysis
            totalcol(nc) = totalcol(nc)+1; % determine the number of colonies (not traces)
            for h = 1:size(traces{j},2)
                [r,~] = find(isfinite(traces{j}(:,h)));                  %
                dat = zeros(tpt,1);
                dat(r,1) = traces{j}(r,h);
                %s = mean(nonzeros(dat(end-N:(end-1))));                        % take the mean of the last chunk of N time points
                if length(nonzeros(dat))>trajmin                                %FILTER OUT SHORT TRAJECTORIES
                    disp(['filter trajectories below' num2str(trajmin)]);
                    disp(['use' num2str(length(nonzeros(dat)))]);
                    figure(nc), plot(dat,'-*','color',colormap(nc,:));hold on          % here plot the traces for specific colony size,no binning
                    notbinned{nc}(:,q+size(dat,2)-1) = dat;
                    figure(nc), ylim([0 2]);
                    %errnotbinned{nc}(:,q+size(dat,2)-1) = std(nonzeros(dat(end-N:(end-1))));  % here store the sd for traces that meat condition
                    %lastT{nc}(:,q+size(dat,2)-1) = s;
                    % disp(q+sz-1)
                end
                
            end
            q = q+size(dat,2);
            
        end
    end
end                       % new segmentation
    
    
%     figure(nc)
%     ylim([0 2.5]);
%     xlim([0 (tpt+10)]);
%     ylabel('mean Nuc/Cyto smad4  ');
%     xlabel('frames');

%% average those trajectories and get the other stats
C = {'c','b','r'};
colormap = prism;
tpt = tpt-1;
vect = (1:tpt)';
ucolmean = zeros(tpt,ucol);
ucolerr= zeros(tpt,ucol);
startpoint = 4;
lastmean = zeros(tpt,ucol);

for j =1:3%size(notbinned,2)% loop over colony sizes;
    for k=1:(size(notbinned{j},1)-1)
        for ii = 1:(size(notbinned{j},2))
            if (notbinned{j}(k,ii) > 1.4) || (notbinned{j}(k,ii) < 0.45);                        % remove outliers
                
                notbinned{j}(k,ii) = 0;
            end
        end
    end
end


% average over time points, each colony separately
for j =1:3%size(notbinned,2)% loop over colony sizes;         
for k=1:(size(notbinned{j},1)-1)
    ucolmean(k,j) = mean(nonzeros(notbinned{j}(k,:)));   % mean over nonzero values of signaling at each time point
    ucolerr(k,j) = std(nonzeros(notbinned{j}(k,:)));     % error over nonzero values of signaling at each time point
end
end
for k=1:3%size(notbinned,2)
    figure(4), plot(ucolmean(startpoint:end,k),'-.','color',colormap(k,:),'Markersize',18,'Linewidth',3), hold on
    ylim([0.4 2]);title('Signaing','FontSize',21);xlabel('Frames')
    legend('1','2','3','Location','northwest');
    figure(6),plot(totalcol,'m-*','Markersize',16,'linewidth',3), hold on,
    ylim([0 130]);
    title('Colonies','FontSize',21);xlabel('Colony Size');
end
%%
for k=1:ucol
    figure(4), plot(ucolmean(startpoint:end,k),'-.','color',colormap(k,:),'Markersize',18,'Linewidth',3), hold on
    ylim([0.4 2]);title('Signaing','FontSize',21);xlabel('Frames')
    legend('1','2','3','Location','northwest');
    figure(5), plot(ucolerr(startpoint:end,k),'*','color',colormap(k,:),'Markersize',16), hold on
    ylim([0 0.5]);title('SD','FontSize',21);xlabel('Frames')
    legend('1','2','3');
    figure(6),plot(totalcol,'m-*','Markersize',16,'linewidth',3), hold on,
    ylim([0 130]);
    title('Colonies','FontSize',21);xlabel('Colony Size');
end
 % see if there is reason to bin
for k=1:ucol
   
    figure(7), subplot(1,ucol,k),histogram((lastT{k}(isfinite(lastT{k}))),10,'normalization','pdf','FaceColor',C{k}), legend([ num2str(k) '-cell colony' ]);xlim([0.5 2]);
    title(['Last ' num2str(N) ' tpts'],'FontSize',21); ylabel('Frequency');
       
end

for k=1:ucol
    x = size(nonzeros((errnotbinned{k}(isfinite(errnotbinned{k})))),1);
    y = nonzeros((errnotbinned{k}(isfinite(errnotbinned{k}))));
    figure(8), scatter((1:x),y(:,1),[],C{k}); hold on,%'*','color',C{k},'Markersize',16
    title(['SD during last ' num2str(N) ' tpts'],'FontSize',21); ylabel('Final signaling');
    xlim([0 50]);
    ylim([0 0.5]);
    legend('1','2','3');
    box on
end

    
    
    
    
    

 