
function [datcell] = AnalyzeCellTraces_AN(dir,col,N,fr_stim,delta_t,flag)

% plot cell traces for specific colonies using the output from the runTracker
% EDS and runcolonygrouping
% dir = direcory with the outfiles, merged for all time groups
% delta_t = interval between the frames , in minuted (converted to hours in
% the code)
% N = size of the colonies that you want to look at
% col = [2 3] ratio taken as 2/3: columns of the cell().fdata  -> 2 -
% mean intensity in the nucleus, 3 - in cyto (all taken from gfp channel)
% region
% flag = if 1, then additional stats will be calculated and plotted 
% for this peaks should have 9 columns: 9th has the size of the colony that
% the cell belongs to
% 
[nums, files]=folderFilesFromKeyword(dir,'Outfile_19totest');%['Outfile_000' num2str(pos) '_tps']

%found_cells = cell(length(peaks),1);
for j=1:length(nums)

matfile = files(j).name;

pp = load(matfile,'peaks','cells'); % AN
peaks = pp.peaks;
cells = pp.cells;
vect{j} = 1:length(peaks);

vect{j} = (vect{j}.*delta_t)./60;% x axis in units of hours  

trN = size(cells,2);
trN =(1:trN);

 %trN = max(peaks{1}(:,4)); % number of trajectories within the frame
%c = {'r','g','b','k','m','c','m','k'};
colors = jet(2*length(nums));
alldata=zeros(length(peaks),length(trN));%length(trN)
 
for xx = 1:length(trN)
    nc = size(cells(xx).data,1);
    for k=1:nc;
        if cells(xx).data(k,5)==N && ~(cells(xx).data(k,4)== -1);
            alldata(k,xx) = cells(xx).fdata(k,col(1))./cells(xx).fdata(k,col(2));
        end
        
    end
end

datcell{j} = alldata;
end

p = fr_stim*delta_t/60;
%bkgsign = zeros(length(datcell),3);
if flag == 1
for j=1:length(datcell)
    
    for k=1:size(datcell{j},2)
        colors2 = hot(size(datcell{j},2));
        if  length(nonzeros(datcell{j}(:,k)))>0%40 % this parameter should be eliminated
            figure(1),plot(vect{j},datcell{j}(:,k),'-*','color',colors2(k,:));
           % avgsign(j) = mean(datcell{j}(:,k));
            legend(['bmp4 added at ' num2str(p) 'hours']);
            hold on
                      
        end
        
    end
    
    mean_before{j} = mean(nonzeros(datcell{j}(1:fr_stim,:)));
    mean_after{j} =  mean(nonzeros(datcell{j}((fr_stim+1):end,:)));  
    %amplitude{j}(:,:) = abs(mean(nonzeros(datcell{j}((fr_stim-1),:)))- mean(nonzeros(datcell{j}((fr_stim+5),:))));%5 for 5min delta_t
    amplitude{j} = abs(mean_before{j}-mean_after{j});
end

    ylim([0 2.4])
    figure(1),title(['CellTraces for colonies of size ' num2str(N) ],'fontsize',20);
    ylabel('nuc/cyto raio');
    xlabel('Time, hours');
    
    mean_before = cell2mat(mean_before);
    mean_after = cell2mat(mean_after);
    amplitude = cell2mat(amplitude);
    figure(2),subplot(1,3,1),plot(mean_before,'*','color',colors(j,:),'markersize',20);
    hold on
    ylim([0 1.8]);
    xlim([0 length(datcell)]);
     legend('before');
     ylabel('nuc/cyto ratio, mean over time');
    xlabel('Positions');
   % title(['microCol of size ' num2str(N) ' ,all subplots'],'fontsize',15);
     figure(2),subplot(1,3,2),plot(mean_after,'*','color',colors(j,:),'markersize',20);
    hold on
    legend('after');
     ylim([0 1.8])
    xlim([0 length(datcell)]);
    ylabel('nuc/cyto ratio, mean over time');
    xlabel('Positions');
    title(['microCol of size ' num2str(N) ' ,all subplots'],'fontsize',15);
    figure(2),subplot(1,3,3),plot(amplitude,'*','color',colors(j,:),'markersize',20);
    hold on
    legend('amplitude');
     ylim([0 0.5])
    xlim([0 length(datcell)]);
    ylabel('nuc/cyto ratio change, mean over cells in frame','fontsize',9);
    xlabel('Positions');
  %  title(['microCol of size ' num2str(N) ],'fontsize',15);
end 
%save(['CellTraces_' num2str(N) ],'datcell','mean_before','mean_after','found_cells','amplitude');


end






