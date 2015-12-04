
function [datcell,bkgsign] = AnalyzeCellTraces_AN(dir,col,N,fr_stim,delta_t,flag)

% plot cell traces for specific colonies using the output from the runTracker
% EDS and runcolonygrouping
% dir = direcory with the outfiles, merged for all time groups
% delta_t = interval between the frames , in minuted (converted to hours in
% the code)
% N = size of the colonies that you want to look at
% col = [6 7] ratio taken as 6/7: columns of the peaks with data -> 6 -
% mean intensity in the nucleus, 7 - in cyto (all taken from gfp channel)
% region
% flag = if 1, then additional stats will be calculated and plotted (NEED
% TO WORK ON, NOT FINAL)
% for this peaks should have 9 columns: 9th has the size of the colony that
% the cell belongs to
% 

[nums, files]=folderFilesFromKeyword(dir,'Outfile_00');%['Outfile_000' num2str(pos) '_t']


for j=1:length(nums)

matfile = files(j).name;

pp = load(matfile,'peaks','cells'); % AN
peaks = pp.peaks;
vect{j} = 1:length(peaks);
vect{j} = (vect{j}.*delta_t)./60;% x axis in units of hours  

for k=1:length(peaks)
if ~isempty(peaks{k})
trN(k) = max(peaks{k}(:,8));
end
end
trN = max(trN);
trN =(1:trN);

 %trN = max(peaks{1}(:,4)); % number of trajectories within the frame
c = {'r','g','b','k','m','c','m','k'};
alldata=zeros(length(peaks),length(trN));%length(trN)
 
for xx = 1:length(trN)
    
    for k=1:length(peaks)
        if ~isempty(peaks{k})
            
            nc = size(peaks{k},1);% number of cells found within each frame k
            for i=1:nc
                
                if peaks{k}(i,9) == N && peaks{k}(i,8) == trN(xx);%
                    alldata(k,xx) = peaks{k}(i,col(1))./peaks{k}(i,col(2));
                    
                end
                
            end
            
        end
    end
end


datcell{j} = alldata;
end
p = fr_stim*delta_t/60;
%bkgsign = zeros(length(datcell),3);
for j=1:length(datcell)
    for k=1:size(datcell{j},2)
        if ~isempty(datcell{j}(:,k))%&& size(datcell{j}(:,k),1)==length(vect{j})
            figure(1),plot(vect{j},datcell{j}(:,k),'*','color',c{j});
            legend(['bmp4 added at ' num2str(p) 'hours']);
            hold on
            
            bkgsign(j,1) = mean(datcell{j}(1:fr_stim,k));hold on% mean signaling before stimulation
            bkgsign(j,2) = mean(datcell{j}(fr_stim:end,k));hold on% mean signaling after stimulation
            %bkgsign(j,3) = abs(datcell{j}(fr_stim-1,k)-mean(datcell{j}(fr_stim:(fr_stim+4),k)));% value of the jump upon stimulation
            bkgsign(j,3) = abs(bkgsign(j,1)-bkgsign(j,2));
            
        end
    end
end
    ylim([0 2.4])
    figure(1),title(['CellTraces for colonies of size ' num2str(N) ],'fontsize',20);
    ylabel('nuc/cyto raio');
    xlabel('Time, hours');
    ylim([0 2.4])
    if flag == 1
    figure(2),subplot(1,2,1),plot(bkgsign(:,1),'r-*');
    legend('Before Stimulation');
    title(['MicroCol Size ' num2str(N) ]);
    ylabel('Mean nuc/cyto raio');
    xlabel('Frames');
    ylim([0 2.4])
    figure(2),subplot(1,2,2),plot(bkgsign(:,2),'b-*');
    legend('After Stimulation');
    title(['MicroCol Size ' num2str(N) ]);
    ylabel('Mean nuc/cyto raio');
    xlabel('Frames');
    ylim([0 2.4])
    figure(2),subplot(1,3,3),plot(bkgsign(:,3),'g-*');
    legend('Upon Stimulation');
    title(['MicroCol Size ' num2str(N) ]);
    ylabel('Mean nuc/cyto raio');
    xlabel('Frames');
    ylim([0 2.4])
    end

    
save(['CellTraces_' num2str(N) ],'datcell','vect','bkgsign');


end





