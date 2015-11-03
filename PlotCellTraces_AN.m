
function [alldata] = PlotCellTraces_AN(matfile, col)
% plot live cell analysis results
pp=load(matfile,'peaks','NucMasks','colonies'); % AN
peaks=pp.peaks;

colonies = pp.colonies;

nc = size(colonies{3},2);% number of colonies found in the image
sz = size(colonies{3},1);% size of the colonies found in the image
%colors = colorcube(ncell);
alldata = zeros(length(colonies),nc);
alldata2 = zeros(length(colonies),nc);
colors = {'r','g','b','m'};
nlines=zeros(length(colonies),nc);
    for k=1:nc


for ii=1:length(colonies)
    if ~isempty(colonies{ii})
                nlines(ii,k)=size(colonies{ii}(k),1);
    end
end
    end


alllines = sum(nlines);
alldata=zeros(alllines(1),nc);
alldata2=zeros(alllines(1),nc);
for k=1:nc


q=1;

for ii=1:length(colonies)
    if ~isempty(colonies{ii})
        
        if length(col)==1
            alldata(q:(q+nlines(ii)-1),k)=colonies{ii}(k).data(:,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
        else
            alldata(q:(q+nlines(ii)-1),k)=mean(colonies{ii}(k).data(:,col(1)))./mean(colonies{ii}(k).data(:,col(2)));
%             if colonies{ii}(k).ncells == 2
%             alldata2(q:(q+nlines(ii)-1),k)=mean(colonies{ii}(k).data(:,col(1)))./mean(colonies{ii}(k).data(:,col(2)));
%             end
        end
        q=q+nlines(ii);
    end
    
end

% figure(4), hold on
% plot(alldata(:,k),'-*','color',colors{k});
% ylim([0 1.5]);
end
end

% nlines=zeros(length(peaks),1);
% for ii=1:length(peaks)
%     nlines(ii)=size(peaks{ii},1);
% end
% alllines = sum(nlines);
% alldata=zeros(alllines,1);
% for k=1:ncell
% q=1;
% 
% for ii=1:length(peaks)
%     if ~isempty(peaks{ii})
%         
%         if length(col)==1
%             alldata(q:(q+nlines(ii)-1))=peaks{ii}(k,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
%         else
%             alldata(q:(q+nlines(ii)-1))=peaks{ii}(k,col(1))./peaks{ii}(k,col(2));
%         end
%         q=q+nlines(ii);
%     end
%     
% end
% 
% figure(4), hold on
% plot(alldata(:,k),'-*','color',colors{k});
% ylim([0 1.5]);
% end

%end