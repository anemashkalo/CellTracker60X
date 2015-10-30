
function PlotCellTraces_AN(peaks, col)

% plot lie cell analysis results
ncell = size(peaks{10},1);
%colors = colorcube(ncell);

colors = {'r','g','b','m'};

nlines=zeros(length(peaks),1);
for ii=1:length(peaks)
    nlines(ii)=size(peaks{ii},1);
end
alllines = sum(nlines);
alldata=zeros(alllines,1);
for k=1:ncell
q=1;

for ii=1:length(peaks)
    if ~isempty(peaks{ii})
        if length(col)==1
            alldata(q:(q+nlines(ii)-1),k)=peaks{ii}(k,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
        else
            alldata(q:(q+nlines(ii)-1),k)=peaks{ii}(k,col(1))./peaks{ii}(:,col(2));
        end
        q=q+nlines(ii);
    end
    
end

figure(4), hold on
plot(alldata(:,k),'-*','color',colors{k});
end
end