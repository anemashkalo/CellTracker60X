function [rawdata1] =  Intensity_vs_ColSize2(nms,nms2,dir,index1,param1,dapimax,chanmax,flag)
% this function is to plot the product of columns of peaks with indexes
% index1 : peaks(:,index1(1))*peaks(:,index1(2))

clear tmp
clear tmp2
clear rawdata
rawdata1 = cell(1,size(nms,2));
[dapimeanall,~] = getmeandapi(nms,dir,index1, dapimax);
for k=1:size(nms,2)
    filename{k} = [dir filesep  nms{k} '.mat'];
    load(filename{k},'plate1');
    colonies{k} = plate1.colonies;
    if ~exist('plate1','var')
        [colonies{k}, ~]=peaksToColonies(filename);
    end
    M(k) = max([colonies{k}.ncells]);
end

for k=1:size(nms,2)
    M = max(M);
    totalcolonies = zeros(M,1);
    rawdata = zeros(M,1);
    tmp2 = zeros(M,1);
    totalcells=zeros(M,1);
    
    col = colonies{k};
    
    for ii=1:length(col)
        a = any(col(ii).data(:,3)>dapimax(1));%%
        b = any(col(ii).data(:,index1(2))>chanmax);
        if ~isempty(col(ii).data) && a==0 && b==0 ;
            nc = col(ii).ncells;
            
            totalcolonies(nc)=totalcolonies(nc)+1;
            % totalcells(nc)=totalcells(nc)+nc;
            tmp = col(ii).data(:,index1(1)).*col(ii).data(:,index1(2)); %assign the value of product to tmp;
            tmp2(nc) = tmp2(nc) + sum(tmp); % add the elements tmp, corresponding to the same colony size, into the tmp2
        end
        
        
    end
    
    for l=1:length(totalcolonies)
        
        totalcells(l)=totalcolonies(l)*l;
    end
    
    for j=1:length(tmp2)
        rawdata(j) = tmp2(j)./totalcells(j); % average intensity of expression ( devide by the total number of cells of each colony size)
    end
    
    if flag == 1
    figure(6);subplot(1,size(nms2,2),k),  plot(rawdata(~isnan(rawdata)),'b*','markersize',15,'linewidth',2); legend(nms2{k});
    xlabel('Colony size');
    ylabel(['Expression of ',(param1),'marker']);
    xlim([0 8]);%size(rawdata(~isnan(rawdata)),1)
    end
    rawdata1{k} = rawdata;
end
end