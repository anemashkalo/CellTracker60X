
function [alldata] = PlotCellTraces_AN(matfile, col,N)

%plot cell traces for specific colonies using the output from the runTracker
%EDS
% in that output last column of peaks = has the trajectory number that the
% cell belongs to
% in the colonies cell array( saved in the otfile) the last column is the
% number of colony within that frma that the cells belong to

pp=load(matfile,'peaks','NucMasks','colonies'); % AN
peaks=pp.peaks;
colonies = pp.colonies;






% pp=load(matfile,'peaks','NucMasks','colonies'); % AN
% peaks=pp.peaks;
% colonies = pp.colonies;
% nc = size(colonies,2);% number of colonies found in the image NEED TO GET THE MAX NUMBER OF COLONIES IN THE IMAGE
% alldata = zeros(length(colonies),nc);
% nlines=zeros(length(colonies),nc);
%     for k=1:nc
% for ii=1:length(colonies)
%     if ~isempty(colonies{ii})&& colonies{ii}(k).ncells == N
%                 nlines(ii,k)=size(colonies{ii}(k),1);
%     end
% end
%     end
% alllines = sum(nlines);
% alldata=zeros(alllines(1),nc);
% for k=1:nc
% q=1;
% for ii=1:length(colonies)
%     if ~isempty(colonies{ii})
%          if length(col)==1
%             alldata(q:(q+nlines(ii)-1),k)=colonies{ii}(k).data(:,col);%make a single column vector from all the data (normalized intensity of col.6 in peaks to dapi (col. 5) in peaks
%         else
%             alldata(q:(q+nlines(ii)-1),k)=mean(colonies{ii}(k).data(:,col(1)))./mean(colonies{ii}(k).data(:,col(2)));
%         end
%         q=q+nlines(ii);
%     end
%     end
% end


end

