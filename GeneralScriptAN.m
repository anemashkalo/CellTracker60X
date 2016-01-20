%%
ilastikfile = ('/Users/warmflashlab/Desktop/IlastikMasks_headless_PluriW0/NucMaskPluri_tg56.h5');
nuc = h5read(ilastikfile,'/exported_data');
  k = 12;
    nuc = nuc(2,:,:,k);% for probabilities exported
    nuc = squeeze(nuc);
    mask1 = nuc;
    
    imshow(mask1);
    Lnuc = im2bw(mask1,0.5);
    figure, imshow(Lnuc);


%%
ff = dir('*jan8set*.mat');
k = 1;
outfile = ff(k).name;
load(outfile);

%%
N  =5;
n = uncompressBinaryImg(imgfiles(N).compressNucMask);
nc = uncompressBinaryImg(imgfilescyto(N).compressNucMask);
%close all
figure,subplot(1,2,1),imshow(n);

hold on
subplot(1,2,2),imshow(nc);
%%
fr_stim = 22;
ff = dir('*jan8set*.mat');
for k=1:size(ff,1)
    outfile = ff(k).name;
    
runTracker(outfile,'newTrackParamAN');
global userParam;
userParam.colonygrouping = 120;
% look at colonies around the stimulation frame (window of couple hours)?
cellsToDynColonies(outfile);
end

%%
% colonies(2).numOfCells(1)
ff = dir('*jan8set*.mat');%jan8set % Pluri_42hrs % 22hr set: k = 3,4,6,16,22 traces k = 10 bad
k =6;
outfile = ff(k).name;
load(outfile);

fr_stim = [];%22
fldat = [2 3];
delta_t = 5; % 12% in minutes
p = fr_stim*delta_t/60;
flag = 1;
resptime = 36;% 36in frames ( converted to hours later)
GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);

%%
for j=1:size(ncells,1)
plot(ncells{j},'-*');hold on
ncells{j}(10)
end

%%
% plot the stats (before, after, amplitude...)
fr_stim = 22;%22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
flag = 1;
resptime =60;% 15 50 36 in frames ( converted to hours later)
coloniestoanalyze = 3; % max size of the colonies that want to look at
for jj = 1:coloniestoanalyze
    
    colSZ = jj;
    
    %GetDynamicColonyTraces(outfile,fr_stim,fldat,delta_t);
    %ff = dir('*10ngmlDifferentiated*.mat');
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs
    for k=1:size(ff,1)
        outfile = ff(k).name; %nms{k};
        load(ff(k).name);
        datafin = GetDynamicColonyStats(outfile,fr_stim,delta_t,flag,colSZ,resptime,coloniestoanalyze);
        
    end
end
if ~isempty(fr_stim)
figure(11)
xx = 0:0.1:2.5;
plot(xx,xx,'-k','linewidth',2);

end
if isempty(fr_stim) 
 figure(11) 
 xx = 0:1:coloniestoanalyze;
 yy = ones(1,coloniestoanalyze+1);
 plot(xx,yy,'--');
 figure(12)
 plot(xx,yy,'--');
end
%%
% get new traces
fr_stim = 22;%22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
userParam.colonygrouping = 120;
flag = 1;
colSZ = 3;
resptime =50;% 15 50 36 in frames ( converted to hours later)
p2 = resptime*delta_t/60;
coloniestoanalyze = 2;
cmap = parula;
C = {'b','g','r','m'};
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs
    %for %k=1:size(ff,1);
    k=5;

        outfile = ff(k).name; %nms{k};
        cellsToDynColonies(outfile);
        load(outfile,'colonies');
        
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        for j = 1:numcol
           % if colSZ == colonies(j).numOfCells(fr_stim); 
            figure(j+10), plot(colonies(j).NucSmadRatio(:),'-*','color',C{j});% cmap(j,:) 'r'traces
            ylim([0 2]);
            ylabel('mean Nuc/Cyto smad4  ');
            xlabel('frames');
        %end
        end
   % end 
%%
% NEW: get the before and after  plots

fr_stim = 22;%22 %38
fldat = [2 3];
delta_t = 12; % 12% in minutes
p = fr_stim*delta_t/60;
global userParam;
userParam.colonygrouping = 120;
flag = 1;

resptime =60;% 15 50 36 in frames ( converted to hours later)
p2 = resptime*delta_t/60;
coloniestoanalyze = 3;
cmap = parula;
C = {'b','g','r','m'};

for jj = 1:coloniestoanalyze
    
    colSZ = jj;
    ff = dir('*jan8set*.mat');%jan8set 10ngmlDifferentiated_22hrs % Pluri_42hrs
    for k=1:size(ff,1)
        
        outfile = ff(k).name; %nms{k};
        %cellsToDynColonies(outfile);
        load(outfile,'colonies','peaks');
        t = length(peaks);
        numcol = size(colonies,2); % how many colonies were grouped within the frame
        for j = 1:numcol
            if colSZ == colonies(j).numOfCells(fr_stim); % how many cells within colony at time fr_stim
                
                stats =  colonies(j).DynNucSmadRatio(t,fr_stim,resptime);%,resptime
                bf = (stats(:,1));
                aft =(stats(:,2));
                
                b = find(isnan(bf));
                bf(b)=0;
                a = find(isnan(aft));
                aft(a)=0;
                hold on,figure(4), plot(bf,aft,'*','color',C{colSZ},'markersize',15);
                hold on,figure(5),subplot(1,2,1), plot(bf,'*','color',C{colSZ},'markersize',15);
                ylim([0 2]);
                hold on,figure(5),subplot(1,2,2), plot(aft,'*','color',C{colSZ},'markersize',15);
                               
                
            end
            
        end
        
        
    end
end
                xx = 0:1:14;
                yy = ones(1,15);
                hold on,figure(4),plot(xx,xx,'-k','linewidth',2);
                ylim([0 2]);
                xlim([0 2]);
                ylabel(['mean Nuc/Cyto smad4  ' num2str(p2) ' hours after stimulation']);
                xlabel('mean Nuc/Cyto smad4 Before stimulation');
                hold on,figure(5),subplot(1,2,2),plot(xx,yy,'-k','linewidth',2);
                ylim([0 2]);
                hold on,figure(5),subplot(1,2,1),plot(xx,yy,'-k','linewidth',2);
                 ylim([0 2]);
%                 xlim([0 2]);
      
      
      
