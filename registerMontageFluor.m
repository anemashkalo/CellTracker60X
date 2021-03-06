%%
direc2 = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/lastTlivecellset2/');
[acdapi, fi]=alignManyPanelsAndorZstackMontage(direc2,[10 4],0,1);% here need to get the fluor, corresponding to the last time points 
[ac2dapi, fi2]=alignManyPanelsAndorZstackMontage('reimageFixedwithDAPIfluor_20160923_73926 PM',[7 11],0,1);% reimaged with dapi

load('tformwithDAPI.mat');% with DAPI
mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');

f2 = fi2;
f1reg = imwarp(fi,mytform);

f1reg_trans = zeros(size(f2));

 
toTranslate = [10,210];   
f1reg_trans =zeros(size(f2,1),size(f2,2));
f1reg_trans(1:(size(f1reg,1)-toTranslate(2)),1:(size(f1reg,2)-toTranslate(1))) = f1reg(toTranslate(2):size(f1reg,1)-1,toTranslate(1):size(f1reg,2)-1);
f1reg_trans2 = zeros(size(f1reg_trans));
img2output = cat(3,uint16(f1reg_trans),uint16(f2),zeros(size(f2)));
%imwrite(img2output,'testwithDAPI.tif');


% red - live cell 4x10
% green - fixed, 11x7
save('alignWithDapi','acdapi','ac2dapi');
%%
load('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/alignWithDapi.mat');
t = imread('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/testwithDAPI.tif');
%t = imread('test4mid.tif');
imshow(t(:,:,1),[]); % 2 - green chanel, 1 - red(4x10)
imcontrast
%%
% align the peaks from the old data again
%direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
direc = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-10-17-projections/nuclear/');

ff = dir(direc);
positions = (0:39);
strdir = '_out.mat';%_40X_imprBGandSegm
N = 0;
l = 100;
xyall = [];
for k=1:size(positions,2) % same as the loop over image numbers
    
    outfile = [ num2str(positions(k)) strdir ];
    
    load(outfile,'peaks');
    
    if ~isempty(peaks{l-N})       
        disp([int2str(k) ': ' int2str(size(peaks{l-N},1)) ' cells.']);
            
    xynow = bsxfun(@plus,peaks{l-N}(:,1:2),acdapi(k).absinds([2 1]));
    xyall = [xyall; xynow];
    end
end
%% translate and rotate peaks
load('tformwithDAPI.mat');% with DAPI
mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');

toTranslate = [-250,50];    

midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);

xyall_rot = bsxfun(@minus,xyall, mean(xyall)); % mean subtract;
xyall_rot = xyall_rot*rotmat; %rotate about center
xyall_rot = bsxfun(@plus,xyall_rot,midpoint([2 1])); %add back

xyall_rot_trans = bsxfun(@plus,xyall_rot,toTranslate);
%%
figure; imshow(t(:,:,1),[0 5e3]); hold on;
plot(xyall_rot(:,1),xyall_rot(:,2),'r.');
plot(xyall(:,1),xyall(:,2),'g.');
plot(xyall_rot_trans(:,1),xyall_rot_trans(:,2),'c.');

%%
% directory with segmented files
fixedsegm = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-09-28-SegmentedFixedData_withDapi';
%strdirfix = '_testBacktocoord';
strdirfix = '_fixedwithDAPI';
ff = dir(fixedsegm);
xyallfix = [];
positions2 = (0:76); % saved outfiles starting from 0


for k=1:size(positions2,2) % same as the loop over image numbers
    
    outfile = [ num2str(positions2(k)) strdirfix ];
    
    load(outfile,'peaks');
       
    if ~isempty(peaks{1})   % only one time point here    
        disp([int2str(k) ': ' int2str(size(peaks{1},1)) ' cells.']);
            
    xynowfix = bsxfun(@plus,peaks{1}(:,1:2),ac2dapi(k).absinds([2 1]));
    xyallfix = [xyallfix; xynowfix];
    end
end
%% plot coordinates of the segmented cells (after fixing)
figure(1), hold on; imshow(t(:,:,1),[0 5e3]); hold on;
plot(xyall_rot_trans(:,1),xyall_rot_trans(:,2),'c.');
figure(1), hold on;
plot(xyallfix(:,1),xyallfix(:,2),'r.');

%xyall rotated and translated peaks (from actual live cell data, last time point)
%xyallfix peaks, moved into the frame of the montge image

%% group the new peaks into colonies (fixed data)
% 
paramfile = 'setUserParamLiveImagingAN_40Xfixed';
run(paramfile);

for k=50:77
outfile = [ num2str(positions2(k)) strdirfix ];
 disp(outfile);
 
load(outfile,'peaks','colonies');

    clear colonies

nc = size(peaks{1},1);
pts = zeros(nc,2);
for ii=1:nc
    pts(ii,:) = (peaks{1}(ii,1:2));
end
allinds=NewColoniesAW(pts);
ngroups = max(allinds);
%ncells = cell(ngroups,1);
%Make colony structure for the single cell algorythm
for ii=1:ngroups;
    cellstouse=allinds==ii;
    colonies(ii)=dynColony(peaks{1}(cellstouse,:));
       
end
if ~isempty(colonies)
save(outfile,'colonies','-append');
end
end

%% check the colony in fixed dataset
%
dirr = fixedsegm;
outfile = '77_fixedwithDAPI';
colormap = colorcube;
load(outfile,'peaks','colonies','imgfiles')
Lnuc = imgfiles(1).NucMask;

figure(2), hold on ,imshow(Lnuc,[])
for k=1:size(colonies,2)
figure(2), hold on, plot(colonies(k).cells(:,1),colonies(k).cells(:,2),'.','color',colormap(k+5,:,:),'markersize',20);
hold on
end
%% plot colonies onto montage
figure(1), hold on
colormap = prism;
str = '_fixedwithDAPI';
for k=1:77
outfile = [ num2str(k-1) str ];
load(outfile,'peaks','colonies','imgfiles')
figure(1), hold on,
if ~isempty(colonies)
for j=1:size(colonies,2)
    
figure(1), hold on, plot(colonies(j).cells(:,1)+ac2dapi(k).absinds(2),colonies(j).cells(:,2)+ac2dapi(k).absinds(1),'.','color',colormap(j+8,:,:));
hold on
end
end
end

%% get the centroids of the fixed colonies
dirfix = '/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-09-28-SegmentedFixedData_withDapi';

colormap = prism;

fluordata = [];
fluordata1 = [];
fluordata2 = [];
fluordata3 = [];
allcell = [];
strfix = '_fixedwithDAPI.mat';
colfixall = [];% will have the centers of the colonies for the fixed data
for k=1:77
outfile = [ dirfix '/' num2str(k-1) strfix ];
load(outfile,'peaks','colonies','imgfiles')
figure(1), hold on,
if ~isempty(colonies)
for j=1:size(colonies,2)
    colfix(:,1) = mean(colonies(j).cells(:,1)+ac2dapi(k).absinds(2)); % mean of the x coord if the clony
    colfix(:,2) = mean(colonies(j).cells(:,2)+ac2dapi(k).absinds(1)); % mean of y
    colSZ = size(colonies(j).cells(:,1),1);
    fluordata1 = mean(colonies(j).cells(:,5));% dapi
    fluordata2 = mean(colonies(j).cells(:,6));% gfp
    fluordata3 = mean(colonies(j).cells(:,8));% cdx2
    colfixall = [colfixall;colfix]; % put data in the matrix
    fluordata = [fluordata;fluordata1 fluordata2 fluordata3 colSZ];
    %allcell = [allcell;colonies(j).cells(:,1)+ac2(k).absinds(2); colonies(j).cells(:,2)+ac2(k).absinds(1)];
%figure(1), hold on, plot(colonies(j).cells(:,1)+ac2(k).absinds(2),colonies(j).cells(:,2)+ac2(k).absinds(1),'.','color',colormap(j+8,:,:));hold on
figure(1), hold on, plot(colfix(:,1),colfix(:,2),'*','color',colormap(j+8,:,:),'markersize',10);hold on

end
end
end
% all data saved in the  matfile 'registeredDAPI'
%% get the centroids of the live data
load('registeredDAPI.mat','colfixall','xyall','fluordata');
load('alignWithDapi.mat');
load('tformwithDAPI.mat');

%dirlive = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/FinalOutfiles_completeTraces');
dirlive = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/2016-07-07-LiveCellTiling_28hr10ngmlBMP4/2016-10-17-projections/new_outfiles_tiling1');

figure(1), imshow(t(:,:,1),[0 5e3]); hold on;

mytform = fitgeotrans(movingPoints, fixedPoints, 'affine');
toTranslate = [-250,50];    
midpoint = floor(size(t)/2);
rotmat = mytform.T(1:2,1:2);

colormap = prism;
positions = (0:39);
strdir = '_out.mat';
N = 0;
l = 100;
q = 1;

for k=1:size(positions,2)
    outfile = [ dirlive '/' num2str(k-1) strdir ];
    load(outfile,'colonies')
    if ~isempty(colonies)
        for j=1:size(colonies,2)
            posnow = [];
            for h=2:size(colonies(j).cells,2) % start from second, since the first one is empty , as returned by the new analysis
                if colonies(j).cells(h).onframes(end) == (100)
                                %figure(q); hold on;
                    posnow = [posnow; colonies(j).cells(h).position(end,:)];
                end
            end
            if ~isempty(posnow)
            datatogether(q).colony = colonies(j);
            posnow = bsxfun(@plus,posnow,acdapi(k).absinds([2 1]));
            posnow = bsxfun(@minus,posnow,mean(xyall));
            posnow = posnow*rotmat;
            posnow = bsxfun(@plus,posnow,midpoint([2 1])+toTranslate);
            
            colnow_center = mean(posnow,1);
            dists = bsxfun(@minus,colfixall,colnow_center);
            dists = sqrt(sum(dists.*dists,2));
            [~, ind]=min(dists);
            datatogether(q).fixedData = fluordata(ind,:);
            figure(1), hold on
            plot(colnow_center(1),colnow_center(2),'.','color',colormap(j+10,:,:),'markersize',10);%
            plot(colfixall(ind,1),colfixall(ind,2),'*','color',colormap(j+10,:,:),'markersize',10);%
            q= q + 1;
            end
        end
        
    end
    
end

%save('registeredDAPInewTraces','datatogether','-append');
%% plot traces with fixed data
load('registeredDAPInewTraces.mat');
%load('registeredDAPI.mat');
positions = (0:39);
clear tracesbycol
last = 100;
trajmin = 40;%50
findata1 = [];
findata2 = [];
fr_stim = 16;%16
traces = [];
cdx2val = [];
colormap = parula;%;%
cdx2todapi = [];
tracesbybin = cell(1,2);% two bins, same colony size
binSZ = [1];%[0.8 1.1];
q = 1;
nc = 1;

clear dat
for k=1:size(datatogether,2)   %  loop over colonies
    %colSZ = datatogether(k).colony.numOfCells(fr_stim);            % old segmentation
   if size(datatogether(k).colony.ncells_actual,1)>fr_stim;          % new 
    % disp(datatogether(k).colony.ncells_actual);
       colSZ =datatogether(k).colony.ncells_actual(fr_stim) ;% new segmentation            % new 
    cdx2todapi(k) = datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1);
    if colSZ == nc && (cdx2todapi(k)<= binSZ(1)) 
       traces{k} = datatogether(k).colony.NucSmadRatio;%                                   % new 
       %traces{k} = datatogether(k).colony.NucSmadRatioOld;                                % old
     
        traces{k}((traces{k} == 0)) = nan;
        sz = size(traces{k},2);
        for h = 1:size(traces{k},2)
            [r,~] = find(isfinite(traces{k}(:,h)));
             dat = zeros(size(traces{k},1),1);
             dat(r,1) = traces{k}(r,h);
            if length(nonzeros(dat))>trajmin 
               disp(['filter trajectories below' num2str(trajmin)]);
               disp(['use' num2str(length(nonzeros(dat)))]);
               figure(1), plot(dat,'-*','color','b');hold on     % here plot the traces
               
                tracesbybin{1}(:,q+sz-1) = dat;      % traces{k}(:,h) here store the traces which meat condition
                cdx2val{1}(:,q+sz-1)= datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1);
                % disp(q+sz-1)
            end
            
        end
        
        q = q+sz;
         text(datatogether(k).colony.cells(h).onframes(end),traces{k}(end,h),num2str(datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1)),'color','b','fontsize',11);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
        figure(1) ,hold on
    ylim([0 2.5]);
    xlim([0 115])
    ylabel('mean Nuc/Cyto smad4  ');
    xlabel('frames');
%     findata1 = [findata1 ; datatogether(k).fixedData(1,3)/datatogether(k).fixedData(1,1)];
%     findata2 = [findata2; colSZ];
    end
    
    if colSZ == nc && (cdx2todapi(k)> binSZ(1))
       traces{k} = datatogether(k).colony.NucSmadRatio;                                     % new
       %traces{k} = datatogether(k).colony.NucSmadRatioOld;
        traces{k}((traces{k} == 0)) = nan;
        sz = size(traces{k},2);
        for h = 1:size(traces{k},2)
            [r,~] = find(isfinite(traces{k}(:,h)));
             dat = zeros(size(traces{k},1),1);
             dat(r,1) = traces{k}(r,h);
            if length(nonzeros(dat))>trajmin 
                figure(2), plot(dat,'-*','color','r');hold on     % here plot the traces
                
                tracesbybin{2}(:,q+sz-1) = dat;      %datatogether(k).colony.NucSmadRatio here store the traces which meat condition
                cdx2val{2}(:,q+sz-1)= datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1);
                %disp(q+sz-1)
            end
            
        end
        
        q = q+sz;
         text(datatogether(k).colony.cells(h).onframes(end),traces{k}(end,h),num2str(datatogether(k).fixedData(:,3)/datatogether(k).fixedData(:,1)),'color','r','fontsize',11);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
        figure(2) ,hold on
    ylim([0 2.5]);
    xlim([0 115])
    ylabel('mean Nuc/Cyto smad4  ');
    xlabel('frames');
    end
   
    
  end                                      % new 
end

%save('registeredDAPI','tracesbybin','binSZ','-append');%% average the trajectories that end up with variaous values of Cdx2
%% 
%load('registeredDAPI.mat','tracesbybin','tracesbybin2','binSZ','datatogether');
load('registeredDAPInewTraces.mat','tracesbybin2','binSZ','datatogether');

% tracesbybin      cell array that contains traces for the one cell
% colonies but separated into cel arrays according to the final Cdx2 value
% that they have 
% binSZ defines the Cdx2 range

%tracesbybin = tracesbybin2;% for the twocellcolonies

vect = (1:100)';
%trajmin = 50;
binmean = zeros(100,2);
err =zeros(100,2); 

for j =1:2                  % remove Nans
for k=1:size(binmean,1)
    for jj=1:size(tracesbybin{j},2)
   if (isfinite(tracesbybin{j}(k,jj))==0) || tracesbybin{j}(k,jj)>1.55  % to remove the noie due to very hugh signal probably in junk that was floting around at least trajmin
       tracesbybin{j}(k,jj) = 0;
   end
    end
end
end

 % average over cells
for j =1:2% loop over the bins;         
for k=1:size(binmean,1)
    
    binmean(k,j) = mean(nonzeros(tracesbybin{j}(k,:)));
    err(k,j) = std(nonzeros(tracesbybin{j}(k,:)));
end
end
%save('registeredDAPI','binmean','-append');
%% plot averaged trajectories
load('registeredDAPInewTraces.mat','tracesbybin','tracesbybin2','binmean','binmean2','err','err2','binSZ','datatogether');
vect = (1:100)';
colormap = colorcube;
C = {'c','r'};
label = {'CDX2 below','CDX2 above'};
b = [binSZ];
nc = 1;
if nc == 2
    binmean = binmean2;
    err = err2;
end
for j = 1:2
figure(5), errorbar(binmean(:,j),err(:,j),'-.','color',C{j},'linewidth',1.5); hold on%colormap(j+5,:)
ylim([0 2.5]);
xlim([0 105])
ylabel('mean Nuc/Cyto smad4  ');
xlabel('frames');
text(vect(end)+0.2*j,binmean(end,j)+0.1*j,[ label(j) num2str(b(1))],'color',C{j},'fontsize',20);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
title('One-cell colonies','fontsize',20);
if nc == 2
title('Two-cell colonies','fontsize',20);
end


figure(6), plot(vect,binmean(:,j),'-.','color',C{j},'linewidth',2); hold on%colormap(j+5,:)
ylim([0.4 1.8]);
xlim([0 115])
ylabel('mean Nuc/Cyto smad4  ');
xlabel('frames');
text(vect(end)+0.2*j,binmean(end,j)+0.1*j,[ label(j) num2str(b(1))],'color',C{j},'fontsize',20);%['mean ColCdx2 ' num2str(colonies2(j).cells(h).fluorData(1,end))]
title('One-cell colonies','fontsize',20);
if nc == 2
title('Two-cell colonies','fontsize',20);
end

end
%% extracs stats from the once cell traces

%load('registeredDAPI.mat','tracesbybin','cdx2val','tracesbybin2','binSZ','datatogether');
load('registeredDAPInewTraces.mat','tracesbybin','tracesbybin2','binmean','binmean2','err','err2','binSZ','datatogether','cdx2val');

base = [];
meanaft = [];
jump = [];
tracesd = [];
signlast = [];
fr_stim = 16;
tresp = 6;
last = [100:100];

for j =1:size(binmean,2)% loop over the bins;         
for k=1:size(tracesbybin{j},2)
    
    base(k,j) = mean((tracesbybin{j}(1:fr_stim,k)));
    meanaft(k,j) = mean((tracesbybin{j}(fr_stim:end,k)));
    jump(k,j) = abs(tracesbybin{j}(fr_stim+tresp,k)-tracesbybin{j}(fr_stim,k));
    tracesd(k,j) = std((tracesbybin{j}(fr_stim:end,k)));
    signlast(k,j)=(tracesbybin{j}(last(1),k));
end
end
%% plot stats
%load('registeredDAPI.mat','base','cdx2val','meanaft','jump','tracesd','signlast');

C = {'b','g','r'};
sym = {'.','*','d','s','o'};
a = 300;

toplot = signlast ;%tracesd base jump signlast;

for k=1:size(binmean,2)
   
        s =  size(nonzeros(toplot(:,k)),1);
        v = nonzeros(cdx2val{k}');
        figure(1),scatter(nonzeros(toplot(:,k)),v(1:s),a,C{k},sym{5});hold on
        
        xlim([0 1.8]);
        ylim([0 1.8]);
    
end
xlabel('Signaling at last time point');%Base level of signaling    mean after
ylabel('Cdx2');
box on
h = figure(1);
h.Children.FontSize = 18;
legend(['<' num2str(binSZ)],['>' num2str(binSZ)]);






