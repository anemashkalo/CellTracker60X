function [datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,img,flag)

ff=readAndorDirectory(direc);
chan = ff.w;

global userParam;
userParam.gaussRadius = 10;
userParam.gaussSigma = 3;
userParam.small_rad = 3;
userParam.presubNucBackground = 1;
userParam.backdiskrad = 300;
userParam.colonygrouping = 120;
areanuclow = 1300;
areanuchi = 9000;

info = h5info(ilastikfile);
info.Datasets;   
infocyto = h5info(ilastikfilecyto);
infocyto.Datasets ;

% check these to make sure the dataset name is correct for the 'h5read' function

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

datacyto = h5read(ilastikfilecyto,'/exported_data');
datacyto = squeeze(datacyto);
%time = size(data,3);   
Lnuc = [];
Lcytofin = [];
    
Lnuc = data(:,:,img) <2;  % <2 need to leave only the nuclei masks and make the image binary
Lnuc =  bwareafilt(Lnuc',[areanuclow areanuchi]);

% while isempty(Lnuc) == 0
%       continue
%  end

stats = regionprops(Lnuc,'Centroid');
xy = [stats.Centroid];
xx = xy(1:2:end);
yy=xy(2:2:end);

vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);

filename = getAndorFileName(ff,pos,ff.t(img),ff.z(zplane),chan(2)); % has to be channel 2 since all the masks should be applied to the gfp channel
filename2 = getAndorFileName(ff,pos,ff.t(img),ff.z(zplane),chan(1)); % to get info from the nuc channel
I2 = imread(filename);
Inuc = imread(filename2);

LcytoIl = datacyto(:,:,img) >1; 
LcytoIl = (LcytoIl');
Lcytonondil = LcytoIl;
LcytoIl = imdilate(LcytoIl,strel('disk',3));


% mechanism to remove random relatively large stuff from cyto channel (if
% it's size is comparable to the size of the cyto area of a small cell
% and couldn't be removed by area filtering

cc = bwconncomp(LcytoIl+Lnuc & ~vImg);%& ~vImg
cnuc = bwconncomp(Lnuc);
st = regionprops(cc,'PixelIdxList');
stnuc = regionprops(cnuc,'PixelIdxList');
goodinds = zeros(length(st),1);

for i = 1:length(stnuc)
    x =stnuc(i).PixelIdxList;
    for k=1:length(st)
        y = st(k).PixelIdxList;
        in = intersect(x,y);
        if ~isempty(in)
            goodinds(k,1) = k;
        end
    end
end
goodindsfin = nonzeros(goodinds);
goodstats = struct();

for i=1:length(goodindsfin)
                goodstats(i).PixelIdxList = st(goodindsfin(i)).PixelIdxList ;
end

 
% here need to leave the PixelIds of the goodinds and then convert back to the binary image by using ind2sub subtract the
% Lnuc to get the final mask of the cytoplasms
%
onebiglist = cat(1,goodstats.PixelIdxList);
Inew = zeros(1024,1024);
Inew(onebiglist) = true;
% open the image a little bit

LcytoIl = (Inew & ~ Lnuc & ~vImg);    % cyto masks initially include both nuclei+cyto, so need to eliminate nuc, use voronoi to divide;
% return back to the non-dilated cyto masks
LcytoIl = bwlabel(LcytoIl,8); 
LcytoIl(Lcytonondil ==0)=0;
Lcytofin = LcytoIl; 



% at this point should have an array of nuc and cyto masks(from ilastik
% and watershed respectively)

I2proc = imopen(I2,strel('disk',userParam.small_rad)); % remove small bright stuff
I2proc = smoothImage(I2proc,userParam.gaussRadius,userParam.gaussSigma); %smooth
I2proc = presubBackground_self(I2proc);


%get the NUCLEAR mean intensity for each labeled object
cc_nuc = bwconncomp(Lnuc,8);
statsnuc = regionprops(cc_nuc,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
statsnucw0 = regionprops(cc_nuc,Inuc,'Area','Centroid','PixelIdxList','MeanIntensity');% these are the stats for the actual nuclear image(rfp)

badinds = [statsnuc.Area] < 1000; 
badinds2 = [statsnucw0.Area] < 1000;
statsnuc(badinds) = [];
statsnucw0(badinds2) = [];


%get the cytoplasmic mean intensity for each labeled object
%cc_cyto = bwconncomp(Lcytofin);
statscyto = regionprops(Lcytofin,I2proc,'Area','Centroid','PixelIdxList','MeanIntensity');
badinds = [statscyto.Area] < 1000; % 
statscyto(badinds) = [];


% ncells = length(statsN);
 xy = stats2xy(statsnucw0);
 nuc_avrw0  = [statsnucw0.MeanIntensity]';%[statsN.NuclearAvr];
 nuc_areaw0  = [statsnucw0.Area]';%[statsN.NuclearArea];
 nuc_avrw1 = [statsnuc.MeanIntensity]';
 cyto_xy  = stats2xy(statscyto);
 cyto_area  = [statscyto.Area]';
 cyto_avrw1  = round([statscyto.MeanIntensity]');
 placeholder = -round(ones(length(xy(:,1)),1));
 % 

if isempty(statscyto)
    cyto_area = zeros(length(nuc_avrw1),1);
    cyto_avrw1 = cyto_area;
end

datacell=[xy(:,1) xy(:,2) nuc_areaw0 placeholder nuc_avrw0 nuc_avrw1 cyto_avrw1];%cyto_area
     

  if flag == 1
      figure, subplot(1,3,1),imshow(I2proc,[]);hold on
      subplot(1,3,2),imshow(Lcytofin);hold on
      subplot(1,3,3),imshow(Lnuc);hold on
  end
 
  
end



