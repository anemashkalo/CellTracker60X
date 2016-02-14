function peaks = runmaskoneANdata(ilastikdir,ilastikdircyto, imagedir,pos,tpt, timegroup,chan,paramfile)
%%
% ilastikdir: directory path of ilastik 2d segmentation probability density maps
% imagedir: directory path of the nuclear channel raw images 


[pnuc, inuc] = readmaskfiles1(ilastikdir,imagedir, pos,tpt, timegroup,chan(1));%
%%[pnuc, inuc] = readmaskfiles1(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);
%%
% 
setUserParam3DsegmentationAN;
global userParam;

pmasks = primaryfilter(pnuc,userParam.logfilter, userParam.bthreshfilter, userParam.diskfilter, userParam.area1filter);
%%

[zrange, smasks] = secondaryfilter(pmasks, userParam.minstartobj, userParam.minsolidity, userParam.diskfilter, userParam.area2filter);
%%

[PILsn,PILsSourcen, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, userParam.matchdistance, zrange, userParam.zmatch);
%%

[nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userParam.overlapthresh);
%%




% plot the distinct nuclei

lb = labelmatrix(masterCCn);
a = label2rgb(lb);
figure, imshow(a);hold on
bw = regionprops(lb,'Centroid','PixelIdxList','Area');
badinds = [bw.Area]< userParam.area2filter ;
bw(badinds) = [];


Inew = zeros(1024,1024);
for j=1:size(bw,1)
Inew(bw(j).PixelIdxList) = 1;
figure, imshow(Inew);
Inew = zeros(1024,1024);
end
%bw has correctly 7 cells, now need to put those objects into a labeled matrix again

N = size(stats,2);
color = colorcube;
C = {'r','g','b','k','y'};
for j=1:N
    if ~isempty(stats{j})
    s = stats{j};    
    xy = [s.Centroid];
    x = xy(1:2:end);
    y =  xy(2:2:end);
    plot(x,y,'*','color',C{j});%color(j,:)
    end
end
%%
%get the zplane numbers and coordinates of the unique nuclei 
figure(5), hold on
for x = 1:size(zrange,2)
    nucinzx = nucleicenter(round(nucleicenter(:,3))==zrange(x),1:3);
    
    if ~ isempty(nucinzx)
        plot(nucinzx(:,1),nucinzx(:,2),'*','color',C{x},'markersize',20);hold on%color(j,:)
    end
    
end
tmp = round(nucleicenter(:,3));
nucleicenter(:,3) = tmp;
nucleicenter_fin =nucleicenter;

goodslices =  unique(nucleicenter_fin(:,3));
size(goodslices,1);
% once the good planes are determined, can get the corresponding
% cytoplasmic masks and raw images

%ilastikdircyto = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/3Dsegmentation_tracking_TrainingSet/masks_zcyto');
[pcyto, icyto] = readmaskfiles1(ilastikdircyto,imagedir, pos,tpt, timegroup,chan(2));%

%%
% run the peaks separately on those slices (nuc and cyto) and combine into
% one peaks
paramfile = 'setUserParamLiveImagingAN';
zdata = cell(1,size(goodslices,1));
n = cell(1,size(goodslices,1));
c = cell(1,size(goodslices,1));
for s=1:size(goodslices,1);
mask1 = pnuc(:,:,goodslices(s))';
mask2 = pcyto(:,:,goodslices(s))';
img_nuc = inuc(:,:,goodslices(s))';
img_cyto = icyto(:,:,goodslices(s))';



[datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks(mask1,mask2,img_nuc,img_cyto,paramfile);
zdata{s} = datacell;
n{s} = Lnuc;
c{s} = Lcytofin;
end
for s=1:(size(goodslices,1)-1)
fulldata = cat(1,zdata{s},zdata{s+1});

end



