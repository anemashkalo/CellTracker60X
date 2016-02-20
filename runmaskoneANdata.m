function [peaks,Lnuc,Lcytofin] = runmaskoneANdata(ilastikdir,ilastikdircyto, imagedir,pos,tpt, timegroup,chan,paramfile)

% ilastikdir: directory path of ilastik 2d segmentation probability density maps
% imagedir: directory path of the nuclear channel raw images 
% 
% obtain the ilastik prob masks from nuclear channel in the form of 3d
% image


%get the prob masks and raw images into the 3d format
[pnuc, inuc] = readmaskfiles1(ilastikdir,imagedir, pos,tpt, timegroup,chan(1));        % nuc
[pcyto, icyto] = readmaskfiles1(ilastikdircyto,imagedir, pos,tpt, timegroup,chan(2));  % cyto processing


setUserParam3DsegmentationAN;
setUserParamLiveImagingAN
global userParam;

% sapna code tracking 
pmasks = primaryfilter(pnuc,userParam.logfilter, userParam.bthreshfilter, userParam.diskfilter, userParam.area1filter);

% zrange: where the nuclei are in z
[zrange, smasks] = secondaryfilter(pmasks, userParam.minstartobj, userParam.minsolidity, userParam.diskfilter, userParam.area2filter);

% get the acual tracking nucleilist = tracked objects, labled 1-N,CC -
% pixelidxlist of all objects tracked in all planes
[PILsn,PILsSourcen, masterCCn, stats, nucleilist, zrange,CC] = traceobjectsz(smasks, userParam.matchdistance, zrange, userParam.zmatch);

% use nucleilist to relabel the tracked objects with unique labels
[newmask_lbl] = lblmask_3Dnuc(CC,nucleilist);

% [nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userParam.overlapthresh);

% leave only planes with cells in the cyto channel too

 for k=1:size(zrange,2)
     
 pmaskscyto(:,:,k) = im2bw(pcyto(:,:,k),userParam.probthresh_cyto);
 icyto_new(:,:,k) =icyto(:,:,zrange(k));                         % inuc, icyto = still have 5 layers instack, need to leave only the ones
 inuc_new(:,:,k) =inuc(:,:,zrange(k));                           % with the nuclei in them ( as was determined by zrange)
 
 end

 % get the labeled cyto masks ( nuclei not subtracted in the
 % GetVoronoiCells3D function since need them filled in order to remove
 % cytos without nuc)
 
 [maskzcyto_lbl] = GetVoronoiCells3D(newmask_lbl,pmaskscyto);
 
 % now maskzcyto and newmask are 3D masks that need to be applied to
 % icyto,inuc, which are all now of the same dim
 % all the masks are labelled
   
 [datacell,Lnuc,Lcytofin] = nucCytoIlastik2peaks_3Dsegm(newmask_lbl,maskzcyto_lbl,inuc_new,icyto_new,paramfile);
 peaks = datacell;
 
 
end




