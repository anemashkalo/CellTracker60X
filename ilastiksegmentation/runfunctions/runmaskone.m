function peaks = runmaskone(segfiledir, rawfiledir, mrnafilepath, paramfile, colonyno, objno)
%%
% segfiledir: directory path of ilastik 2d segmentation probability density maps
% rawfiledir: directory path of the nuclear channel raw images (the same
% that was fed to FISH)
% nzslices = no. of z slices
% colonyno: image no (1 colony implies I image or 1 unique position imaged)
% objno: cell no. for which you want to check no. of mRNA's assigned (for
% eg: if your raw image has around 20 different cells, then objno could be
% anywhere between 1 and 20.

global userparam
eval(paramfile)

nzslices = userparam.nzslices;

[dirinfo, start] = readdirectory(segfiledir);
[dirinfo1, start1] = readdirectory(rawfiledir);


maskno = start + (colonyno-1)*nzslices; % ilastik mask
imageno = start1 + (colonyno-1)*nzslices ; %nuclear channel raw image
[pnuc, inuc] = readmaskfiles(maskno, segfiledir, rawfiledir, dirinfo, dirinfo1, nzslices, imageno);

%%

pmasks = primaryfilter(pnuc,userparam.logfilter, userparam.bthreshfilter, userparam.diskfilter, userparam.area1filter);
%%

[zrange, smasks] = secondaryfilter(pmasks, userparam.minstartobj, userparam.minsolidity, userparam.diskfilter, userparam.area2filter);
%%
[PILsn,PILsSourcen, CC, masterCCn, stats, nucleilist, zrange] = traceobjectsz(smasks, userparam.matchdistance, zrange, userparam.zmatch);
%%
overlapthresh = 80;
colonyno = (maskno - start)/nzslices + 1;
imviews = 0;
[nucleilist, masterCC] =  overlapfilter(PILsn, PILsSourcen, masterCCn, nucleilist, inuc, zrange, userparam.overlapthresh, userparam.imviews);
%%
% nuclei3dmasks = masksin3d(CC, nucleilist, nzslices);

%%
clear peaks;
peaks{colonyno} = mrnapercells(nucleilist, stats, mrnafilepath, colonyno, zrange, userparam.channels, userparam.cmcenter);

%%
if (~exist('objno'))
    objno = [5];
end
nucleimrnacheck(masterCC, inuc, zrange, peaks, colonyno, objno, channels, mrnafilepath, cmcenter);
