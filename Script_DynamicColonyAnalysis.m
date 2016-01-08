% add the loop over positions
% save every 4th file as a separate position

zplane = [];

chan = [0 1];
[nums, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec,'CytoMask');%num2str(pos)]
timegroups = 4;
positions = length(nums)/timegroups; % number of separate position numbers (start from 0)
for jj = 1:positions
pos = jj-1;
outfile = ([ num2str(pos) '_' num2str(outfile)]);
end
paramfile = 'setUserParamLiveImagingAN';
outfile = 'testheadless.mat';% basic name for all positions
ilastikDirec = ('/Users/warmflashlab/Desktop/IlastikMasks_headlessW1');
imgDirec1 = ('/Users/warmflashlab/Desktop/MaxProjectionsLiveImg_DiffCondition(nov12data)/Nov12ImagingMaxProj_W0');% already max projections
imgDirec2 = ('/Users/warmflashlab/Desktop/MaxProjectionsLiveImg_DiffCondition(nov12data)/Nov12ImagingMaxProj_W1');% already max projections

peaks = nucCytoIlastik2peaksLoop(ilastikDirec,imgDirec1,imgDirec2,zplane,pos,chan,paramfile,outfile);% tsted
outfile = ([ num2str(pos) '_' num2str(outfile)]);
addShiftToPeaks(outfile,fr_stim);
runTracker(outfile,'newTrackParam');
global userParam;
userParam.colonygrouping = 120;
cellsToDynColonies(outfile);



fr_stim = 38;
fldat = [2 3];
delta_t = 5; 
p = fr_stim*delta_t/60;
colSZ = 1;
flag = 1;
 % add the loop oved positions here
GetDynamicColonyTraces(matfile,fr_stim,fldat,delta_t,colSZ);
datafin = GetDynamicColonyStats(matfile,fr_stim,delta_t,flag,colSZ);