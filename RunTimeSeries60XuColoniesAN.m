function [peaks,dims,NucMasks,CytoMasks,colonies,imgfiles] = RunTimeSeries60XuColoniesAN(ilastikfile,ilastikfilecyto,pos,zplane,direc,flag)

info = h5info(ilastikfile);
info.Datasets;   

data = h5read(ilastikfile,'/exported_data');
data = squeeze(data);

time = size(data,3);   

peaks = cell(1,time); 
peakscyto = cell(1,time);
NucMasks = cell(1,time);
CytoMasks = cell(1,time);
imgfiles = struct([]);
for k=1:time
    
[datacell,Lnuc,Lcytofin] = IlastikplusWatershed_AW(ilastikfile,ilastikfilecyto,pos,zplane,direc,k,flag);

peaks{k} = datacell;

%colonies=peaksToMicroColoniesAN(peaks{k});
%plate1=plate(colonies,dims,[],ff.w,[],[], matfile);
imgfiles(k).compressNucMask = compressBinaryImg(Lnuc, size(Lnuc) );

NucMasks{k} = Lnuc;
CytoMasks{k} = Lcytofin;


dims = size(peaks);

end
colonies=peaksToMicroColoniesAN(peaks);% for each time frame % here the colonies is a cell array : each cell is a colony object

 %save('Outfile','peaks','NucMasks','CytoMasks','colonies');

 
end

