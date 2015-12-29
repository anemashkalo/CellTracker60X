function peaks = nucCytoIlastik2peaksLoop(ilastikDirec,imageDirec,zplane,pos,chan,paramfile,outfile)

% dir - has the ilastik masks, the  .h5 files
% dir = ('/Users/warmflashlab/Desktop/A_NEMASHKALO_Data_and_stuff/9_LiveCllImaging/Nov12ImaginfResults/');
% zplane =  need to specify which image in z plane from the raw data to take
% direc =  directory with actual images from live-cell imaging experiment
% pos = the number of frame processed TO DO: need to loop over those
% dt =  data type: if 1 - means that all the deta is separate and the code
% will not use bioformats to extract the images
% if td = 0, the n the code will decompose gtouped data frames and extract
% the correct images
% tg = time group, this should be a vector, its length should be the
% number of time froups the data was devided into upon saving
% TO DO:

[~, ilastikCytoAll]=folderFilesFromKeyword(ilastikDirec,'Cyto',['{00' num2str(pos) '}']);% all cyto masks for the frame 0 ( four time groups)'Cyto','{0002}'['Outfile_000' num2str(pos) '_t']
[~, ilastikNucAll]=folderFilesFromKeyword(ilastikDirec,'Nuc',['{00' num2str(pos) '}']);% all nuc masks for the frame 0 ( four time groups)

for j = 1:length(ilastikCytoAll)
    
    %get the ilastik masks
    ilastikNuc = fullfile(ilastikDirec,ilastikNucAll(j).name);
    ilastikCyto= fullfile(ilastikDirec,ilastikCytoAll(j).name);
    
    nuc_mask_all = h5read(ilastikNuc,'/exported_data');
    nuc_mask_all = squeeze(nuc_mask_all);
    
    cyto_mask_all = h5read(ilastikCyto,'/exported_data');
    cyto_mask_all = squeeze(cyto_mask_all);
    
    %get the image readers
    ff = readAndorDirectory(imageDirec);
    filename1 = getAndorFileName(ff,pos,0,zplane,chan(1));
    filename2 = getAndorFileName(ff,pos,0,zplane,chan(2));
    img_nuc_reader = bfGetReader(filename1);
    img_cyto_reader = bfGetReader(filename2);
    
    nT = img_nuc_reader.getSizeT;
    
    for k = 1:nT
        
        plane1 = img_nuc_reader.getIndex(0,0, k - 1) + 1;
        nuc_img = bfGetPlane(img_nuc_reader,plane1);
        
        plane1 = img_cyto_reader.getIndex(0,0, k - 1) + 1;
        nuc_cyto = bfGetPlane(img_cyto_reader,plane1);
        
        outdat = nucCytoIlastik2peaks(nuc_mask_all(:,:,k),cyto_mask_all(:,:,k),nuc_img,nuc_cyto,paramfile);
        
        peaks{k+(j-1)*nT} = outdat;
        if k == 1 || mod(k,10) == 0
            save(outfile,'peaks');
        end
    end
    
    save(outfile,'peaks');
end
