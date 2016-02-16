
%function [maskz] = 3Dlblmask(CC,nucleilist)

%%

[PILsn,PILsSourcen, masterCCn, stats, nucleilist, zrange,CC] = traceobjectsz(smasks, userParam.matchdistance, zrange, userParam.zmatch);
%%

size(CC,2); % number of zplanes
size(nucleilist,2); % oved how many plane the niclei are spread
size(nucleilist,1); % howmany objects were found in plane 1
badind = cellfun(@isempty,CC);
CC(badind) = [];

 %%
   % from AW
   % get the masks based on nuclei list
  for ii = 1:length(CC)
      newplane = zeros(1024);
      for jj = 1:length(CC{ii}.PixelIdxList)
            newplane(CC{ii}.PixelIdxList{jj}) = jj;
      end
      trymask(:,:,ii) = newplane;
  end
  %  plot these
for k =1:4
figure(1),subplot(2,2,k); showMaskWithNumber(trymask(:,:,k));
end
  newmask = getGoodMask(trymask,nucleilist);
for k =1:4
figure(2),subplot(2,2,k); showMaskWithNumber(newmask(:,:,k));
end
  
 %%
 % get the mean intensitied from regionprops
 % newmask = the 3dlabel mask
 stats3d = regionprops(newmask,inuc(:,:,1:4),'MeanIntensity','Centroid');
 xyz_fin = cat(1,stats3d.Centroid);
% stats = regionprops(Lnuc,'Centroid','PixelIdxList');   % from 2D code
% xy = [stats.Centroid];
% xx = xy(1:2:end);
% yy=xy(2:2:end);

%vImg = mkVoronoiImageFromPts([xx' yy'],[1024 1024]);
%[V,C] = voronoin(X)
[V,C] = voronoin(xyz_fin);
K = convhulln(V);
trisurf(K,V(2:end,1),V(2:end,2),V(2:end,3));


%cc = bwconncomp(LcytoIl+Lnuc & ~vImg);    % from 2D code
 
 
 
 
 
 %%
% global userParam;
% maskz = zeros(1024,1024,size(zrange,2));
% Inew = zeros(1024,1024);
% 
% for k=1:size(nucleilist,2)
%     bw = regionprops(CC{k},'PixelIdxList')
%     for xx = 1:size(bw,1)
%         
%        % if isfinite(nucleilist(xx,k))
%             ind = find(~isnan(nucleilist(:,k)));
%             Inew(bw(xx).PixelIdxList) =nucleilist(ind(xx),k); % label of the object as it is in the nuclei list, same label for this object in different planes
%             
%         %end
%     end
%     maskz(:,:,k) = Inew;
%     Inew = zeros(1024,1024);
% end
% %%
% for l=1:size(maskz,3)
%     figure, imshow(maskz(:,:,l),[])
% end
% 
% 
% %%
% 
% for jj=1:size(zrange,2)
% z1 = regionprops(maskz(:,:,jj),'PixelIdxList');
% numobj(jj) = size(z1,1);
% 
% end
% maxobj = max(numobj); % number of object found in all planes, will be the number of separate labels
% maxobjN=find(numobj==maxobj);
% % now need to get the intersections of all the elements within the
% % mask(:,:,maxobjN(1)) with all the elements of the remaining masks and
% % assign same labels to the intersecting pixels, if not matched with
% % anything, give it a new label
% lbl = 1:maxobj;
% match = cell(maxobj,1);
% ref = regionprops(maskz(:,:,maxobjN(1)),'PixelIdxList','Centroid'); % get the pixel id listis of the mask where the max number of separate objects was found
% for i=1:maxobj
%     
%     for j=1:size(zrange,2)
%         if j ~=maxobjN(1)     % get the intersection of each object in ref with each object in the remaining planes, excluding the comparison btw the same planes
%             test = regionprops(maskz(:,:,j),'PixelIdxList','Centroid');
%             for h=1:size(test,1)
%                 [overlap] = intersect(ref(i).PixelIdxList,test(h).PixelIdxList);
%                 if ~isempty(overlap)
%                     xy1 = [ref(i).Centroid];
%                     xy2 = [test(h).Centroid];
%                     D = sqrt(power((xy1(1)-xy2(1)),2)+power((xy1(2)-xy2(2)),2));
%                     if D < userParam.matchdistance;
%                         
%                         match{i}(j,1)= j;
%                         match{i}(j,2)=h;
%                         match{i}(j,3)=D;
%                         %match{i}(j,3)=overlap;
%                         
%                     end
%                     %                     else
%                     %                     match{i}(:,:) = 0;
%                 end
%             end
%             
%         end
%     end
% end
% %%
%    lbl = 1:maxobj;
%    maxobjN(1); % reference plane ( with max number of objects)
%    ref = regionprops(maskz(:,:,maxobjN(1)),'PixelIdxList','Centroid');
%    maskznew = zeros(1024,1024,size(zrange,2));
%    Inew = zeros(1024,1024);
%    Inew2 = zeros(1024,1024);
%    
%    for i = 1:size(match,1)
%    Inew(ref(i).PixelIdxList)=lbl(i);
%    maskznew(:,:,maxobjN(1)) = Inew;
%    lbl(i)
%    %Inew = zeros(1024,1024);
%    
%    inotherplane = nonzeros(match{i}(:,1));
%    pxlinotherplane = nonzeros(match{i}(:,2));
%    
%    for j=1:size(inotherplane,1)
%        otherplane = regionprops(maskz(:,:,inotherplane(j)),'PixelIdxList','Centroid');
%        Inew2(otherplane(pxlinotherplane(j)).PixelIdxList)=lbl(i);
%        maskznew(:,:,inotherplane(j)) = Inew2;
%        lbl(i)
%        inotherplane(j);
%       % Inew2 = zeros(1024,1024);
%    end
%    
%    end
   
  
      