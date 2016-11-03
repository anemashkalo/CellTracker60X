function [I2_bgsubtract] = simplebg(Lcytofin,Lnuc,I2)

if sum(sum(sum(Lcytofin))) == 0   % if there is no Lcytofin, make i an empty image of the same size as Lnuc
    Lcytofin = zeros(size(Lnuc,1),size(Lnuc,2),size(Lnuc,3));
end
allmask = Lcytofin | Lnuc;
allmask = imdilate(allmask,strel('disk',1));

%I2(I2>1700)=0; % make the junk also holes so it will be filled with local pixels as the cells themselves
imgNoCells = I2;
imgNoCells(allmask) = 0;

for k=1:size(Lnuc,3)
imgNoCells2(:,:,k) = imfill(imgNoCells(:,:,k));
end
%imshow(imgNoCells)
%imshow(imgNoCells,[])

I2_bgsubtract = imsubtract(I2,imgNoCells2);

%imshow(I2_bgsubtract(:,:,1),[])
end