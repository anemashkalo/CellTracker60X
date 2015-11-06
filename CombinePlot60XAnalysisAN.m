% combined analysis of 60X images so far
% frame 0: two one-cell colonies
% frame 19: two 0ne-cell colonies and one 2-cell colony
% frame 10: three one-cell colonies 
% frame 1: two 2-cell colonies 
% frame 13: one cell dividing into 2 (frames 12-16 show division)

dir = '.';
load('Frame00_Analysis.mat');
fr0col1 = cell2mat(onecellcol1)';
fr0col2 = cell2mat(onecellcol2)';
load('Frame10_Analysis.mat');
fr10col1 = cell2mat(onecellcol1)';
fr10col2 = cell2mat(onecellcol2)';
fr10col3 = cell2mat(onecellcol3)';
load('Frame19_Analysis.mat');
fr19col1 = cell2mat(onecellcol1)';
fr19col2 = cell2mat(onecellcol2)';
% load('Frame13_Analysis.mat');
% for i=1:11
% fr13col1{i} = nucmeanInt{i}./cytomeanInt{i};
% end
% for i=17:35
% fr13col2{i} = nucmeanInt{i}./cytomeanInt{i};
% end

   vect = 1:length(fr0col1);
    plot(vect',fr0col1,'r--*'),hold on
   vect = 1:length(fr0col2); 
    plot(fr0col2,'g--*'),hold on
   vect = 1:length(fr10col1); 
    plot(vect,fr10col1,'b--*'),hold on
   vect = 1:length(fr10col2); 
    plot(vect,fr10col2,'y--*'),hold on
   vect = 1:length(fr10col3); 
    plot(vect,fr10col3,'m--*'),hold on
    vect = 1:length(fr19col1);
    plot(vect,fr19col1,'k--*'),hold on
    vect = 1:length(fr19col2);
    plot(vect,fr19col2,'b--*');hold on
    
ylim([0.1 1.9]);

xlabel('time, frames');
 ylabel('nuc/cyto mean for cells in the frame');
 title('GFPsmad4RFPh2b cells, 10ng/ml bmp4 added after frame 16, ~ 9 hours imaging time');

 
 
% hold on
% plot(fr13col1,'*y');
% hold on
% plot(fr13col2,'*y');