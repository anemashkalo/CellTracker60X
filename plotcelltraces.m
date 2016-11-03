 % plot cell tracesfor the specific outfile
 
 function plotcelltraces(outfile,trajmin)
 colormap = prism;
 C = {'b','r','g','m','b','c','g','r','b','g','m'};
 load(outfile,'colonies');
 % tps = length(peaks);
 numcol = size(colonies,2); % how many colonies were grouped within the frame
 traces = cell(1,numcol);
 for j = 1:numcol
     
     traces{j} = colonies(j).NucSmadRatio;%colonies(j).NucSmadRatio(:)
     
     traces{j}((traces{j} == 0)) = nan;
     for h = 1:size(traces{j},2)
         [r,~] = find(isfinite(traces{j}(:,h)));                  %
                dat = zeros(size(traces{j}(:,h),1),1);
                dat(r,1) = traces{j}(r,h);
         if length(nonzeros(dat))>trajmin
             disp(['filter trajectories below' num2str(trajmin)]);
             disp(['use' num2str(length(nonzeros(dat)))]);
             hold on;figure(j+11),plot(dat,'-o','color',C{j});hold on% 
             %text(traces{j}(:,h)+0.1,num2str(colSZ1(:)));
             ylim([0 2.5]);
             ylabel('mean Nuc/Cyto smad4  ');
             xlabel('frames');
         end
     end
 end
 end
             