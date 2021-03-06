#!/bin/bash
# ilastik segmentation output to peaks: input- multiple time points, multichannel, multitiff files.
#ndir: new directory where segmentation results are saved
#projectpath: path corresponding to ilastik project that was usecd to train a sample dataset
#savepath: path corresponding to ilastik output export directory
#segchannel: channel used for segmenting images
#samplepath: path corresponding to sample images that are yet to be segmented - serves as an input to the matlab function
#files: image files in the samplepath that are to be segmented 


i=0
num=1

ndir=/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2016-09-08-cytomask

projectpath=/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2016-09-08-Retrain/Retraining_cyto.ilp

savepath=/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2016-09-08-cytomask/jan8datacytomask_11.h5

#segchannel=1

#mkdir $ndir

samplepath=/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2D_ToProcess_Jan8Data/W1/;

files=/Users/warmflashlab/Desktop/JANYARY_8_DATA_ilasik/2D_ToProcess_Jan8Data/W1/*;

for f in $files
do
/Applications/ilastik-1.1.8-OSX.app/Contents/MacOS/python /Applications/ilastik-1.1.8-OSX.app/Contents/Resources/ilastik.py --headless  --project=$projectpath $f

i=$((i+num))

scp $savepath ${ndir}cytomasks2D${i}.h5

echo $f $i
done

exit;



