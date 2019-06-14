### How to create a DICOM-aware workstation solution for image processing

DICOM is the standard data exchange format in a clinical workflow. Here is an example of how to read and write DICOM files that are compatible with most imaging systems. That is they can be imported back into the PACS system and will be assigned to the previous patient record.

There are two basic steps, reading in the data and writing the data back as a series of individual DICOM images.

#### Reading in a series of DICOM files

```
file_list = recursive_dir(masterdicom);
totfiles = length(file_list);

% we should sort by DICOM series first, for each series we should do...

totdicomfiles = 0;
file_list_dicom = {};
instancenumbervec = [];
slicelocationvec = [];
for i = 1:totfiles;
    try
        dcminfotmp = dicominfo(file_list{i});
        totdicomfiles = totdicomfiles+1;
        file_list_dicom{totdicomfiles} = file_list{i};
        instancenumbervec(totdicomfiles) = dcminfotmp.InstanceNumber;
        slicelocationvec(totdicomfiles) = dcminfotmp.SliceLocation;
    catch
        continue
    end
end
```

#### Writing a series of DICOM files to disk

In order to get a unique new image series we want to keep the study instance UID to have our new series be assigned to the existing study. But we want to create a new series instance UID to avoid creating a copy of our existing series. We also need to create new SOP instance UIDs, which are the identifiers that make each DICOM file unique. We will leave this to the dicomwrite.

```
%% sort the slices by instance number to get them into the righ order
[dummy sortindx] = sort(instancenumbervec);
file_list_dicom = file_list_dicom(sortindx);

%% now save that file as a series of DICOM images again
SeriesInstanceUID = dicomuid;

for i =1:size(data.img,3),
  metadata = dicominfo(file_list_dicom{i});
  metadata.SeriesInstanceUID = SeriesInstanceUID;
  SeriesDescription = sprintf('%s (nifti to dcm)', metadata.SeriesDescription);
  metadata.SeriesDescription = SeriesDescription;

  % need to adjust the voxel size
  metadata.PixelSpacing(1,1) = spacing(1);
  metadata.PixelSpacing(2,1) = spacing(2);
  metadata.SliceLocation = (i-1)*slice_location_delta;
  %%fprintf('%f pixel spacing in DICOM', metadata.PixelSpacing);
  metadata.ImagePositionPatient = [metadata.SliceLocation, 0, 0];
  metadata.InstanceNumber = i;
  metadata.ImageOrientationPatient = [0, 1, 0, 0, 0, -1]';
  
  fname_out = sprintf('%s/im%3.4i.dcm',output, i);
  fprintf('%s - Writing DICOM images, %d of %d\r',fname_out,i,size(data.img,3));
  dicomwrite(int16(rot90(data.img(:,:,i))), fname_out, metadata, 'CreateMode', 'Copy');
end
```

Writing out color encoded DICOM images (cannot be read by all workstations, OsiriX/Horos is fine).
```
outputdir = sprintf('%s',outdir);
% fname = sprintf('%s/volDV_RGB.mgz',outdir);
dicomdir = dicomsT0;
tag = 'MMIV T1 + DV';
SeriesNumber = 2000;
Mreg = M_hatl_to_MPR;

SeriesInstanceUID = dicom_generate_uid('instance');
pause(1);
SeriesDescription = sprintf('%s',tag);


dat = uint8(255*vol.imgs);
% try to fix front-back issue
dat = flipdim(dat,3);
T1_1_atl.imgs = flipdim(T1_1_atl.imgs, 3);
T1_2_atl.imgs = flipdim(T1_2_atl.imgs, 3);

dout = sprintf('%s/DV_RGB',outputdir);
mkdir(dout);
for i = 1:vol.dimd
    fprintf('%s - Writing DICOM images, %d of %d\r',mfilename,i,vol.dimd);
    metadata = dcminfotmp;
    metadata.PixelSpacing = [vol.vx vol.vy]';
    metadata.Columns = vol.dimc;
    LPHcoord = Mvx2lph*[0 0 i-1 1]';% Need to check this...
    metadata.SliceLocation = LPHcoord(3);
    metadata.ImagePositionPatient = LPHcoord(1:3)';
    metadata.ImageOrientationPatient = [Mvx2lph(1:3,1);Mvx2lph(1:3,2)]';
    metadata.Rows = vol.dimr;
    metadata.AcquisitionMatrix = [vol.dimc vol.dimr 0 0]';
    metadata.ColorType = 'truecolor';
    metadata.SamplesPerPixel = 3;
    metadata.InstanceNumber = i;
    metadata.PhotometricInterpretation = 'RGB';
    metadata.SmallestImagePixelValue = 0;
    metadata.LargestImagePixelValue = 255;
    metadata.BitsAllocated = 8;
    metadata.BitsStored = 8;
    metadata.HighBit = 7;
    metadata.BitDepth = 8;
    metadata.SeriesDescription = SeriesDescription;
    metadata.SeriesNumber = SeriesNumber;
    metadata.SeriesInstanceUID = SeriesInstanceUID;
    metadata.SoftwareVersions = version;
    fname_out = sprintf('%s/rgb_im%3.4i.dcm',dout,i);
    %slice = permute(squeeze(dat(:,:,vol.dimd-i+1,:)),[2 1 3]);
    slice = permute(squeeze(dat(:,:,i,:)),[2 1 3]);
    slice = flipdim(slice,2);
    dicomwrite(slice, fname_out, metadata, 'CreateMode', 'copy');
end
fprintf('\n%s - Finished T1 + DV\n',mfilename);
```

Here is a bit that will calculate a good initial window level for the data. This is probably one of the most useful pieces of code for medical imaging. The idea is to use the proportion of body against background as a means to calculate two thresholds for setting the images brightness and contrast. This information can be added to the DICOM tags and DICOM aware display programs will be able to use this information for the initial display. If you write novel image information you should calculate this information.

The calculation is parametrized by two thresholds as the proportion of low intensity pixel and the proportion of high intensity pixel. Usual values are 0.01 and 0.999. This varies by the body part and modality displayed. It is appropriate if the object of interest has intermediate intensities (not too bright such as in metal or contrast in vessels, and not too dark such as in the background).
```
[hc hv] = hist(T1_1_atl.imgs(find(T1_1_atl.imgs>0)),1000);
cdf_hc = cumsum(hc/sum(hc));
val0 = hv(min(find(cdf_hc>.99)));
im0 = max(0,min(1,T1_1_atl.imgs/val0));
```

In order to show a map (activity, segmented region etc.) overlayed on structural data you can use a simple alpha blending procedure. A linear weighting of the overlap map (in color) onto the gray-scale structural underlay map. The result of the procedure is a color image (RGB) that can be written as a DICOM using the above procedure.
```
% weighting f between overlay and underlay 
f = .6;
volR(inds)  = colors(round(aseg_atl.imgs(inds))+1,1);
a = max(0, min(1, real(smooth3d(volR, 5,5,5))));
volASEG_RGB(:,:,:,1) = f*a + (1-f)*im1(:,:,:,1);

volG(inds) = colors(round(aseg_atl.imgs(inds))+1,2);
a = max(0, min(1, real(smooth3d(volG, 5,5,5))));
volASEG_RGB(:,:,:,2) = f*a + (1-f)*im1(:,:,:,2);

volB(inds) = colors(round(aseg_atl.imgs(inds))+1,3);
a = max(0, min(1, real(smooth3d(volB, 5,5,5))));
volASEG_RGB(:,:,:,3) = f*a + (1-f)*im1(:,:,:,3);
```

## Detecting the type of a scan

Image sequences produce different contrast images. Some processing pipelines only work with a specific contrast distribution per tissue type and body region, and therefore DICOM series have to be filtered by type before processing. The problematic bit is that the sequence depends on a number of parameters and that therefore a concise description of the scanning sequences is difficult. Here an example of a dumb-method for detecting the type of the scan.

```
#!/bin/bash
# enter each directory (assume one series per directory and classify the data)

# limit output to one type by calling
# FILTER=T1 ./classifyData.sh
# FILTER=T2 ./classifyData.sh
# FILTER=PD ./classifyData.sh

## Detect the type of a DICOM image
#
#### Rule for T1, T2, and PD
#
#TE < TR
#short TR: usually lower than 500ms, equal to the average T1
#long TR: 2 times the short TR, usually greater than 1500ms
#short TE: is usually lower than 30ms
#long TE: is 3 times the short TE usually greater than 90ms
#
#T1: TR short 300-600ms
#T1: TE short 10-30ms
#T2: TR long: 2000ms
#T2: TE long: 90-140ms
#PD: TR long: 1000-3000ms
#PD: TE short: 15ms
#
#Proton density-weighted: long TR and short TE
#T1-weighted: short TR and short TE
#T2-weighted: long TR and long TE
#
# This might be outside the scope: https://radiopaedia.org/articles/mri-sequence-parameters
#

dirs=()
numFiles=()
input=$1

# filter for either T1, T2, or PD, or ALL
t=${FILTER}
if [ -z "$t" ]; then
    t="ALL"
fi

if [ -z "${input}" ]; then
    while IFS= read -r -d $'\0'; do
	# how many files in this directory?
	nfiles=$(find $REPLY -mindepth 1 -maxdepth 1 -type f -print | wc -l)
	if [ ${nfiles} -ge 10 ]; then
	    #echo "Classify this directory ($REPLY).. ${nfiles} files inside"
	    dirs+=("$REPLY")
	    numFiles+=($nfiles)
	#else
	    #echo "Ignore this directory ($REPLY).. less than 3 files ${nfiles}  inside"
	fi
    done < <(find . -type d -print0)
else
    nfiles=$(find "${input}" -mindepth 1 -maxdepth 1 -type f -print | wc -l)    
    dirs+=("${input}")
    numFiles+=($nfiles)
fi

classify () {
    local run=$1
    local numFiles=$2
    local t=$3

    # get the TE (echo time, 0018,0081) and TR (repetition time, 0018, 0080)
    f=$(ls ${run}/* |head -1)
    #echo "test file : \"$f\""
    TI=$(dcmdump +P "0018,0082" "$f" | head -1 | cut -d"[" -f2 | cut -d"]" -f1 | head -1 | cut -d'.' -f1)
    TE=$(dcmdump +P "0018,0081" "$f" | head -1 | cut -d"[" -f2 | cut -d"]" -f1 | head -1 | cut -d'.' -f1)
    TR=$(dcmdump +P "0018,0080" "$f" | head -1 | cut -d"[" -f2 | cut -d"]" -f1 | head -1 | cut -d'.' -f1)
    CO=$(dcmdump +P "0018,1048" "$f" | head -1 | cut -d"[" -f2 | cut -d"]" -f1 | head -1 | cut -d'.' -f1)
    SEQUENCENAME=$(dcmdump +P "0018,0024" "$f" | head -1 | cut -d"[" -f2 | cut -d"]" -f1 | head -1 | cut -d'.' -f1)
    if [[ $SEQUENCENAME == *"no value available"* ]]; then
	SEQUENCENAME="no sequence name"
    fi
    if [[ "$TI" == *"no value available"* ]]; then
	TI="no value"
    fi
    if [[ "$CO" == *"no value available"* ]]; then
	CO="no contrast"
    fi    
    if [ -z "$CO" ]; then
	CO="no contrast"
    fi
    if [ -z "$TI" ]; then
	TI="no value"
    fi
    #echo "TE: $TE, TR: $TR"
    if [ -z "$TR" ]; then
	TR=-1
	(>&2 echo "no TR in $f")
    fi
    if [ -z "$TE" ]; then
	TE=-1
	(>&2 echo "no TE in $f")
    fi
    
    # small large TE, TR based on tissue properties in brain
    TRshort=0
    TRlong=0
    TEshort=0
    TElong=0
    if [ ! "$TR" -eq "-1" -a "$TR" -le 600 ]; then
	TRshort=1
	#echo "TR short"
    elif [ "$TR" -ge 600 ]; then
	TRlong=1
	#echo "TR long"
    fi
    if [ ! "$TE" -eq "-1" -a "$TE" -le 30 ]; then
	TEshort=1
	#echo "TE short"
    elif [ "$TE" -ge 30 ]; then
	TElong=1
	#echo "TE long"
    else
	(>&2 echo "TE strange?")
    fi
    # weighted scan
    T1=0
    T2=0
    PD=0
    FL=0
    if [ "${TEshort}" -eq 1 -a "${TRshort}" -eq 1 ]; then
	T1=1
	if [ $t == "T1" -o $t == "ALL" ]; then
	    echo -e "$run"
	    (>&2 echo -e "T1 weighted image series [#$numFiles, $SEQUENCENAME] (TR: $TR - short, TE: $TE - short, TI: $TI, contrast: $CO)")
	fi
    elif [ "${TElong}" -eq 1 -a "${TRlong}" -eq 1 ]; then
	# flip angle should be 90 (inversion time between 1700 and 2200 is FLAIR)
	T2=1
	if [ $t == "T2" -o $t == "ALL" ]; then
	    echo -e "$run"
	    (>&2 echo -e "T2 weighted image series [#$numFiles, $SEQUENCENAME] (TR: $TR - long, TE: $TE - long, TI: $TI, contrast: $CO)")
	fi
    elif [ "${TRlong}" -eq 1 -a "${TEshort}" -eq 1 ]; then
	PD=1
	if [ $t == "PD" -o $t ==  "ALL" ]; then	
	    echo -e "$run"
	    (>&2 echo -e "PD weigthed image series [#$numFiles, $SEQUENCENAME] (TR: $TR - long, TE: $TE - short, TI: $TI, contrast: $CO)")
	fi
    fi    
}


# iterate over the array and export 4 folders at the same time
i=1
N=8
(
for ((j = 0; j < ${#dirs[@]}; j++))
do
    ((i=i%N)); ((i++==0)) && wait
    #echo "a directory \"${dirs[$j]}\""
    classify "${dirs[$j]}" "${numFiles[$j]}" "$t" &
done
)		   

```
