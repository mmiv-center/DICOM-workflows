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
tag = 'MMIL T1 + DV';
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

Here is a bit that will calculate a good initial window level for the data
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
