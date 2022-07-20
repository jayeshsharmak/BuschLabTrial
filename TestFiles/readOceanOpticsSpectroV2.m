% 180327 D Busch
% read in spectar from ocean optics spectrophotometers
% data (wavelength, abs) begins on line after:
% >>>>>Begin Spectral Data<<<<<

function datOut=readOceanOpticsSpectro(fileIn);
goOn=1;
    fid=fopen(fileIn);
    tind=1;
    while goOn>=1;

    tline=fgetl(fid);
    tbool=strcmp(tline,'>>>>>Begin Spectral Data<<<<<');
    if tbool==0
        tind=tind+1;
    elseif tbool==1
        goOn=0;
    end
end
fclose(fid);

datOut=dlmread(fileIn,'\t',tind,0);
% datOut=reshape(datOut,2,[])';