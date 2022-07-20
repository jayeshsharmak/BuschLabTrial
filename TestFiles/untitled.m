%% Enter per experiment spectral data
par.currentDIR=pwd;
par.nigrosinFile=uigetfile([par.datDir '/*'],'C:\Users\rekha\Downloads\20220708_690nm.txt');


%% read in spectral data
% read function for new Ocean Optics software
data.tmpNig=readOceanOpticsSpectroV2([par.datDir par.nigrosinFile]);
data.tmpFhi=readOceanOpticsSpectroV2([par.datDir par.fhiFile]);
data.tmpWat=readOceanOpticsSpectroV2([par.datDir par.waterFile]);
% clip negative values, typically in UV
data.tmpNig(data.tmpNig(:,2)<0,2)=NaN;
data.tmpFhi(data.tmpFhi(:,2)<0,2)=NaN;
data.tmpWat(data.tmpWat(:,2)<0,2)=NaN;