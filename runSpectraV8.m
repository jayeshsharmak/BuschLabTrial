% code rewrite 20220505 drb
% discrepencies in how spectra were read in- multiple hemoglobin/water
% spectra.

'note that this script does not provide perfect alignment of the input saturations and expected output values'

clear all;
path(path,'SubFunctions/')
% these are per-measurement editable parameters
par.isAbs=0; % is this absorption data or transmission?
par.debugPlot=1;
%% rarely changed params
par.issLbsFile='Parameters/issWavelengths.csv';
par.lbLim=5; % how far, in nm, from spectra ext. coef are allowed.
% scattering prefactors
par.A=40000;
par.b=1.2;

%par.C= par.A*par.b > par.b/ par.A 

% MUST match phantom/tissue type.
par.watFbio=0.756;
par.fatFbio=0.117;
par.watF=0.97832; %0.756;
par.fatF=0.01828; %0.117;

%% reference spectra
% read in spectra
% % hemoglobin spectra from prahl
% par.wholeBloodFactor=NaN- we aren't useing whole blood and measuring in hb in grams 150/64500;
par.fileHb='Spectra/hemoglobin-Prahl-baseE.dat';
par.fileWat='Spectra/segelstein-Water-AbsorptionInvCMFull.dat';
par.fileFat='Spectra/veen-fat-absInvINMETERSfull.dat';

% hemoglobin
data.tmpHbarr=dlmread(par.fileHb);
data.spec.HbO=data.tmpHbarr(:,1:2);
data.spec.Hb =data.tmpHbarr(:,[1 3]);

% water
data.spec.Wat=dlmread(par.fileWat);
% datWat(:,3)=datWat(:,2)*par.watF; %phantom
% datWat(:,4)=datWat(:,2)*par.watFbio; %bio

% fat
data.spec.Fat=dlmread(par.fileFat);
data.spec.Fat(:,2)=data.spec.Fat(:,2)/100; % spectra in meters
% datFat(:,3)=datFat(:,2)*par.fatF; % phantom
% datFat(:,4)=datFat(:,2)*par.fatFbio; % bio

%% Enter per experiment spectral data
par.currentDIR=pwd;
if ispc
    par.datDir=[uigetdir(par.currentDIR,'Select data directory') '\'];
    % local data directory- this is where the spectra files are placed
    % datDir='Data/Phantoms/20220504_dyespectra/';
else
    par.datDir=[uigetdir(par.currentDIR,'Select data directory') '/'];
end
par.nigrosinFile=uigetfile([par.datDir '/*'],'Select nigrosin file');
par.dilute.Nig=input('enter nigrosin dilution factor\n');

par.fhiFile=uigetfile([par.datDir '/*'],'Select FHI file');
par.dilute.Fhi=input('enter FHI dilution factor\n');

par.waterFile=uigetfile([par.datDir '/*'],'Select water/buffer file');
%% simiple params
% goal total hemoglobin
par.hbt=input('enter target total hemoglobin concentration in uM (e.g., 37)\n');
par.hbt=par.hbt*1e-6;
'converting hbt to M from uM'
% goal saturation
par.sat=input('enter target saturation in percent, bracketed (e.g., [70 80])\n');
par.sat=par.sat/100;% convert to sensible fractions
% sources
par.sourcesused=input('enter ISS sources used, bracketed (e.g., [1 2 3])\n');

% are the spectra transmission or absorption?
par.isAbs=input('enter 1 if spectra are absorption, 0 if transmission\n');
%
par.issLbs=csvread(par.issLbsFile,1,0); % skip first line
'using iss wavelengths'
par.lbs=unique(par.issLbs(:,2));par.lbs=sort(par.lbs);
% saving
par.figSave=0; %auto save figures?
par.timeVec=timestring;

% plotting
par.plotLimLb=[650 850];
par.cmap=jet(length(par.sourcesused));
% wavelength limits for plot
par.xl=[350 900];
par.nirWind=[625 875];
par.visWind=[400 700];

% fitting
% dyes, INPUT dilution fraction starting elements for fit.
% NOT measured dilution fractions for spectra
par.fhiF=0;
par.nigF=5e-5;


par.op=optimset('fminsearch');
par.options=optimset(par.op,'MaxIter',1e4,'MaxFunEvals',1e4,...
    'TolFun',1.000e-16,'TolX',1.000e-16,'Display','none');


%% derived variables
% identify output name- name of directory spectra are stored in.
if ispc==0
tind=strfind(par.datDir,'/');
elseif ispc==1
    tind=strfind(par.datDir,'\');
else
    error('no macs allowed')
end
par.outNamePre=[par.datDir par.datDir((tind(end-1)+1):(tind(end)-1))];
% par.outNamePre
clear tind
%% read in spectral data
% read function for new Ocean Optics software
data.tmpNig=readOceanOpticsSpectroV2([par.datDir par.nigrosinFile]);
data.tmpFhi=readOceanOpticsSpectroV2([par.datDir par.fhiFile]);
data.tmpWat=readOceanOpticsSpectroV2([par.datDir par.waterFile]);
% clip negative values, typically in UV
data.tmpNig(data.tmpNig(:,2)<0,2)=NaN;
data.tmpFhi(data.tmpFhi(:,2)<0,2)=NaN;
data.tmpWat(data.tmpWat(:,2)<0,2)=NaN;
% if spectral data are transmission, convert to absorption
% columns: wavelength nm, measured abs spectra diluted, measured abs
% spectra
if par.isAbs==0
    'transmission data'
    data.spec.Nig(:,1)=data.tmpNig(:,1);
    data.spec.Nig(:,3)=-log(data.tmpNig(:,2)./data.tmpWat(:,2)); % alreadyin base e
    data.spec.Nig(:,2)=data.spec.Nig(:,3)*par.dilute.Nig;
    
    data.spec.Fhi(:,1)=data.tmpFhi(:,1);
    data.spec.Fhi(:,3)=-log(data.tmpFhi(:,2)./data.tmpWat(:,2)); % alreadyin base e
    data.spec.Fhi(:,2)=data.spec.Fhi(:,3)*par.dilute.Fhi;
elseif param.isAbs==1
    error('20220505 have not tested absorption data')
    data.spec.Nig(:,1)=data.tmpNig(:,1);
    data.spec.Nig(:,3)=data.tmpNig(:,2);
    data.spec.Nig(:,2)=data.spec.Nig(:,3)*par.dilute.Nig;
    
    data.spec.Fhi(:,1)=data.tmpFhi(:,1);
    data.spec.Fhi(:,3)=data.tmpFhi(:,2);
    data.spec.Fhi(:,2)=data.spec.Fhi(:,3)*par.dilute.Fhi;
    
else
    param.isAbs
    error('haven''t dealt with non-absorption or transmission spectra')
end

%% test plot
if par.debugPlot==1
    figure(33);clf
    subplot(2,2,1)
    if par.isAbs==0
        plot(data.tmpWat(:,1),data.tmpWat(:,2),'c')
        hold on
        plot(data.tmpNig(:,1),data.tmpNig(:,2),'r')
        plot(data.tmpFhi(:,1),data.tmpFhi(:,2),'b')
        xlim([par.plotLimLb])
        xlabel('\lambda [nm]')
        ylabel('Transmission')
    end
    subplot(2,2,2)
    plot(data.spec.Nig(:,1),data.spec.Nig(:,2),'r')
    hold on
    plot(data.spec.Fhi(:,1),data.spec.Fhi(:,2),'b')
    xlim([par.plotLimLb])
    xlabel('\lambda [nm]')
    ylabel('Absorbance')
    
    subplot(2,2,3)
    plot(data.spec.Nig(:,1),data.spec.Nig(:,3),'r')
    hold on
    plot(data.spec.Fhi(:,1),data.spec.Fhi(:,3),'b')
    xlim([par.plotLimLb])
    xlabel('\lambda [nm]')
    ylabel('\mu_a soln')
    
    subplot(2,2,4)
    plot(data.spec.Nig(:,1),data.spec.Nig(:,2),'r')
    hold on
    plot(data.spec.Fhi(:,1),data.spec.Fhi(:,2),'b')
    xlim([par.plotLimLb])
    xlabel('\lambda [nm]')
    ylabel('\mu_a stock')
end

%% identify extinction coefficients at the wavelengths of interest
tnames=fieldnames(data.spec);
for ii=1:length(tnames)
    tstr=['tdat=data.spec.' tnames{ii} ';'];
    eval(tstr)
    for ilb=1:length(par.lbs)
        [tval, tind]=min(abs(tdat(:,1)-par.lbs(ilb)));
        tstr2=['data.specL. ' tnames{ii} '(ilb,:)=[' num2str(par.lbs(ilb)) ' ' num2str(tdat(tind,1)) ' ' num2str(tdat(tind,2)) '];'];
        eval(tstr2)
    end;clear ilb
    
end;clear ii tstr tnames tdat tstr2 tind tval
%% target spectra, from biology
data.muaBio=[];
for iSat=1:length(par.sat)
    ['targetting Hbt=' num2str(par.hbt) ' StO2=' num2str(par.sat(iSat))]
    tmat=[data.specL.HbO(:,3)*par.hbt*par.sat(iSat) ...
        data.specL.Hb(:,3)*par.hbt*(1-par.sat(iSat)) ...
        data.specL.Wat(:,3)*par.watFbio...
        data.specL.Fat(:,3)*par.fatFbio ];
    data.muaBio(iSat,:)=sum(tmat,2);
    clear tmat
end;clear iSat

%% water/fat spectra from phantom
data.fatPhant=data.specL.Fat(:,3)*par.fatF;
data.watPhant=data.specL.Wat(:,3)*par.watF;

%% calculations
% muaBio has the water and fat fractions for tissue embedded in it.
% here, we want to acknowledge the absorption of water (fat) in the
% phantom and NOT add dyes to compensate for this absorption.
data.calcFracArr=[];
for iSat=1:length(par.sat)
    [calcFrac,fval,exitflag]=...
        fminsearchbnd(@fitDyeSpectralDiffV2,[par.nigF par.fhiF],...
        [0 0],[1 1],par.options,[data.specL.Nig(:,3) data.specL.Fhi(:,3)],...
        data.muaBio(iSat,:)'-data.watPhant-data.fatPhant);
    data.calcFracArr(iSat,:)=calcFrac;
end;clear iSat fval exitflag calcFrac

% calculate the mua for the phantom at the relevent wavelengths
data.muaPhant=[];
for iSat=1:length(par.sat)
    data.muaPhant(iSat,:)=data.calcFracArr(iSat,1)*data.specL.Nig(:,3)+...
        data.calcFracArr(iSat,2)*data.specL.Fhi(:,3)+...
        data.watPhant+data.fatPhant;
end;clear iSat
data.muaPhant
data.muaBio

%% sanity check- fit the planned phantom mua
data.expectedChromoOut=[];
for iSat=1:length(par.sat)
    
    muaIn=data.muaPhant(iSat,:)'-data.watPhant-data.fatPhant
%     sigmaS=fitHbOHbFixedWatFatV1(concHbs,muaIn,muaHbO,muaHb)
    
    [concHbs,fval,exitflag]=...
        fminsearchbnd(@fitHbOHbFixedWatFatV1,[1 1]*10e-6,...
        [0 0],[1 1],par.options,muaIn,...
        data.specL.HbO(:,3),...
    data.specL.Hb(:,3))
    
data.HbConcentrationsDye(iSat,:)=concHbs;
['HbO ' num2str(par.hbt*par.sat(iSat)*1e6) '\mu{M}']
        ['Hb ' num2str(par.hbt*(1-par.sat(iSat))*1e6) '\mu{M}']
data.expectedChromoOut(iSat,:)=[concHbs*1e6 sum(concHbs)*1e6 concHbs(1)/sum(concHbs)];    
end;clear iSat concHbs fval exitflag muaIn
data.HbConcentrationsDye*1e6

{'HbO','Hb','Hbt','StO2','Nominal StO2'}
([data.expectedChromoOut par.sat'])

%% make output table
data.out=[];
data.out.satTarg=par.sat';
data.out.satExp=round(data.expectedChromoOut(:,4),3);
data.out.hbtExp=round(data.expectedChromoOut(:,3),3);
data.out.hbExp=round(data.expectedChromoOut(:,2),3);
data.out.hboExp=round(data.expectedChromoOut(:,1),3);
data.out.nigVolFrac=round(data.calcFracArr(:,1),6);
data.out.fhiVolFrac=round(data.calcFracArr(:,2),6);
data.outTab=struct2table(data.out);

writetable(data.outTab,[par.outNamePre '-table.csv'])
%% plot 
fid=figure(35);clf
[h, w]=hwcalc(length(par.sat));
for iSat=1:length(par.sat)
    subplot(h,w,iSat)
    %     setAxisProp(36,6)
    lWidth=4;
    % plot the biology
    plot(data.spec.HbO(:,1),data.spec.HbO(:,2)*par.hbt*(par.sat(iSat)),'r','LineWidth',lWidth); hold on
    plot(data.spec.Hb(:,1),data.spec.Hb(:,2)*par.hbt*(1-par.sat(iSat)),'b','LineWidth',lWidth)
    plot(data.spec.Wat(:,1),data.spec.Wat(:,2)*par.watFbio,'c','LineWidth',lWidth)
    plot(data.spec.Fat(:,1),data.spec.Fat(:,2)*par.fatFbio,'g','LineWidth',lWidth)
    plot(par.lbs,data.muaBio(iSat,:),'o','MarkerSize',40,'Color',[0.8 0.8 0.8],'MarkerFaceColor',0.4* [1 1 1])
    
    % plot the phantom
    
    plot(data.spec.Fhi(:,1),data.spec.Fhi(:,2)*data.calcFracArr(iSat,2),'--r','LineWidth',lWidth); hold on
    plot(data.spec.Nig(:,1),data.spec.Nig(:,2)*data.calcFracArr(iSat,1),'--b','LineWidth',lWidth)
    plot(data.spec.Wat(:,1),data.spec.Wat(:,2)*par.watF,'--c','LineWidth',lWidth)
    plot(data.spec.Fat(:,1),data.spec.Fat(:,2)*par.fatF,'--g','LineWidth',lWidth)
    
    
    plot(par.lbs,data.muaPhant(iSat,:),'d','MarkerSize',40,'Color','k','MarkerFaceColor','k')
    %     ,[0.8 0.8 0.8],'MarkerFaceColor',0.95* [1 1 1])
    xlim([par.plotLimLb])
    xlabel('\lambda [nm]')
    ylabel('\mu_a/cm')
    
    legend({['HbO ' num2str(par.hbt*par.sat(iSat)*1e6) '\mu{M}'];...
        ['Hb ' num2str(par.hbt*(1-par.sat(iSat))*1e6) '\mu{M}'];...
        ['H2O (' num2str(par.watFbio,3) ')'];...
        ['Fat (' num2str(par.fatFbio,3) ')'];...
        'Total';...
                'FHI';...
                'Nigrosin';...
        ['H2O phant (' num2str(par.watF,3) ')'];...
        ['Fat phant (' num2str(par.fatF,3) ')'];...
        'Total phant'},'Location','EastOutside')
    
    
    set(gca,'LineWidth',lWidth,'FontSize',18)
    xlabel('Wavelength [nm]')
    ylabel('\mu_a/cm')
    ylim([0 0.15])
    %    vh=vline(par.lbs);set(vh,'LineWidth',2*lWidth);
    title({['target sat ' num2str(par.sat(iSat)*100) '%, expected ' num2str(data.expectedChromoOut(iSat,4)*100,3) '%'];...
        ['fNig.=' num2str(data.calcFracArr(iSat,1),3) ', fFHI=' num2str(data.calcFracArr(iSat,2),3)]})
    
end;clear iSat lWidth h w fid


%% 
% 'this needs files written out'
save([par.outNamePre '.mat'])

