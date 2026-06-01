function ImageProcessingGUI

% =========================================================================
% ImageProcessingGUI
% =========================================================================
%
% Author:
%   Edcer Laguda
%
% Affiliation:
%   Department of Physics and Astronomy
%   McMaster University
%   Hamilton, Ontario, Canada
%
% Developed:
%   2025
%
% Description:
%   ImageProcessingGUI is a MATLAB-based graphical user interface developed
%   for quantitative and qualitative neutron imaging using a CZT-based gamma
%   camera system. The GUI supports energy-resolved neutron imaging workflows,
%   including detector characterization, transmission imaging, attenuation
%   mapping, uniformity correction, detector efficiency estimation, and
%   macroscopic cross-section analysis.
%
% Main Features:
%   - Open-beam and object data loading
%   - Energy-corrected and non-energy-corrected spectrum visualization
%   - Flat-field correction map generation
%   - Uniformity correction
%   - Transmission image generation
%   - Attenuation mapping using A = -ln(I/I0)
%   - Macroscopic cross-section estimation using Sigma = -ln(I/I0)/d
%   - ROI-based quantitative image analysis
%   - Windowed energy image generation
%   - Detector efficiency estimation
%   - Energy-dependent count-rate analysis
%   - Reference model comparison
%   - Interactive image inspection and spectral analysis
%
% MATLAB Requirements:
%   Required:
%     - MATLAB R2022b or newer recommended
%     - Image Processing Toolbox
%
%   Recommended:
%     - Statistics and Machine Learning Toolbox
%
% Expected Input Variables:
%   The input .mat files should contain:
%     - LMmap_Ecorrect
%     - LMmap_noEcorrect
%     - keV, Energy_keV, E_keV, energy_keV, Eaxis_keV, or E
%
% Detector Configuration:
%   - CZT detector matrix: 32 x 32 pixels
%   - Pixel pitch: 2.46 mm
%   - Active field of view: approximately 7.9 cm x 7.9 cm
%   - Circular neutron beam diameter: 7.3 cm
%
% Quantitative Assumptions:
%   - Open-beam and object acquisitions use the same exposure time.
%   - Transmission is calculated as:
%
%         T = I/I0
%
%   - Attenuation is calculated as:
%
%         A = -ln(T)
%
%   - Macroscopic cross section is calculated as:
%
%         Sigma = -ln(T)/d
%
%     where d is the object thickness in cm.
%
% Research Application:
%   This software was developed as part of doctoral research on neutron
%   imaging, neutron attenuation measurements, detector characterization,
%   and quantitative neutron radiography using a CZT-based gamma camera
%   system.
%
% Disclaimer:
%   This software is provided for research and educational use. Users are
%   responsible for validating all calculations, assumptions, and results
%   before use in publications, clinical work, engineering design, or
%   regulatory applications.
%
% License:
%   Copyright (c) 2025 Edcer Laguda
%
%   Permission is hereby granted, free of charge, to any person obtaining a
%   copy of this software and associated documentation files to use, copy,
%   modify, merge, publish, distribute, sublicense, and/or sell copies of
%   the software, subject to inclusion of this copyright notice in all copies
%   or substantial portions of the software.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
%
% =========================================================================


%% === CONSTANTS / STATE ===
NX = 32; NY = 32;
PIX_CM = 0.246;
BEAM_DIAM_CM = 7.3;
CENTER = (NX+1)/2;
R_PIX = (BEAM_DIAM_CM/2)/PIX_CM;
NA = 6.02214076e23;

COL_NOEC = [0 0 1];
COL_EC   = [1 0 0];

data = struct();
data.beamMode = 'LEM'; 
data.NX = NX; data.NY = NY; data.PIX_CM = PIX_CM;
[data.X,data.Y] = meshgrid(1:NX,1:NY);
data.mask = (data.X - CENTER).^2 + (data.Y - CENTER).^2 <= R_PIX^2;
data.roi_idx = builtin('find', data.mask(:));
data.beam_area_cm2 = pi*(BEAM_DIAM_CM/2)^2;
data.lastOpenPath = pwd;
data.lastObjPath  = pwd;
data.lastCsvPath  = pwd;
data.filepath     = pwd;
data.prefix       = 'Object';
data.results = struct('has_roi',false,'sigma_exp_cm1',NaN,'t_mean',NaN,'t_pred',NaN, ...
                      'sigma_ref_cm1',NaN,'t_ref',NaN,'note_ref','');
data.thickness_cm = 1.0;

% Reference constants
refSigmaThermal_b = struct('H_nat',0.3326,'C_nat',0.0035,'O_nat',0.00019,'Ni_nat',4.49, ...
    'B_nat',767,'B10',3835,'Cd_nat',2520,'Cd113',20600,'Gd_nat',49700,'Zn_nat',1.11,'Te_nat',4.3, ...
    'Li_nat',70.5,'Mn_nat',13.3);   % <-- added Li, Mn

refProps = struct('H_nat',struct('A',1.008,'rho',1.0), ...
    'C_nat',struct('A',12.011,'rho',2.267), ...
    'O_nat',struct('A',15.999,'rho',1.0), ...
    'Ni_nat',struct('A',58.6934,'rho',8.90), ...
    'B_nat',struct('A',10.81,'rho',2.34), ...
    'B10',struct('A',10.0129,'rho',2.34), ...
    'Cd_nat',struct('A',112.414,'rho',8.65), ...
    'Cd113',struct('A',112.904,'rho',8.65), ...
    'Gd_nat',struct('A',157.25,'rho',7.90), ...
    'Zn_nat',struct('A',65.38,'rho',7.14), ...
    'Te_nat',struct('A',127.60,'rho',6.24), ...
    'Li_nat',struct('A',6.94,'rho',0.534), ...   
    'Mn_nat',struct('A',54.938,'rho',7.21));     

data.endf = struct(); 
data.spectrum = [];

%% === UI ROOT ===
fig = uifigure('Name','Image Processing GUI','Position',[60 40 1400 920]);
fig.AutoResizeChildren = 'on';    
fig.SizeChangedFcn     = [];      


% Tabs
tabgp = uitabgroup(fig,'Position',[10 10 fig.Position(3)-20 fig.Position(4)-20]); 
tabgp = uitabgroup(fig,'Position',[10 10 1380 900]); 
tabLoad = uitab(tabgp,'Title','File Loader'); 
tabProcess = uitab(tabgp,'Title','Correction Map Generator'); 
tabObjectSpectra = uitab(tabgp,'Title','Object Energy Spectra'); 
tabUniformity = uitab(tabgp,'Title','Uniformity & Attenuation'); 
tabInspect = uitab(tabgp,'Title','Inspection');
tabWindowedImages = uitab(tabgp,'Title','Windowed Energy Images'); 
tabWindowRate = uitab(tabgp,'Title','Windowed Count-Rate Plot');
tabTrans           = uitab(tabgp,'Title','Transmission & Sigma (ROI)');      
tabMacroXS         = uitab(tabgp,'Title','Reference Model & Comparison');      

%% ------------------ Tab 1: File Loader ------------------
uibutton(tabLoad,'Text','Load Open-Beam .mat', ...
    'Position',[30 640 250 30], 'ButtonPushedFcn', @(~,~)loadData());
loadStatus = uilabel(tabLoad,'Text','Status: No open-beam file loaded','Position',[30 610 1100 22]);

uibutton(tabLoad,'Text','Load Object .mat', ...
    'Position',[30 560 250 30], 'ButtonPushedFcn', @(~,~)loadObject());
objStatus = uilabel(tabLoad,'Text','Object Status: No file loaded','Position',[30 530 1100 22]);

uilabel(tabLoad,'Position',[30 480 300 22],'Text','Detected keV vector variable:');
keVLabel = uilabel(tabLoad,'Position',[30 455 1100 22],'Text','(none)');

uilabel(tabLoad,'Text','Prefix name:','Position',[30 720 100 22]);
prefixField1 = uieditfield(tabLoad,'text','Position',[140 720 240 24], ...
    'Value',data.prefix,'ValueChangedFcn',@(e,~)setPrefix(e.Value));

uilabel(tabLoad,'Text','Detector mode:','Position',[420 720 100 22]);
ddBeamMode = uidropdown(tabLoad,'Items',{'LEM','HEM'}, ...
    'Position',[520 720 120 24], 'Value',data.beamMode, ...
    'ValueChangedFcn', @(~,~)setBeamMode());


%% ------------------ Tab 2: Correction Map Generator ------------------
xminField = uieditfield(tabProcess,'numeric','Position',[30 700 90 30],'Value',65);
xmaxField = uieditfield(tabProcess,'numeric','Position',[140 700 90 30],'Value',195);
uibutton(tabProcess,'Text','Set Energy Range (zoom)','Position',[30 660 200 30], ...
    'ButtonPushedFcn', @(~,~)updateZoom());

uilabel(tabProcess,'Text','Flat-field percentile (for correction map):', ...
    'Position',[30 620 280 22]);
ddPxx2 = uidropdown(tabProcess,'Items',{'P100','P95','P85'}, ...
    'Position',[30 595 120 24],'Value','P95');

uibutton(tabProcess,'Text','Generate & Save Correction Map','Position',[30 560 220 30], ...
    'ButtonPushedFcn', @(~,~)saveCorrectionMap());
uibutton(tabProcess,'Text','Export Map to PNG','Position',[30 520 220 30], ...
    'ButtonPushedFcn', @(~,~)exportCorrectionMapPNG());
saveStatus = uilabel(tabProcess,'Text','Output: Not saved','Position',[30 485 1100 22], 'WordWrap','on');

ax1 = uiaxes(tabProcess,'Position',[360 640 640 180]);
ax2 = uiaxes(tabProcess,'Position',[360 420 640 180]);
xlabel(ax1,'Energy (keV)'); ylabel(ax1,'Counts'); grid(ax1,'on');
xlabel(ax2,'Energy (keV)'); ylabel(ax2,'Counts'); grid(ax2,'on');
title(ax1,'Full Beam Spectrum (full energy range)');
title(ax2,'Full Beam Spectrum (selected energy range)');

axOB = uiaxes(tabProcess,'Position',[120 20 520 380]);   % original OB
ax3  = uiaxes(tabProcess,'Position',[720 20 520 380]);   % correction map

%% ------------------ Tab 3: Object Energy Spectra ------------------
axObjFull  = uiaxes(tabObjectSpectra,'Position',[300 640 640 180]);
axObjTrunc = uiaxes(tabObjectSpectra,'Position',[300 420 640 180]);
xlabel(axObjFull,'Energy (keV)'); ylabel(axObjFull,'Counts'); grid(axObjFull,'on');
xlabel(axObjTrunc,'Energy (keV)'); ylabel(axObjTrunc,'Counts'); grid(axObjTrunc,'on');

uilabel(tabObjectSpectra,'Position',[30 360 160 22],'Text','Flux (n/cm^2/s):');
fluxField = uieditfield(tabObjectSpectra,'numeric','Position',[190 360 120 22],'Value',1.2e4);
uilabel(tabObjectSpectra,'Position',[30 330 160 22],'Text','Exposure Time (s):');
timeField = uieditfield(tabObjectSpectra,'numeric','Position',[190 330 120 22],'Value',30);
uibutton(tabObjectSpectra,'Text','Compute ROI Intrinsic + Throughput','Position',[30 290 260 30], ...
    'ButtonPushedFcn', @(~,~)calcEfficiency());
effLabel = uilabel(tabObjectSpectra,'Position',[30 260 1100 22], 'Text','Efficiencies: not calculated');

%% ------------------ Tab 4: Uniformity & Attenuation ------------------
% Removed redundant prefix edit; single centered header reflecting Tab 1 prefix
labelBigUniform = uilabel(tabUniformity,'Text',data.prefix,'Position',[420 828 870 26], ...
    'FontSize',18,'FontWeight','bold','HorizontalAlignment','center');

uibutton(tabUniformity,'Text','Generate Images (no save)','Position',[30 812 240 30], ...
    'ButtonPushedFcn', @(~,~)generateUniformityImages());

uilabel(tabUniformity,'Text','Flat-field percentile:', ...
    'Position',[30 780 140 22],'HorizontalAlignment','left');
ddPercentile = uidropdown(tabUniformity,'Items',{'P100','P95','P85'}, ...
    'Position',[180 778 120 24],'Value','P100','ValueChangedFcn',@(~,~)onPercentileChanged());

uilabel(tabUniformity,'Text','Grayscale view:', ...
    'Position',[30 748 140 22],'HorizontalAlignment','left');
ddGrayMode = uidropdown(tabUniformity,'Items',{'Grayscale (Uniformity)','Grayscale of Amap','Transmission (T)'}, ...
    'Position',[180 746 220 24],'Value','Grayscale (Uniformity)','ValueChangedFcn',@(~,~)refreshGrayPreview());

uilabel(tabUniformity,'Text','Minimum Resolvable Feature Size: 4.92 mm', ...
    'Position',[30 716 420 22],'FontWeight','bold');

roiStatsPanel = uipanel(tabUniformity,'Title','ROI Stats','Position',[20 670 388 68]);
roiStatsLabel = uilabel(roiStatsPanel,'Position',[8 8 368 40], ...
    'Text','ROI stats: (populate after generation)','WordWrap','on');

roiPanel = uipanel(tabUniformity,'Title','ROI Tools','Position',[20 440 388 220]);
activeImgLabel = uilabel(roiPanel,'Text','Active image: (click a preview)','Position',[10 180 360 22]);
uibutton(roiPanel,'Text','Add ROI (rectangle)','Position',[10 146 180 26], ...
    'ButtonPushedFcn',@(~,~)addRoi('rect'));
uibutton(roiPanel,'Text','Add Ignore Region (rect -> 0)','Position',[200 146 180 26], ...
    'ButtonPushedFcn',@(~,~)ignoreRectToZero());
uibutton(roiPanel,'Text','Repopulate (pick pixel via small rect)','Position',[10 114 370 26], ...
    'ButtonPushedFcn',@(~,~)repopulateZerosFromPoint());

roiList = uilistbox(roiPanel,'Position',[10 10 180 92],'Items',{},'Multiselect','off');
uibutton(roiPanel,'Text','Delete ROI','Position',[200 76 180 26], ...
    'ButtonPushedFcn',@(~,~)deleteSelectedRoi());
uibutton(roiPanel,'Text','Clear All','Position',[200 44 180 26], ...
    'ButtonPushedFcn',@(~,~)clearAllRois());

roiOut = uitextarea(tabUniformity,'Position',[20 210 388 210],'Editable','off', ...
    'Value',{'ROI results will appear here.'});

exportPanel = uipanel(tabUniformity,'Title','Export','Position',[20 20 388 170]);
uilabel(exportPanel,'Text','Export mode:', 'Position',[12 100 100 22],'HorizontalAlignment','left');
ddExportMode = uidropdown(exportPanel,'Items',{'Separate PNGs','Single panel (2x2)'}, ...
    'Position',[120 98 200 24],'Value','Separate PNGs');
uibutton(exportPanel,'Text','Export PNGs','Position',[12 60 140 30], ...
    'ButtonPushedFcn',@(~,~)exportUniformityPNGs());
uibutton(exportPanel,'Text','Save Attenuation (.mat)','Position',[162 60 170 30], ...
    'ButtonPushedFcn',@(~,~)saveAmapMat());
uibutton(exportPanel,'Text','Save Uniformity (.mat)','Position',[12 20 170 30], ...
    'ButtonPushedFcn',@(~,~)saveUniformMat());

axOrig   = uiaxes(tabUniformity,'Position',[420 520 420 300]);
axUniform= uiaxes(tabUniformity,'Position',[870 520 420 300]);
axAtten  = uiaxes(tabUniformity,'Position',[420 120 420 300]);
axGray   = uiaxes(tabUniformity,'Position',[870 120 420 300]);

data.hImg = struct('orig',[],'uniform',[],'atten',[],'gray',[]);
data.lastActive = '';
data.rois = struct('name',{},'obj',{},'type',{},'on',{},'mask',{},'axesTag',{}); % uniformity tab ROIs
if ~isfield(data,'lastUniformity') || ~isstruct(data.lastUniformity), data.lastUniformity = struct(); end
if ~isfield(data.lastUniformity,'flat_percentile'), data.lastUniformity.flat_percentile = 95; end


SL_X_MIN = 1298; SL_X_MAX = 1336; SL_Y = 520; SL_H = 300; SL_W = 3;
data.winLimits = [0 1];
sMin = uislider(tabUniformity,'Orientation','vertical', ...
    'Position',[SL_X_MIN SL_Y SL_W SL_H],'Limits',[0 1],'Value',0, ...
    'MajorTicks',[],'MinorTicks',[], ...
    'ValueChangingFcn',@(~,ev)onMinChanging(ev.Value), ...
    'ValueChangedFcn', @(~,~)onWinRangeChanged());
sMax = uislider(tabUniformity,'Orientation','vertical', ...
    'Position',[SL_X_MAX SL_Y SL_W SL_H],'Limits',[0 1],'Value',1, ...
    'MajorTicks',[],'MinorTicks',[], ...
    'ValueChangingFcn',@(~,ev)onMaxChanging(ev.Value), ...
    'ValueChangedFcn', @(~,~)onWinRangeChanged());
uistack(sMin,'top'); uistack(sMax,'top');    

uilabel(tabUniformity,'Text','Window (min / max)','Position',[SL_X_MIN-10 SL_Y+SL_H+14 140 22]);
lbl0  = uilabel(tabUniformity,'Text','0','Position',[SL_X_MAX+22 SL_Y-2 120 18]);
lblMid= uilabel(tabUniformity,'Text','mid','Position',[SL_X_MAX+22 SL_Y+SL_H/2-9 120 18]);
lblMax= uilabel(tabUniformity,'Text','max','Position',[SL_X_MAX+22 SL_Y+SL_H-8 120 18]);
getWinValues = @() sort([sMin.Value sMax.Value]);

%% ------------------ Tab 5: Inspection (with smoothing controls) ------------------
leftW = 360; axW = 940;

uilabel(tabInspect,'Text','Inspect image:','Position',[24 750 110 24]);
ddInspectSrc = uidropdown(tabInspect,'Items',{'Uniformity','Attenuation','Grayscale Uniformity'}, ...
    'Position',[130 750 210 24],'Value','Uniformity');

btnShowImg = uibutton(tabInspect,'Text','Show Selected Image','Position',[24 715 316 28], ...
    'ButtonPushedFcn', @(~,~)showInspectImage());

btnArm     = uibutton(tabInspect,'Text','Start Picking (3 pts)','Position',[24 680 152 28], ...
    'ButtonPushedFcn', @(~,~)armInspectPicking());
uibutton(tabInspect,'Text','Clear points','Position',[188 680 152 28], ...
    'ButtonPushedFcn', @(~,~)clearInspectPoints());

uilabel(tabInspect,'Text','Range mode:','Position',[24 642 90 22]);
ddRangeMode = uidropdown(tabInspect,'Items',{'Manual','LEM (auto)','HEM (auto)'}, ...
    'Position',[115 642 225 24],'Value','LEM (auto)', ...
    'ValueChangedFcn', @(~,~)onRangeModeChanged());

uilabel(tabInspect,'Text','E(min) (keV):','Position',[24 610 90 22]);
eMinField = uieditfield(tabInspect,'numeric','Position',[115 610 90 24],'Value',70);
uilabel(tabInspect,'Text','E(max) (keV):','Position',[212 610 90 22]);
eMaxField = uieditfield(tabInspect,'numeric','Position',[300 610 40 24],'Value',120);

btnPlot3 = uibutton(tabInspect,'Text','Plot All Three Points','Position',[24 572 316 30], ...
    'ButtonPushedFcn', @(~,~)plotAllThreePoints());

lblPickStatus = uilabel(tabInspect,'Text','Points: 0/3','Position',[24 542 316 24]);


axInspectImg = uiaxes(tabInspect,'Position',[40 190 leftW-60 320]);

axInspectSpec1 = uiaxes(tabInspect,'Position',[leftW+40 540 axW 180]);
axInspectSpec2 = uiaxes(tabInspect,'Position',[leftW+40 335 axW 180]);
axInspectSpec3 = uiaxes(tabInspect,'Position',[leftW+40 130 axW 180]);
for a = [axInspectSpec1 axInspectSpec2 axInspectSpec3]
    xlabel(a,'Energy (keV)'); ylabel(a,'Counts'); grid(a,'on'); box(a,'on');
end

uilabel(tabInspect,'Text','Smoothing method:','Position',[leftW+40 90 140 22]);
ddSmooth = uidropdown(tabInspect,'Items',{'LOESS','moving','sgolay'}, ...
    'Position',[leftW+180 90 120 24],'Value','LOESS');
uilabel(tabInspect,'Text','Fractional span (0-1):','Position',[leftW+320 90 160 22]);
spanField = uieditfield(tabInspect,'numeric','Position',[leftW+480 90 70 24],'Value',0.2,'Limits',[0.01 0.9]);
uibutton(tabInspect,'Text','Apply smoothing to plotted spectra','Position',[leftW+570 88 260 28], ...
    'ButtonPushedFcn', @(~,~)applySmoothing());

data.inspect = struct('pts',zeros(0,3),'armed',false,'maxPts',3, ...
    'hIm',[],'markObjs',gobjects(0),'textObjs',gobjects(0), ...
    'lastPlots',{[]},'lastX',{[]},'lastRawY',{[]});

%% ------------------ Tab 6: Windowed Energy Images ------------------
uilabel(tabWindowedImages,'Position',[30 740 220 22],'Text','Energy Range 1 (keV):');
minE1 = uieditfield(tabWindowedImages,'numeric','Position',[30 715 90 30],'Value',70);
maxE1 = uieditfield(tabWindowedImages,'numeric','Position',[140 715 90 30],'Value',90);
uilabel(tabWindowedImages,'Position',[30 680 220 22],'Text','Energy Range 2 (keV):');
minE2 = uieditfield(tabWindowedImages,'numeric','Position',[30 655 90 30],'Value',91);
maxE2 = uieditfield(tabWindowedImages,'numeric','Position',[140 655 90 30],'Value',110);
uilabel(tabWindowedImages,'Position',[30 620 220 22],'Text','Energy Range 3 (keV):');
minE3 = uieditfield(tabWindowedImages,'numeric','Position',[30 595 90 30],'Value',111);
maxE3 = uieditfield(tabWindowedImages,'numeric','Position',[140 595 90 30],'Value',130);
uilabel(tabWindowedImages,'Position',[30 560 220 22],'Text','Energy Range 4 (keV):');
minE4 = uieditfield(tabWindowedImages,'numeric','Position',[30 535 90 30],'Value',131);
maxE4 = uieditfield(tabWindowedImages,'numeric','Position',[140 535 90 30],'Value',150);
uibutton(tabWindowedImages,'Text','Generate Windowed Images','Position',[30 490 200 30], ...
    'ButtonPushedFcn', @(~,~)generateWindowedImages());
uibutton(tabWindowedImages,'Text','Export Windowed PNGs','Position',[30 450 200 30], ...
    'ButtonPushedFcn', @(~,~)exportWindowedPNGs());
uilabel(tabWindowedImages,'Text','Minimum Resolvable Feature Size: 4.92 mm', ...
    'Position',[30 420 400 22],'FontWeight','bold');
axWin1 = uiaxes(tabWindowedImages,'Position',[300 520 440 300]);
axWin2 = uiaxes(tabWindowedImages,'Position',[780 520 440 300]);
axWin3 = uiaxes(tabWindowedImages,'Position',[300 120 440 300]);
axWin4 = uiaxes(tabWindowedImages,'Position',[780 120 440 300]);

%% ------------------ Tab 7: Windowed Count-Rate Plot ------------------
uilabel(tabWindowRate,'Text','Mode:', 'Position',[30 740 60 22]);
ddMode = uidropdown(tabWindowRate,'Items',{'Count rate','Transmission/A','Efficiency'}, ...
     'Position',[90 738 160 24],'Value','Count rate');
cbShowA = uicheckbox(tabWindowRate,'Text','Plot as A = -ln(T)', ...
    'Position',[270 738 160 24],'Value',false);

uilabel(tabWindowRate,'Text','Flux:', 'Position',[30 706 40 22]);
ddFluxKind = uidropdown(tabWindowRate,'Items',{'Fixed 1.2e4 n/cm^2/s','Custom'}, ...
    'Position',[75 704 190 24],'Value','Fixed 1.2e4 n/cm^2/s','ValueChangedFcn',@(~,~)toggleFluxField());
fluxCustomField = uieditfield(tabWindowRate,'numeric','Position',[275 704 120 24],'Value',1.2e4,'Enable','off');
cbEffThroughput = uicheckbox(tabWindowRate,'Text','Use object counts (throughput)', ...
    'Position',[410 704 220 24],'Value',false);

uibutton(tabWindowRate,'Text','Make Plot','Position',[30 670 200 30], ...
    'ButtonPushedFcn', @(~,~)plotWindowRates());
axEffPlot = uiaxes(tabWindowRate,'Position',[100 120 1160 520]);

%% ------------------ Tab 8: Transmission & Sigma (ROI) ------------------
uilabel(tabTrans,'Text','Thickness d (cm):','Position',[20 840 120 22]);
thickFieldCM = uieditfield(tabTrans,'numeric','Position',[140 840 100 22],'Value',data.thickness_cm,'Limits',[eps Inf], ...
    'ValueChangedFcn',@(e,~)setThickness(e.Value));
uilabel(tabTrans,'Text','Material:','Position',[260 840 80 22]);
ddTransMat = uidropdown(tabTrans,'Items', { ...
    'Lucite (PMMA)','Cadmium (metal)','Boric Acid (H3BO3)','Boric Acid Solution', ...
    'Gadovist Solution','Gadolinium (metal)','Custom', ...
    'Boron carbide (B4C)', ...
    'Lighter (butane fuel)', ...
    'Lighter (PMMA shell + butane fuel)', ...
    'Battery (Li anode only)', ...
    'Battery (entire)'}, ...
    'Position',[340 840 260 24], ...
    'Value','Lucite (PMMA)','ValueChangedFcn',@(~,~)togglePPM());
uilabel(tabTrans,'Text','PPM / Conc:','Position',[610 840 100 22]);
ppmField = uieditfield(tabTrans,'numeric','Position',[710 840 100 22],'Value',0,'Enable','off');
uibutton(tabTrans,'Text','Compute Transmission / Sigma','Position',[830 836 220 28], ...
    'ButtonPushedFcn',@(~,~)computeTransmissionTab());
uibutton(tabTrans,'Text','Export PNGs','Position',[1060 836 120 28], ...
    'ButtonPushedFcn',@(~,~)exportTransPNGs());
uilabel(tabTrans,'Text','Energy:', 'Position',[1190 840 60 22]);
ddTransE = uidropdown(tabTrans,'Items',{'Thermal 25.3 meV','Thermal 12.1 meV','Cold 5.1 meV','Custom (meV)'}, ...
    'Position',[1245 840 120 24],'Value','Thermal 25.3 meV','ValueChangedFcn',@(~,~)toggleTransE());
enTrans = uieditfield(tabTrans,'numeric','Position',[1245 812 120 22],'Value',25.3,'Enable','off');

panelComp = uipanel(tabTrans, ...
    'Title','Composite mix (fractions)', ...
    'Units','pixels', ...
    'Position',[0 0 480 54], ...  
    'Visible','off');

lblL_PMMA = uilabel(panelComp,'Text','PMMA:',     'Position',[12 16 48 22],'Visible','off');
fL_PMMA   = uieditfield(panelComp,'numeric',      'Position',[60 16 60 22],'Value',0.30,'Limits',[0 1],'Visible','off');
lblL_But  = uilabel(panelComp,'Text','Butane:',   'Position',[130 16 56 22],'Visible','off');
fL_But    = uieditfield(panelComp,'numeric',      'Position',[186 16 60 22],'Value',0.70,'Limits',[0 1],'Visible','off');

lblB_Shell = uilabel(panelComp,'Text','Fe shell:',      'Position',[256 16 64 22],'Visible','off');
fB_Shell   = uieditfield(panelComp,'numeric',           'Position',[320 16 54 22],'Value',0.10,'Limits',[0 1],'Visible','off');
lblB_Cath  = uilabel(panelComp,'Text','MnO2:',          'Position',[380 16 44 22],'Visible','off');
fB_Cath    = uieditfield(panelComp,'numeric',           'Position',[424 16 44 22],'Value',0.50,'Limits',[0 1],'Visible','off');

axTrans = uiaxes(tabTrans,'Position',[60 480 560 320]);
axSigm  = uiaxes(tabTrans,'Position',[760 480 560 320]);
axTHist = uiaxes(tabTrans,'Position',[60 140 560 280]);
axSHist = uiaxes(tabTrans,'Position',[760 140 560 280]);
lblTransInfo = uilabel(tabTrans,'Text','(Load open-beam & object, then Compute)','Position',[20 800 1320 22],'WordWrap','on');

transRoiPanel = uipanel(tabTrans,'Title','Sigma-ROI (draw on Sigma image)','Position',[20 20 1320 100]);
uibutton(transRoiPanel,'Text','Add ROI on Sigma Map','Position',[10 40 180 28], ...
    'ButtonPushedFcn',@(~,~)addSigmaRoi());
uibutton(transRoiPanel,'Text','Clear ROI','Position',[200 40 120 28], ...
    'ButtonPushedFcn',@(~,~)clearSigmaRoi());
txtExp = uitextarea(transRoiPanel,'Position',[340 10 960 70],'Editable','off','Value',{'(Press "Compute Transmission / Sigma", then draw an ROI on the Sigma image.)'});

data.trans = struct('roi',[],'roiMask',[],'hSigmaImg',[]);

%% ------------------ Tab 9 (now): Reference Model & Comparison ------------------
ddMaterial = uidropdown(tabMacroXS,'Items',{'Boron (nat)','Boron-10','Cadmium (nat)','Cadmium-113','Gadolinium (nat)','CdZnTe (CZT)','Custom'}, ...
    'Position',[160 830 220 24], 'Value','Cadmium-113', 'ValueChangedFcn',@(~,~)onMaterialChanged(), ...
    'Enable','off');  

ddPreset = uidropdown(tabMacroXS,'Items',{'Thermal 25 meV','Thermal 12 meV','Cold 8 meV','Custom (meV)'}, ...
    'Position',[160 800 160 24], 'Value','Thermal 12 meV', 'ValueChangedFcn',@(~,~)onPreset(), ...
    'Enable','off'); 

enField = uieditfield(tabMacroXS,'numeric','Position',[440 800 80 22], 'Value',16, ...
    'Enable','off'); 

thicknessField = uieditfield(tabMacroXS,'numeric','Position',[160 770 120 22], ...
    'Value',data.thickness_cm, 'Enable','off');  

uilabel(tabMacroXS,'Text','Material:','Position',[20 830 120 22]);
uilabel(tabMacroXS,'Text','Energy preset:','Position',[20 800 120 22]);
uilabel(tabMacroXS,'Text','Energy (meV):','Position',[360 800 90 22]);
uilabel(tabMacroXS,'Text','Thickness d (cm):','Position',[20 770 130 22]);

uilabel(tabMacroXS,'Text','Density rho (g/cm^3):','Position',[20 740 140 22]);
rhoField = uieditfield(tabMacroXS,'numeric','Position',[160 740 120 22],'Value',8.65);
uilabel(tabMacroXS,'Text','Atomic weight A (g/mol):','Position',[20 710 140 22]);
awField = uieditfield(tabMacroXS,'numeric','Position',[160 710 120 22],'Value',112.904);

uipanelCZT = uipanel(tabMacroXS,'Title','CZT parameters','Position',[20 620 520 80]);
uilabel(uipanelCZT,'Text','Zn mole fraction x in Cd_{1-x}Zn_xTe:','Position',[10 30 220 22]);
czt_x_Field = uieditfield(uipanelCZT,'numeric','Position',[240 30 60 22],'Value',0.10,'Limits',[0 1]);
uilabel(uipanelCZT,'Text','Density rho (g/cm^3):','Position',[320 30 120 22]);
czt_rho_Field = uieditfield(uipanelCZT,'numeric','Position',[430 30 60 22],'Value',5.78);

refPanel = uipanel(tabMacroXS,'Title','Reference Sigma options','Position',[20 470 520 140]);
useWestcott = uicheckbox(refPanel,'Text','Apply Westcott g-factor (on top of 1/sqrt(E))', ...
    'Position',[10 90 320 22],'Value',false);
uilabel(refPanel,'Text','g-factor:', 'Position',[340 90 60 22]);
gField = uieditfield(refPanel,'numeric','Position',[400 90 80 22],'Value',1.0);
useENDF = uicheckbox(refPanel,'Text','Use ENDF sigma(E) + Spectrum phi(E) (overrides presets)', ...
    'Position',[10 60 360 22],'Value',false);
uilabel(refPanel,'Text','CSV energy units:', 'Position',[10 30 120 22]);
ddUnits = uidropdown(refPanel,'Items',{'meV','eV'},'Position',[130 30 80 22],'Value','meV');
uilabel(refPanel,'Text','Emin (meV):','Position',[230 30 80 22]);
EminField = uieditfield(refPanel,'numeric','Position',[310 30 80 22],'Value',14);
uilabel(refPanel,'Text','Emax (meV):','Position',[400 30 80 22]);
EmaxField = uieditfield(refPanel,'numeric','Position',[480 30 30 22],'Value',18);

uibutton(tabMacroXS,'Text','Load sigma(E) CSV for selected material','Position',[20 430 260 28], ...
    'ButtonPushedFcn',@(~,~)loadENDFcsvForMaterial());
uibutton(tabMacroXS,'Text','Load Spectrum phi(E) CSV','Position',[290 430 200 28], ...
    'ButtonPushedFcn',@(~,~)loadSpectrumCsv());

refInfo = uitextarea(tabMacroXS,'Position',[20 390 520 34],'Editable','off', ...
    'Value',{'Mode: Preset sigma0 + 1/sqrt(E) (set g if needed). ENDF overrides if checked & files loaded.'});
uilabel(tabMacroXS,'Text','Known Sigma (cm^-1, optional manual):','Position',[20 360 220 22]);
knownXSField = uieditfield(tabMacroXS,'text','Position',[250 360 120 22],'Value','');

uibutton(tabMacroXS,'Text','Compute Reference & Compare','Position',[20 320 260 30], ...
    'ButtonPushedFcn',@(~,~)computeRefCompare());

xsResultLabel = uilabel(tabMacroXS,'Text','Result: Not calculated', ...
    'Position',[20 260 520 60],'WordWrap','on');

axXSPlot = uiaxes(tabMacroXS,'Position',[560 180 800 620]);

onMaterialChanged(); onPreset(); toggleFluxField(); toggleTransE();
thicknessField.Value = data.thickness_cm;
syncRefTabFromTrans();
%% ====================== FUNCTIONS ======================
    function setThickness(val)
        if ~(isfinite(val) && val>0), thickFieldCM.Value = data.thickness_cm; return; end
        data.thickness_cm = val;
        thicknessField.Value = val; % mirror to Tab 9
        syncRefTabFromTrans();  
        if isgraphics(axSigm) && isfield(data,'lastTrans') && ~isempty(data.lastTrans)
            computeTransmissionTab(); % recompute maps with new thickness
        end
    end

    function onPreset()
        val = ddPreset.Value;
        switch val
            case 'Thermal 25 meV'
                enField.Value = 25; enField.Enable = 'off';
            case 'Thermal 12 meV'
                enField.Value = 12; enField.Enable = 'off';
            case 'Cold 8 meV'
                enField.Value = 8;  enField.Enable = 'off';
            otherwise
                enField.Enable = 'on';
        end
    end

    function setPrefix(val)
        if ~(ischar(val) || isstring(val)), return; end
        val = char(string(val));
        if isempty(strtrim(val)), val = 'Object'; end
        data.prefix = val;
        if ~strcmp(prefixField1.Value,val), prefixField1.Value = val; end
        title(axObjFull,  sprintf('%s Spectra (Full energy range)', data.prefix));
        title(axObjTrunc, sprintf('%s Spectra (Selected energy range)', data.prefix));
        labelBigUniform.Text = data.prefix; 
    end

    function onMaterialChanged()
        mat = ddMaterial.Value;
        switch mat
            case 'Boron (nat)';    key = 'B_nat';
            case 'Boron-10';       key = 'B10';
            case 'Cadmium (nat)';  key = 'Cd_nat';
            case 'Cadmium-113';    key = 'Cd113';
            case 'Gadolinium (nat)'; key = 'Gd_nat';
            case 'CdZnTe (CZT)';   key = 'CZT';
            otherwise;             key = '';
        end
        isCZT = strcmp(key,'CZT');
        uipanelCZT.Enable = iff(isCZT,'on','off');
        rhoField.Enable   = iff(~isCZT,'on','off');
        awField.Enable    = iff(~isCZT,'on','off');
        if ~isCZT && ~isempty(key) && isfield(refProps,key)
            rhoField.Value = refProps.(key).rho;
            awField.Value  = refProps.(key).A;
        end
    end

   function togglePPM()
    m = ddTransMat.Value;

    % Old ppm enable for solutions:
    if any(strcmp(m,{'Boric Acid Solution','Gadovist Solution'}))
        ppmField.Enable = 'on';
    else
        ppmField.Enable = 'off';
    end

    % Show/hide composite panel & fields
    isLightComposite = strcmp(m,'Lighter (PMMA shell + butane fuel)');
    isBatteryEntire  = strcmp(m,'Battery (entire)');

    panelComp.Visible = iff(isLightComposite || isBatteryEntire,'on','off');

    % Lighter composite fields
    set([lblL_PMMA fL_PMMA lblL_But fL_But], 'Visible', iff(isLightComposite,'on','off'));

    % Battery composite fields
    set([lblB_Shell fB_Shell lblB_Cath fB_Cath], 'Visible', iff(isBatteryEntire,'on','off'));
    layoutCompositePanel();
    syncRefTabFromTrans();
end


    function syncRefTabFromTrans()
    if ~isvalid(tabMacroXS), return; end

    % 1) Thickness mirrors Tab 8
    thicknessField.Value = data.thickness_cm;

    % 2) Energy mirrors Tab 8 exactly; also pick a preset if close
    e = max(enTrans.Value, eps);
    enField.Value  = e;

    % try to show a meaningful preset name if close to common values
    if abs(e-25.3) < 0.6 || abs(e-25) < 0.6
        ddPreset.Value = 'Thermal 25 meV';
    elseif abs(e-12.1) < 0.6 || abs(e-12) < 0.6
        ddPreset.Value = 'Thermal 12 meV';
    elseif abs(e-8.0) < 1.2 || strcmp(ddTransE.Value,'Cold 5.1 meV')
        % cold values vary a lot across sources; show a cold preset if near
        ddPreset.Value = 'Cold 8 meV';
    else
        ddPreset.Value = 'Custom (meV)';
    end

    % Keep Tab 9 read-only by design
    ddPreset.Enable = 'off';
    enField.Enable  = 'off';

    % 3) Material mapping: Tab 8 -> Tab 9
    m8 = ddTransMat.Value;
    switch m8
        case 'Cadmium (metal)'
            ddMaterial.Value = 'Cadmium (nat)';
        case 'Gadolinium (metal)'
            ddMaterial.Value = 'Gadolinium (nat)';
        case {'Boric Acid (H3BO3)','Boric Acid Solution'}
            ddMaterial.Value = 'Boron (nat)';
        case 'Gadovist Solution'
            ddMaterial.Value = 'Gadolinium (nat)';
        case 'Lucite (PMMA)'
            ddMaterial.Value = 'Custom';   % no elemental match
        otherwise
            if any(strcmp(ddMaterial.Items, m8))
                ddMaterial.Value = m8;
            else
                ddMaterial.Value = 'Custom';
            end
    end
    ddMaterial.Enable = 'off';

    % Re-run material side-effects (updates rho/A if applicable)
    onMaterialChanged();
end


function toggleTransE()
        if strcmp(ddTransE.Value,'Custom (meV)')
            enTrans.Enable = 'on';
        else
            enTrans.Enable = 'off';
            switch ddTransE.Value
                case 'Thermal 25.3 meV', enTrans.Value = 25.3;
                case 'Thermal 12.1 meV', enTrans.Value = 12.1;
                case 'Cold 5.1 meV',     enTrans.Value = 5.1;
            end
        end
        syncRefTabFromTrans();
    end

    function [ok,msg,keV] = pick_keV(S)
        ok = false; msg=''; keV=[];
        cand = {'keV','Energy_keV','E_keV','energy_keV','Eaxis_keV','E'};
        for i=1:numel(cand)
            if isfield(S,cand{i})
                keV = S.(cand{i}); ok = true; return;
            end
        end
        msg = 'Could not find keV vector in file (looked for keV/Energy_keV/E_keV/etc.).';
    end

    function updateZoom()
        xmin = xminField.Value; xmax = xmaxField.Value;
        if xmin>=xmax
            uialert(fig,'X-min must be less than X-max.','Invalid Input'); return;
        end
        ax2.XLim = [xmin xmax];
        axObjTrunc.XLim = [xmin xmax];
    end

    function setBeamMode()
    data.beamMode = ddBeamMode.Value;
    switch data.beamMode
        case 'LEM'
            if xminField.Value>=xmaxField.Value || xminField.Value<60 || xmaxField.Value<200
                xminField.Value = 65; xmaxField.Value = 195; updateZoom();
            end
            % If user chose auto mode on Inspect tab, sync the fields
            if strcmp(ddRangeMode.Value,'LEM (auto)') || strcmp(ddRangeMode.Value,'HEM (auto)')
                ddRangeMode.Value = 'LEM (auto)'; onRangeModeChanged();
            end
        case 'HEM'
            if xminField.Value>=xmaxField.Value || xminField.Value<170 || xmaxField.Value<455
                xminField.Value = 180; xmaxField.Value = 450; updateZoom();
            end
            if strcmp(ddRangeMode.Value,'LEM (auto)') || strcmp(ddRangeMode.Value,'HEM (auto)')
                ddRangeMode.Value = 'HEM (auto)'; onRangeModeChanged();
            end
    end
    end

    function loadData()
        [file,path] = uigetfile('*.mat','Select Open-Beam .mat', data.lastOpenPath);
        if isequal(file,0), return; end
        S = load(fullfile(path,file));
        [ok,msg,keV] = pick_keV(S);
        if ~ok, uialert(fig,msg,'Error'); return; end
        if ~isfield(S,'LMmap_Ecorrect') || ~isfield(S,'LMmap_noEcorrect')
            uialert(fig,'Missing LMmap_Ecorrect / LMmap_noEcorrect in file.','Error'); return; end
        data.lastOpenPath = path; data.filepath = path;
        data.keV = keV(:).';
        data.LMmap_Ecorrect   = S.LMmap_Ecorrect;
        data.LMmap_noEcorrect = S.LMmap_noEcorrect;
        data.open_beam_image  = reshape(sum(data.LMmap_Ecorrect,2), NX, NY);

        cla(ax1); cla(ax2);
        sum_corr = sum(data.LMmap_Ecorrect,1);
        sum_unc  = sum(data.LMmap_noEcorrect,1);
        plot(ax1, data.keV, sum_unc, 'Color',COL_NOEC, 'DisplayName','No Energy Correction'); hold(ax1,'on');
        plot(ax1, data.keV, sum_corr,'Color',COL_EC,   'DisplayName','Energy Corrected'); hold(ax1,'off');
        legend(ax1,'Location','best'); grid(ax1,'on'); xlabel(ax1,'Energy (keV)'); ylabel(ax1,'Counts');
        title(ax1,'Full Beam Spectrum (full energy range)');

        plot(ax2, data.keV, sum_unc, 'Color',COL_NOEC, 'DisplayName','No Energy Correction'); hold(ax2,'on');
        plot(ax2, data.keV, sum_corr,'Color',COL_EC,   'DisplayName','Energy Corrected'); hold(ax2,'off');
        legend(ax2,'Location','best'); grid(ax2,'on'); xlabel(ax2,'Energy (keV)'); ylabel(ax2,'Counts');
        title(ax2,'Full Beam Spectrum (selected energy range)');
        ax2.XLim = [xminField.Value xmaxField.Value];

        imagesc(axOB, flipud(data.open_beam_image)); axis(axOB,'image'); colormap(axOB,'jet'); colorbar(axOB);
        title(axOB,'Original Image (Open-Beam, summed)');

        keVLabel.Text = sprintf('Found keV vector (%d channels)', numel(data.keV));
        loadStatus.Text = ['Loaded: ' file];
    end

    function loadObject()
        startDir = data.lastObjPath; 
        if isfield(data,'filepath'), startDir = data.filepath; end
        [file,path] = uigetfile('*.mat','Select the Object Data File', startDir);
        if isequal(file,0), return; end
        S = load(fullfile(path,file));
        [ok,msg,keV] = pick_keV(S);
        if ~ok, uialert(fig,msg,'Error'); return; end
        if ~isfield(S,'LMmap_Ecorrect') || ~isfield(S,'LMmap_noEcorrect')
            uialert(fig,'Missing LMmap_Ecorrect or LMmap_noEcorrect in object file.','Error'); return;
        end
        data.lastObjPath = path; data.objectPath = path;
        data.object_Ecorrect   = S.LMmap_Ecorrect;
        data.object_noEcorrect = S.LMmap_noEcorrect;
        data.objectImage       = reshape(sum(S.LMmap_Ecorrect,2), NX, NY);
        data.keV = keV(:).';
        objStatus.Text = ['Object Loaded: ' file];
        plotObjectSpectra();
    end

    function saveCorrectionMap()
        if ~isfield(data,'LMmap_Ecorrect')
            uialert(fig,'Load open-beam data first.','Error'); return;
        end
        I0 = data.open_beam_image;
        I0_roi = I0; I0_roi(~data.mask) = NaN;
        val = ddPxx2.Value;
        pct = str2double(regexprep(val,'P',''));
        if ~isfinite(pct), pct = 95; end
        norm_fac = robust_percentile(I0_roi(:), pct);
        if ~isfinite(norm_fac) || norm_fac<=0, norm_fac = 1; end
        correction_map = I0 ./ norm_fac;
        correction_map(~data.mask) = NaN;
        inside = data.mask & ~isnan(correction_map);
        correction_map(inside) = max(correction_map(inside), eps);
        data.correction_map = correction_map;

        save(fullfile(data.filepath,'Correction_Map.mat'),'correction_map');
        imagesc(ax3, flipud(correction_map)); axis(ax3,'image'); colormap(ax3,'jet'); colorbar(ax3);
        title(ax3, sprintf('Correction Map (Open-beam normalized to ROI %s)', val));
        saveStatus.Text = 'Correction Map saved.';
    end

    function exportCorrectionMapPNG()
        if ~isfield(data,'correction_map')
            uialert(fig,'No correction map available.','Error'); return;
        end
        exportgraphics(ax3, fullfile(data.filepath,'Correction_Map.png'));
        saveStatus.Text = 'Correction Map PNG exported.';
    end

    function plotObjectSpectra()
        if ~isfield(data,'object_Ecorrect') || ~isfield(data,'object_noEcorrect') || ~isfield(data,'keV')
            return;
        end
        keV = data.keV;
        sum_corr = sum(data.object_Ecorrect,1);
        sum_unc  = sum(data.object_noEcorrect,1);
        cla(axObjFull); cla(axObjTrunc);
        plot(axObjFull, keV, sum_unc, 'Color',COL_NOEC,'DisplayName','No Energy Correction'); hold(axObjFull,'on');
        plot(axObjFull, keV, sum_corr,'Color',COL_EC,  'DisplayName','Energy Corrected');    hold(axObjFull,'off');
        legend(axObjFull,'Location','best'); grid(axObjFull,'on');
        xlabel(axObjFull,'Energy (keV)'); ylabel(axObjFull,'Counts'); xlim(axObjFull,[keV(1) keV(end)]);
        title(axObjFull, sprintf('%s Spectra (Full energy range)', data.prefix));

        plot(axObjTrunc, keV, sum_unc, 'Color',COL_NOEC,'DisplayName','No Energy Correction'); hold(axObjTrunc,'on');
        plot(axObjTrunc, keV, sum_corr,'Color',COL_EC,  'DisplayName','Energy Corrected');    hold(axObjTrunc,'off');
        legend(axObjTrunc,'Location','best'); grid(axObjTrunc,'on');
        xlim(axObjTrunc,[xminField.Value xmaxField.Value]);
        xlabel(axObjTrunc,'Energy (keV)'); ylabel(axObjTrunc,'Counts');
        title(axObjTrunc, sprintf('%s Spectra (Selected energy range)', data.prefix));
    end

    function calcEfficiency()
        if ~isfield(data,'keV') || (~isfield(data,'LMmap_Ecorrect') && ~isfield(data,'open_beam_image'))
            uialert(fig,'Load open-beam data first.','Error'); return;
        end
        if ~isfield(data,'object_Ecorrect')
            uialert(fig,'(Optional) Load object if you want throughput too.','Info');
        end
        Phi = fluxField.Value; 
        Texp = timeField.Value;
        Abeam = data.beam_area_cm2;

        Nopen = sum( data.LMmap_Ecorrect(data.roi_idx, :), 'all' );
        Nobj  = NaN;
        if isfield(data,'object_Ecorrect')
            Nobj = sum( data.object_Ecorrect(data.roi_idx, :), 'all' );
        end

        eta_int = (Nopen/max(Texp,eps))/max(Phi*Abeam,eps);
        eta_sys = NaN;
        if ~isnan(Nobj)
            eta_sys = (Nobj/max(Texp,eps))/max(Phi*Abeam,eps);
        end

        if isnan(eta_sys)
            effLabel.Text = sprintf('Intrinsic efficiency (open-beam): %.4e   |  Counts_open: %d   |  Incident: %.3e', ...
                eta_int, round(Nopen), Phi*Abeam*Texp);
        else
            effLabel.Text = sprintf(['Intrinsic efficiency (open-beam): %.4e   |  Throughput (with object): %.4e   ' ...
                ' |  Counts_open: %d, Counts_obj: %d   |  Incident: %.3e'], ...
                eta_int, eta_sys, round(Nopen), round(Nobj), Phi*Abeam*Texp);
        end
    end

    %% ---------- UNIFORMITY TAB ----------
    function onPercentileChanged()
        val = ddPercentile.Value; % 'P100' | 'P95' | 'P85'
        p = str2double(regexprep(val,'P',''));
        if ~isfinite(p), p = 95; end
        data.lastUniformity.flat_percentile = p;
        if isfield(data,'objectImage') && isfield(data,'open_beam_image')
            generateUniformityImages();
        end
    end

    function generateUniformityImages()
        if ~isfield(data,'objectImage') || ~isfield(data,'open_beam_image')
            uialert(fig,'Load object image and open-beam first (Tab 1).','Error'); return;
        end
        I0_raw = data.open_beam_image;
        obj    = data.objectImage;

        if ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,'flat_percentile')
            data.lastUniformity.flat_percentile = 95;
        end
        pxx = data.lastUniformity.flat_percentile;
        I0_roi = I0_raw; I0_roi(~data.mask) = NaN;
        nf = robust_percentile(I0_roi(:), pxx);
        if ~isfinite(nf) || nf<=0, nf = 1; end

        cm = I0_raw ./ nf; cm(~data.mask) = NaN; cm(data.mask) = max(cm(data.mask),eps);
        uniform = obj ./ cm; uniform(~data.mask) = NaN;

        % --- Flat-field corrected transmission (recommended) ---
% cm was already computed above as: cm = I0_raw ./ nf;

I0corr = I0_raw ./ cm;     % open-beam corrected
Icorr  = obj    ./ cm;     % object corrected

Tmap = Icorr ./ max(I0corr, eps);
Tmap(~data.mask) = NaN;
Tmap(data.mask)  = min(max(Tmap(data.mask), eps), 1);

Amap = -log(Tmap);


        roi_vals = uniform(data.mask);
        mu = mean(roi_vals,'omitnan'); sd = std(roi_vals,'omitnan');
        p5 = percentile_no_tb(roi_vals,5);
        p50= percentile_no_tb(roi_vals,50);
        p95= percentile_no_tb(roi_vals,95);
        roiStatsLabel.Text = sprintf('ROI stats (Uniformity, P%d): mean=%.3g, SD=%.3g, P5=%.3g, P50=%.3g, P95=%.3g', pxx,mu,sd,p5,p50,p95);

        cla(axOrig); cla(axUniform); cla(axAtten); cla(axGray);
        data.hImg.orig    = imagesc(axOrig,   flipud(obj));     axis(axOrig,'image');   colormap(axOrig,'jet'); colorbar(axOrig);  title(axOrig,'Original Image');
        data.hImg.uniform = imagesc(axUniform,flipud(uniform)); axis(axUniform,'image');colormap(axUniform,'jet');colorbar(axUniform);title(axUniform,sprintf('Uniformity Corrected (P%d)',pxx));
        data.hImg.atten   = imagesc(axAtten,  flipud(Amap));    axis(axAtten,'image');  colormap(axAtten,'jet'); colorbar(axAtten); title(axAtten,'Attenuation Map (-ln(I/I_0))');

        data.lastUniformity.obj     = obj;
        data.lastUniformity.uniform = uniform;
        data.lastUniformity.atten   = Amap;
        data.lastUniformity.Tmap    = Tmap;
        data.lastUniformity.I0_raw = I0_raw;
        data.lastUniformity.gray    = [];
        data.lastUniformity.flat_n95= nf;
        data.lastUniformity.note    = sprintf('P%d flat-field; A=-ln(I/I0); arrays are NON-flipped',pxx);

        setClickable(data.hImg.orig,   'orig');
        setClickable(data.hImg.uniform,'uniform');
        setClickable(data.hImg.atten,  'atten');

        data.lastActive = 'uniform';
        activeImgLabel.Text = 'Active image: Uniformity';

        refreshGrayPreview();
        computeRoiCounts();
        initWinRange();
    end

    function refreshGrayPreview()
        if ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,'uniform') || isempty(data.lastUniformity)
            return;
        end
        mode = ddGrayMode.Value;
        cla(axGray);
        switch mode
            case 'Grayscale (Uniformity)'
                img = data.lastUniformity.uniform;
                data.hImg.gray = imagesc(axGray, flipud(img));
                axis(axGray,'image'); colormap(axGray,'gray'); colorbar(axGray);
                axGray.CLimMode = 'manual';               
                title(axGray,'Grayscale Uniformity');
            case 'Grayscale of Amap'
                img = data.lastUniformity.atten;
                data.hImg.gray = imagesc(axGray, flipud(img));
                axis(axGray,'image'); colormap(axGray,'gray'); colorbar(axGray);
                axGray.CLimMode = 'manual';                 
                title(axGray,'Grayscale Attenuation');
            otherwise % Transmission (T)
                img = data.lastUniformity.Tmap;
                data.hImg.gray = imagesc(axGray, flipud(img), [0 1]);
                axis(axGray,'image'); colormap(axGray,'gray'); colorbar(axGray);
                axGray.CLimMode = 'manual';                 
                title(axGray,'Grayscale Transmission (T)');
        end
        onWinRangeChanged();              
        drawnow limitrate                 
    end
    function exportUniformityPNGs()
    if ~isfield(data,'lastUniformity') || isempty(fieldnames(data.lastUniformity))
        uialert(fig,'Generate images first.','Info'); return;
    end
    if ~isfield(data,'objectPath') || ~isfolder(data.objectPath)
        uialert(fig,'Object path not set. Load object file first (Tab 1).','Error'); return;
    end

    modeNow = ddExportMode.Value;   
    prefix  = data.prefix;

    obj     = data.lastUniformity.obj;
    uniform = data.lastUniformity.uniform;
    Amap    = data.lastUniformity.atten;
    gray    = data.lastUniformity.gray; if isempty(gray), gray = uniform; end

    switch modeNow
        case 'Separate PNGs'
            save_png(flipud(obj),     fullfile(data.objectPath,[prefix '_original.png']),            [prefix ' - Original Image'],'jet',true);
            save_png(flipud(uniform), fullfile(data.objectPath,[prefix '_uniformity_corrected.png']),[prefix ' - Uniformity Corrected'],'jet',true);
            save_png(flipud(Amap),    fullfile(data.objectPath,[prefix '_attenuation_map.png']),     [prefix ' - Attenuation Map (-ln(I/I_0))'],'jet',true);
            save_png(flipud(gray),    fullfile(data.objectPath,[prefix '_grayscale.png']),           [prefix ' - Gray (' ddGrayMode.Value ')'],'gray',true);
            uialert(fig,sprintf('Exported 4 PNGs to:\n%s',data.objectPath),'Success');

        case 'Single panel (2x2)'
            f = figure('Visible','off','Color','w','Position',[100 100 1200 900]);
            tl = tiledlayout(f,2,2,'Padding','compact','TileSpacing','compact');

            nexttile(tl); imagesc(flipud(obj));     axis image off; colormap(gca,'jet');  colorbar; title([prefix ' - Original']);
            nexttile(tl); imagesc(flipud(uniform)); axis image off; colormap(gca,'jet');  colorbar; title([prefix ' - Uniformity']);
            nexttile(tl); imagesc(flipud(Amap));    axis image off; colormap(gca,'jet');  colorbar; title([prefix ' - A = -ln(I/I_0)']);
            nexttile(tl); imagesc(flipud(gray));    axis image off; colormap(gca,'gray'); colorbar; title([prefix ' - ' char(ddGrayMode.Value)]);

            [file,path] = uiputfile( ...
                {'*.png','PNG Image (*.png)';'*.pdf','PDF (*.pdf)';'*.tif','TIFF (*.tif)'}, ...
                'Save 2x2 Panel As', fullfile(data.objectPath,[prefix '_uniformity_panel.png']));
            if isequal(file,0), delete(f); return; end
            outPath = fullfile(path,file);
            try
                if endsWith(lower(outPath),'.pdf')
                    exportgraphics(tl, outPath, 'ContentType','vector');
                else
                    exportgraphics(tl, outPath, 'Resolution',300);
                end
                delete(f);
                try clipboard('copy', outPath); end 
                uialert(fig, sprintf('Saved panel to:\n%s\n(Path copied to clipboard)', outPath), 'Export complete');
            catch ME
                delete(f);
                uialert(fig, sprintf('Export failed:\n%s', ME.message), 'Export failed');
            end
    end
end


    function saveAmapMat()
        if ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,'atten')
            uialert(fig,'Generate images first.','Info'); return;
        end
        if ~isfield(data,'objectPath') || ~isfolder(data.objectPath)
            uialert(fig,'Object path not set. Load object file first (Tab 1).','Error'); return;
        end
        S = struct();
        S.Amap = data.lastUniformity.atten; 
        S.mask = data.mask;
        S.pixel_pitch_cm = data.PIX_CM;
        S.flat_norm_Pxx  = data.lastUniformity.flat_percentile;
        S.note = data.lastUniformity.note;
        S.NX = data.NX; S.NY = data.NY;
        outFile = fullfile(data.objectPath, [data.prefix '_Amap.mat']);
        save(outFile,'-struct','S');
        uialert(fig,['Saved: ' outFile],'Success');
    end

    function saveUniformMat()
        if ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,'uniform')
            uialert(fig,'Generate images first.','Info'); return;
        end
        if ~isfield(data,'objectPath') || ~isfolder(data.objectPath)
            uialert(fig,'Object path not set. Load object file first (Tab 1).','Error'); return;
        end
        S = struct();
        S.uniform = data.lastUniformity.uniform; % NON-flipped
        S.mask = data.mask;
        S.pixel_pitch_cm = data.PIX_CM;
        S.flat_norm_Pxx  = data.lastUniformity.flat_percentile;
        S.note = data.lastUniformity.note;
        S.NX = data.NX; S.NY = data.NY;
        outFile = fullfile(data.objectPath, [data.prefix '_uniform.mat']);
        save(outFile,'-struct','S');
        uialert(fig,['Saved: ' outFile],'Success');
    end

    function setClickable(hImg, tagname)
        hImg.HitTest = 'on'; 
        hImg.PickableParts = 'all';
        hImg.ButtonDownFcn = @(~,~)setActive(tagname);
    end

    function setActive(tagname)
        data.lastActive = tagname;
        activeImgLabel.Text = ['Active image: ' prettyTag(tagname)];
    end

    function addRoi(kind)
        if isempty(data.lastActive)
            uialert(fig,'Click on a preview image first to set it active.','Info'); return;
        end
        ax = axesFromTag(data.lastActive);
        if isempty(ax) || isempty(data.hImg.(data.lastActive))
            uialert(fig,'Generate the images first (Use "Generate Images").','Info'); return;
        end
        if strcmp(kind,'rect')
            obj = drawrectangle(ax,'StripeColor','y');
        else
            return;
        end
        n = numel(data.rois)+1;
        name = sprintf('ROI %d',n);
        entry = struct('name',name,'obj',obj,'type','rect','on',true,'mask',[],'axesTag',data.lastActive);
        data.rois(end+1) = entry; 
        roiList.Items = {data.rois.name};
        roiList.Value = name;
        addlistener(obj,'MovingROI', @(~,~) computeRoiCounts());
        addlistener(obj,'ROIMoved',  @(~,~) computeRoiCounts());
        computeRoiCounts();
    end

    function ax = axesFromTag(tagname)
        switch tagname
            case 'orig',    ax = axOrig;
            case 'uniform', ax = axUniform;
            case 'atten',   ax = axAtten;
            case 'gray',    ax = axGray;
            otherwise,      ax = [];
        end
    end

    function deleteSelectedRoi()
        if isempty(data.rois), return; end
        idx = builtin('find', strcmp({data.rois.name}, roiList.Value), 1);
        if isempty(idx), return; end
        try
            if isvalid(data.rois(idx).obj), delete(data.rois(idx).obj); end
        end
        data.rois(idx) = [];
        roiList.Items = {data.rois.name};
        if isempty(roiList.Items), roiList.Value = {}; else, roiList.Value = roiList.Items{1}; end
        computeRoiCounts();
    end

    function clearAllRois()
        for k=1:numel(data.rois)
            try
                if isvalid(data.rois(k).obj), delete(data.rois(k).obj); end
            end
        end
        data.rois = struct('name',{},'obj',{},'type',{},'on',{},'mask',{},'axesTag',{});
        roiList.Items = {};
        roiList.Value = {};
        roiOut.Value = {'ROI results cleared.'};
    end

    function computeRoiCounts()
        if ~isfield(data,'lastUniformity') || isempty(fieldnames(data.lastUniformity))
            roiOut.Value = {'(Generate images first)'}; 
            return;
        end
        if ~isfield(data,'rois') || isempty(data.rois)
            roiOut.Value = {'(No ROIs yet. Use Add ROI.)'}; 
            return;
        end
        txt = {};
        for k=1:numel(data.rois)
            entry = data.rois(k);
            if ~isvalid(entry.obj), continue; end
            imh = data.hImg.(entry.axesTag);
            maskDisplay = createMask(entry.obj, imh);
            maskData = flipud(maskDisplay);
            arr = data.lastUniformity.(entry.axesTag);
            if isempty(arr), arr = data.lastUniformity.uniform; end
            m = maskData & data.mask;
            vals = arr(m);
            s  = sum(vals,'omitnan');
            mu = mean(vals,'omitnan');
            sd = std(vals,'omitnan');
            good = vals(isfinite(vals));
            if isempty(good), mx = NaN; else, mx = max(good); end
            txt{end+1} = sprintf('%s on %s: sum=%.3g, mean=%.3g, SD=%.3g, max=%.3g, N=%d', ...
                 entry.name, prettyTag(entry.axesTag), s, mu, sd, mx, numel(vals)); %#ok<AGROW>
        end
        if isempty(txt), txt = {'(No ROIs yet. Use Add ROI.)'}; end
        roiOut.Value = txt;
    end

    function ignoreRectToZero()
        if isempty(data.lastActive)
            uialert(fig,'Click a preview image to make it active first.','Info'); return;
        end
        tag = data.lastActive;
        ax  = axesFromTag(tag);
        imh = data.hImg.(tag);
        if isempty(imh) || ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,tag)
            uialert(fig,'Generate images first.','Info'); return;
        end
        try
            r = drawrectangle(ax,'StripeColor','r','FaceAlpha',0.1); wait(r);
            mkDisp = createMask(r, imh); delete(r);
            mkData = flipud(mkDisp);
            arr = data.lastUniformity.(tag);
            arr(mkData & data.mask) = 0;
            data.lastUniformity.(tag) = arr;

            if strcmp(tag,'uniform') || strcmp(ddGrayMode.Value,'Grayscale (Uniformity)')
                if strcmp(ddGrayMode.Value,'Grayscale (Uniformity)'), data.lastUniformity.gray = arr; end
                data.hImg.gray.CData = flipud(data.lastUniformity.gray);
            end
            data.hImg.(tag).CData = flipud(arr);
            computeRoiCounts();
        catch ME
            uialert(fig, sprintf('Ignore failed: %s', ME.message), 'Error');
        end
    end

    function repopulateZerosFromPoint()
        if isempty(data.lastActive)
            uialert(fig,'Click a preview image to make it active first.','Info'); return;
        end
        tag = data.lastActive;
        ax  = axesFromTag(tag);
        imh = data.hImg.(tag);
        if isempty(imh) || ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,tag)
            uialert(fig,'Generate images first.','Info'); return;
        end
        try
            r = drawrectangle(ax,'StripeColor','g','FaceAlpha',0.1); wait(r);
            pos = r.Position; delete(r);
            if any(~isfinite(pos)) || pos(3)<=0 || pos(4)<=0
                uialert(fig,'Selection cancelled or invalid.','Info'); return;
            end
            cx = round(pos(1) + pos(3)/2);
            cyDisp = round(pos(2) + pos(4)/2);
            cx = max(1,min(NX,cx));
            cyDisp = max(1,min(NY,cyDisp));
            cyData = NY - cyDisp + 1;

            arr  = data.lastUniformity.(tag);
            seed = arr(cyData,cx);
            if ~isfinite(seed)
                uialert(fig,'Picked pixel is NaN/Inf. Choose another point.','Info'); return;
            end
            zmask = (arr==0);
            if ~any(zmask(:))
                uialert(fig,'No 0-valued pixels to repopulate in this image.','Info'); return;
            end
            arr(zmask) = seed;
            data.lastUniformity.(tag) = arr;

            if strcmp(tag,'uniform') || strcmp(ddGrayMode.Value,'Grayscale (Uniformity)')
                if strcmp(ddGrayMode.Value,'Grayscale (Uniformity)'), data.lastUniformity.gray = arr; end
                data.hImg.gray.CData = flipud(data.lastUniformity.gray);
            end
            data.hImg.(tag).CData = flipud(arr);
            computeRoiCounts();
        catch ME
            uialert(fig,sprintf('Repopulate failed: %s',ME.message),'Error');
        end
    end

    function t = prettyTag(tagname)
        switch tagname
            case 'orig',    t = 'Original';
            case 'uniform', t = 'Uniformity';
            case 'atten',   t = 'Attenuation';
            case 'gray',    t = 'Grayscale';
            otherwise,      t = tagname;
        end
    end

    %% ---------- Windowed & Plots ----------
    function generateWindowedImages()
    if ~isfield(data,'object_Ecorrect') || ~isfield(data,'LMmap_Ecorrect') || ~isfield(data,'keV')
        uialert(fig,'Required data (object, open-beam) not loaded.','Error'); return;
    end

    % --- NEW: auto ranges for HEM (4 equal, integer bins from 140 to 400) ---
    if isfield(data,'beamMode') && strcmp(data.beamMode,'HEM')
        edges = round(linspace(140,400,5));    % [140 205 270 335 400]
        % push to UI fields (so user sees them)
        minE1.Value = edges(1); maxE1.Value = edges(2);
        minE2.Value = edges(2); maxE2.Value = edges(3);
        minE3.Value = edges(3); maxE3.Value = edges(4);
        minE4.Value = edges(4); maxE4.Value = edges(5);
    end
    keV = data.keV;
    ranges = [minE1.Value maxE1.Value; minE2.Value maxE2.Value; ...
              minE3.Value maxE3.Value; minE4.Value maxE4.Value];
    axesList = {axWin1, axWin2, axWin3, axWin4};
   
        for i=1:4
            eMin = ranges(i,1); eMax = ranges(i,2);
            idx = keV>=eMin & keV<=eMax;
            if ~any(idx)
                uialert(fig, sprintf('No channels found for %g-%g keV.', eMin, eMax), 'Warning'); continue;
            end
            partial_obj = reshape(sum(data.object_Ecorrect(:,idx), 2), NX, NY);
            partial_ob  = reshape(sum(data.LMmap_Ecorrect(:,idx),  2), NX, NY);

            ob_roi = partial_ob; ob_roi(~data.mask) = NaN;
            nf = robust_percentile(ob_roi(:), data.lastUniformity.flat_percentile);
            if ~isfinite(nf) || nf<=0, nf = 1; end
            cm_win = partial_ob ./ nf; cm_win(~data.mask) = NaN; cm_win(data.mask) = max(cm_win(data.mask),eps);

            corrected = partial_obj ./ cm_win; corrected(~data.mask) = NaN;

            ax = axesList{i};
            imagesc(ax, flipud(corrected)); axis(ax,'image'); colormap(ax,'jet'); colorbar(ax);
            title(ax, sprintf('%g-%g keV (per-window flat-field, P%d)', eMin, eMax, data.lastUniformity.flat_percentile));
        end
    end

    function exportWindowedPNGs()
        if ~isfield(data,'object_Ecorrect')
            uialert(fig,'Load object data first.','Error'); return;
        end
        prefix = 'Windowed';
        exportgraphics(axWin1, fullfile(data.objectPath,[prefix '_Image1.png']), 'Resolution',300);
        exportgraphics(axWin2, fullfile(data.objectPath,[prefix '_Image2.png']), 'Resolution',300);
        exportgraphics(axWin3, fullfile(data.objectPath,[prefix '_Image3.png']), 'Resolution',300);
        exportgraphics(axWin4, fullfile(data.objectPath,[prefix '_Image4.png']), 'Resolution',300);
        uialert(fig,'Windowed images exported.','Success');
    end

    function toggleFluxField()
        if strcmp(ddFluxKind.Value,'Custom')
            fluxCustomField.Enable = 'on';
        else
            fluxCustomField.Enable = 'off';
        end
    end

    function plotWindowRates()
        if ~isfield(data,'object_Ecorrect') || ~isfield(data,'keV')
            uialert(fig,'Load object data first.','Error'); return;
        end
        if ~isfield(data,'LMmap_Ecorrect')
            uialert(fig,'Load open-beam data first.','Error'); return;
        end

        keV = data.keV;
        ranges = [minE1.Value maxE1.Value; minE2.Value maxE2.Value; minE3.Value maxE3.Value; minE4.Value maxE4.Value];
        labels = strings(1,4);
        for i=1:4, labels(i) = sprintf('%g-%g',ranges(i,1),ranges(i,2)); end
        Texp = max(timeField.Value, eps);
        modeSel = ddMode.Value;

        Nobj = zeros(1,4); Nopen = zeros(1,4);
        for i=1:4
            idx = keV>=ranges(i,1) & keV<=ranges(i,2);
            if any(idx)
                Nobj(i)  = sum( data.object_Ecorrect(data.roi_idx, idx), 'all' );
                Nopen(i) = sum( data.LMmap_Ecorrect(data.roi_idx,  idx), 'all' );
            else
                Nobj(i) = 0; Nopen(i) = 0;
            end
        end

        cla(axEffPlot); hold(axEffPlot,'on'); grid(axEffPlot,'on');

        switch modeSel
            case 'Count rate'
                rates = Nobj ./ Texp;
                errs  = sqrt(max(Nobj,0)) ./ Texp;
                bar(axEffPlot, 1:4, rates);
                errorbar(axEffPlot, 1:4, rates, errs, 'k.', 'LineWidth',1);
                ylabel(axEffPlot,'ROI Count-Rate (s^{-1})'); title(axEffPlot,'Bandwise Count Rate'); ylim(axEffPlot,'auto');
            case 'Transmission/A'
                T  = Nobj ./ max(Nopen, eps);
                relVar = (1./max(Nobj,eps)) + (1./max(Nopen,eps));
                dT = T .* sqrt(max(relVar,0));
                if cbShowA.Value
                    A  = -log(max(T, eps));
                    dA = dT ./ max(T,eps);
                    bar(axEffPlot, 1:4, A);
                    errorbar(axEffPlot, 1:4, A, dA, 'k.', 'LineWidth',1);
                    ylabel(axEffPlot,'Attenuation A = -ln(T)'); title(axEffPlot,'Bandwise Attenuation'); ylim(axEffPlot,'auto');
                else
                    bar(axEffPlot, 1:4, T);
                    errorbar(axEffPlot, 1:4, T, dT, 'k.', 'LineWidth',1);
                    ylabel(axEffPlot,'Transmission T = N_{obj}/N_{open}'); title(axEffPlot,'Bandwise Transmission'); ylim(axEffPlot,[0 1]);
                end
            otherwise % Efficiency
                if strcmp(ddFluxKind.Value,'Custom'), Phi = max(fluxCustomField.Value, eps); else, Phi = 1.2e4; end
                Abeam = data.beam_area_cm2;
                if cbEffThroughput.Value
                    rate_use = Nobj ./ Texp;       
                    labelY   = 'Throughput \eta_{sys}'; ttl = 'Bandwise Throughput (with object)';
                else
                    rate_use = Nopen ./ Texp;      
                    labelY   = 'Intrinsic efficiency \eta'; ttl = 'Bandwise Intrinsic Efficiency';
                end
                eta  = rate_use ./ max(Phi*Abeam, eps);
                derr = sqrt(max(rate_use*Texp,0))./Texp ./ max(Phi*Abeam, eps);
                bar(axEffPlot, 1:4, eta);
                errorbar(axEffPlot, 1:4, eta, derr, 'k.', 'LineWidth',1);
                ylabel(axEffPlot,labelY); title(axEffPlot,ttl); ylim(axEffPlot,'auto');
        end
        set(axEffPlot,'XTick',1:4,'XTickLabel',labels);
        hold(axEffPlot,'off');
    end

    %% ---------- Transmission (Tab 8) ----------
function computeTransmissionTab()
    % Requires: thickFieldCM, ddTransMat, ppmField, enTrans,
    %           axTrans, axSigm, axTHist, axSHist, txtExp,
    %           refProps, refSigmaThermal_b, NA, data.mask, NX, NY
    if ~isfield(data,'open_beam_image') || ~isfield(data,'object_Ecorrect')
        uialert(fig,'Load open-beam & object data first (Tab 1).','Error'); return;
    end

    % Sum object counts over channels into image
    I0 = data.open_beam_image;
    I  = reshape(sum(data.object_Ecorrect,2), NX, NY);
    I0(~data.mask) = NaN; I(~data.mask) = NaN;

   % --- Build the SAME correction map cm as Tab 4 ---
I0_roi = I0; I0_roi(~data.mask) = NaN;
nf = robust_percentile(I0_roi(:), data.lastUniformity.flat_percentile);
if ~isfinite(nf) || nf<=0, nf = 1; end

cm = I0 ./ nf;
cm(~data.mask) = NaN;
cm(data.mask)  = max(cm(data.mask), eps);

% --- Flat-field corrected transmission ---
I0corr = I0 ./ cm;
Icorr  = I  ./ cm;

Tmap = Icorr ./ max(I0corr, eps);
Tmap(~data.mask) = NaN;
Tmap(data.mask)  = min(max(Tmap(data.mask), eps), 1);

% --- Sigma map ---
t_cm = max(thickFieldCM.Value, eps);
Sigmap = -log(Tmap) ./ t_cm;
Sigmap(~data.mask) = NaN;


    % Draw maps + histograms (plain text, no LaTeX)
cla(axTrans); cla(axSigm); cla(axTHist); cla(axSHist);

% Transmission image
imagesc(axTrans, flipud(Tmap), [0 1]);
axis(axTrans,'image'); colormap(axTrans,'gray'); colorbar(axTrans);
title(axTrans,'Transmission T = I/I0','Interpreter','none');

% Sigma image
data.trans.hSigmaImg = imagesc(axSigm, flipud(Sigmap));
axis(axSigm,'image'); colormap(axSigm,'jet'); colorbar(axSigm);
title(axSigm,'Macroscopic Sigma (cm^-1)','Interpreter','none');

% --- Histograms ---
% Flatten, mask, and remove NaNs/Infs
Tvals = Tmap(data.mask);  Tvals = Tvals(isfinite(Tvals));
Svals = Sigmap(data.mask); Svals = Svals(isfinite(Svals));

% T histogram
histogram(axTHist, Tvals, 50);  % 50 bins; tweak if you like
grid(axTHist,'on'); box(axTHist,'on');
xlabel(axTHist,'T','Interpreter','none');
ylabel(axTHist,'Count','Interpreter','none');
title(axTHist,'T histogram','Interpreter','none');

% Sigma histogram
histogram(axSHist, Svals, 50);
grid(axSHist,'on'); box(axSHist,'on');
xlabel(axSHist,'Sigma (cm^-1)','Interpreter','none');
ylabel(axSHist,'Count','Interpreter','none');
title(axSHist,'Sigma histogram','Interpreter','none');

    % Stash for ROI workflow
    data.lastTrans = struct('Tmap',Tmap,'Sigmap',Sigmap);
    if ~isfield(data,'trans'), data.trans = struct(); end
    data.trans.roi = []; data.trans.roiMask = [];

    % Reference for current selection
    E_meV   = max(enTrans.Value, eps);
    matName = ddTransMat.Value;
    ppm     = ppmField.Value;
    [Sigma_ref, refNote] = getSigmaRef_E(matName, ppm, E_meV);
    data.results = struct();
    data.results.sigma_ref_cm1 = Sigma_ref;
    data.results.t_ref         = exp(-Sigma_ref * t_cm);
    data.results.note_ref      = refNote;
    data.results.has_roi       = false;

    txtExp.Value = { ...
    sprintf('Energy: %.3g meV    Thickness: %.3f cm', E_meV, t_cm), ...
    'Maps ready. Click "Add ROI on Sigma Map", then drag a rectangle to compute ROI statistics.', ...
    sprintf('Reference: Sigma_ref = %s cm^{-1}; T_ref = %s', ...
          numOrNaN(Sigma_ref), numOrNaN(data.results.t_ref)) ...
    };
end

function applySmoothing()
    % Need object_Ecorrect (per-pixel spectra), keV axis, and picked points
    if ~isfield(data,'object_Ecorrect') || isempty(data.object_Ecorrect)
        uialert(fig,'Load the OBJECT file first.','Error'); return;
    end
    if ~isfield(data,'keV') || isempty(data.keV)
        uialert(fig,'keV axis not found. Load data first.','Error'); return;
    end
    if ~isfield(data,'inspect') || isempty(data.inspect.pts)
        uialert(fig,'Pick up to 3 points, plot them, then click smoothing.','Info'); return;
    end

    keV = data.keV(:);
    [eMin, eMax, ~] = getInspectRange();
    inwin = keV>=eMin & keV<=eMax;
    if ~any(inwin)
        uialert(fig,'Selected energy window has no points.','Info'); return;
    end
    x = keV(inwin);

    method = string(ddSmooth.Value);     
    span   = spanField.Value;           
    if ~isfinite(span) || span<=0, span = 0.02; end

    mkWin = @(n) max(3, 2*floor(max(3, round(max(0.01,min(0.9,span))*n))/2)+1);

    axs = [axInspectSpec1, axInspectSpec2, axInspectSpec3];
    nPick = min(size(data.inspect.pts,1), 3);

    for k = 1:nPick
        a = axs(k);
        % Grab the pixel's raw spectrum
        pidx = data.inspect.pts(k,3);
        yfull = data.object_Ecorrect(pidx,:).';
        y = yfull(inwin);

        % Compute smoothed y2
        switch lower(method)
            case 'loess'
                frac = min(max(span,0.005),0.9);
                y2 = smooth(x, y, frac, 'loess');
            case 'moving'
                win = mkWin(numel(x));
                y2 = movmean(y, win, 'omitnan');
            otherwise  % 'sgolay'
                win = mkWin(numel(x));
                ord = min(3, win-2);
                y2  = sgolayfilt(y, ord, win);
        end

        cla(a);
        plot(a, x, y2, 'LineWidth',1.4, 'DisplayName','Smoothed');
        grid(a,'on'); xlabel(a,'Energy (keV)'); ylabel(a,'Counts');
        xlim(a,[eMin eMax]);
        title(a, sprintf('Pixel %d Spectra (smoothed)', k), 'Interpreter','none');
        legend(a,'Location','northeast');
    end
end

function addSigmaRoi()
    if ~isfield(data,'lastTrans') || isempty(data.lastTrans) || ~isfield(data.trans,'hSigmaImg') || isempty(data.trans.hSigmaImg)
        uialert(fig,'Compute Transmission/Sigma first.','Info'); return;
    end
    clearSigmaRoi();
    ax = axSigm;
    obj = drawrectangle(ax,'StripeColor','y');
    data.trans.roi = obj;
    addlistener(obj,'MovingROI', @(~,~) recomputeTransRoi());
    addlistener(obj,'ROIMoved',  @(~,~) recomputeTransRoi());
    recomputeTransRoi();
end

function clearSigmaRoi()
    if isfield(data,'trans') && isfield(data.trans,'roi') && ~isempty(data.trans.roi)
        try
            if isvalid(data.trans.roi), delete(data.trans.roi); end
        end
    end
    if ~isfield(data,'trans'), data.trans = struct(); end
    data.trans.roi = [];
    data.trans.roiMask = [];
    if ~isfield(data,'results'), data.results = struct(); end
    data.results.has_roi = false;
end

function recomputeTransRoi()
    if ~isfield(data.trans,'roi') || isempty(data.trans.roi) || ~isvalid(data.trans.roi), return; end
    if ~isfield(data,'lastTrans') || isempty(data.lastTrans), return; end

    % Display mask -> data mask
    mkDisp = createMask(data.trans.roi, data.trans.hSigmaImg);
    mkData = flipud(mkDisp) & data.mask;
    data.trans.roiMask = mkData;

    Tmap   = data.lastTrans.Tmap;
    Sigmap = data.lastTrans.Sigmap;

    valsS = Sigmap(mkData);
    valsT = Tmap(mkData);

    sigma_exp = mean(valsS, 'omitnan');
    t_mean    = mean(valsT, 'omitnan');

    t_cm = max(thickFieldCM.Value, eps);
    t_pred = exp(-sigma_exp * t_cm);

    data.results.has_roi        = true;
    data.results.sigma_exp_cm1  = sigma_exp;
    data.results.t_mean         = t_mean;
    data.results.t_pred         = t_pred;
   
    Sigma_ref = data.results.sigma_ref_cm1;
    T_ref     = data.results.t_ref;

    dSigmaPct = pctDiff(sigma_exp, Sigma_ref);
    dTpredPct = pctDiff(t_pred, T_ref);
    dTmeanPct = pctDiff(t_mean, T_ref);

    txtExp.Value = formatTransReport(enTrans.Value, t_cm, ...
        sigma_exp, t_mean, t_pred, Sigma_ref, T_ref, ...
        dSigmaPct, dTpredPct, dTmeanPct, data.results.note_ref);
end

function exportTransPNGs()
    if ~isfield(data,'open_beam_image') || ~isfield(data,'objectPath')
        uialert(fig,'Load data first.','Error'); return;
    end
    try
        exportgraphics(axTrans, fullfile(data.objectPath,'Transmission_Tmap.png'), 'Resolution',300);
        exportgraphics(axSigm,  fullfile(data.objectPath,'Transmission_SigmaMap.png'), 'Resolution',300);
        exportgraphics(axTHist, fullfile(data.objectPath,'Transmission_T_hist.png'), 'Resolution',300);
        exportgraphics(axSHist, fullfile(data.objectPath,'Transmission_Sigma_hist.png'), 'Resolution',300);
        uialert(fig,'Exported Transmission PNGs.','Success');
    catch ME
        uialert(fig, ['Export failed: ' ME.message], 'Error');
    end
end

function [Sigma_ref, note] = getSigmaRef_E(matName, ppm, E_meV)
    sE = sqrt(25.3 / max(E_meV, eps)); % 1/sqrt(E) scaling

    function S = Sigma_element_cm1(key, rho_override)
        if ~isfield(refProps,key) || ~isfield(refSigmaThermal_b,key), S = NaN; return; end
        rho = refProps.(key).rho; if nargin>=2 && ~isempty(rho_override), rho = rho_override; end
        A   = refProps.(key).A;
        N   = rho * NA / A;                     % atoms/cm^3
        sigma_b = refSigmaThermal_b.(key) * sE; % barns
        S = N * sigma_b * 1e-24;                % cm^-1
    end

    Sigma_ref = NaN; note = '';
    switch matName
        case 'Lucite (PMMA)'
            rho_PMMA = 1.18;
            M  = 5*refProps.C_nat.A + 2*refProps.O_nat.A + 8*refProps.H_nat.A;
            Nf = rho_PMMA * NA / M;
            sigmaMol_b = (5*refSigmaThermal_b.C_nat + 2*refSigmaThermal_b.O_nat + 8*refSigmaThermal_b.H_nat) * sE;
            Sigma_ref = Nf * sigmaMol_b * 1e-24;
            note = sprintf(' (PMMA, rho=%.2f g/cc, %.1f meV, 1/sqrt(E))', rho_PMMA, E_meV);

        case 'Cadmium (metal)'
            Sigma_ref = Sigma_element_cm1('Cd_nat');
            note = sprintf(' (Cd nat metal, %.1f meV, 1/sqrt(E))', E_meV);

        case 'Gadolinium (metal)'
            Sigma_ref = Sigma_element_cm1('Gd_nat');
            note = sprintf(' (Gd nat metal, %.1f meV, 1/sqrt(E))', E_meV);

        case 'Boric Acid (H3BO3)'
            rho = 1.435;
            M  = 3*refProps.H_nat.A + refProps.B_nat.A + 3*refProps.O_nat.A;
            Nf = rho * NA / M;
            sigmaMol_b = (3*refSigmaThermal_b.H_nat + refSigmaThermal_b.B_nat + 3*refSigmaThermal_b.O_nat) * sE;
            Sigma_ref = Nf * sigmaMol_b * 1e-24;
            note = sprintf(' (solid H3BO3, rho=%.3g g/cc, %.1f meV, 1/sqrt(E))', rho, E_meV);

        case 'Boric Acid Solution'
            f = max(min(ppm/1e6, 1), 0);        
            S_water = Sigma_water(sE);
            [S_h3bo3, ~] = getSigmaRef_E('Boric Acid (H3BO3)', 0, E_meV);
            if ~isfinite(S_h3bo3), S_h3bo3 = 0; end
            Sigma_ref = (1-f)*S_water + f*S_h3bo3;
            note = sprintf(' (solution approx, ppm=%.0f, %.1f meV, 1/sqrt(E))', ppm, E_meV);

        case 'Gadovist Solution'
            f = max(min(ppm/1e6, 1), 0);
            S_water = Sigma_water(sE);
            S_gd = Sigma_element_cm1('Gd_nat'); if ~isfinite(S_gd), S_gd = 0; end
            Sigma_ref = (1-f)*S_water + f*S_gd;
            note = sprintf(' (solution approx, ppm=%.0f as Gd, %.1f meV, 1/sqrt(E))', ppm, E_meV);

        case 'Boron carbide (B4C)'   
             rho = 2.52;
             M   = 4*refProps.B_nat.A + refProps.C_nat.A;
             Nf  = rho * NA / M;  % molecules/cm^3
             sigmaMol_b = (4*refSigmaThermal_b.B_nat + refSigmaThermal_b.C_nat) * sE;
             Sigma_ref = Nf * sigmaMol_b * 1e-24;
             note = sprintf(' (B4C, rho=%.2f g/cc, %.1f meV, 1/sqrt(E))', rho, E_meV);

        case 'Lighter (butane fuel)'
             rho = 0.58;
             M   = 4*refProps.C_nat.A + 10*refProps.H_nat.A;
             Nf  = rho * NA / M;
             sigmaMol_b = (4*refSigmaThermal_b.C_nat + 10*refSigmaThermal_b.H_nat) * sE;
             Sigma_ref = Nf * sigmaMol_b * 1e-24;
             note = sprintf(' (C4H10, rho=%.2f g/cc, %.1f meV, 1/sqrt(E))', rho, E_meV);

        case 'Lighter (PMMA shell + butane fuel)'
            f_PMMA = min(max(fL_PMMA.Value,0),1);
            f_But  = min(max(fL_But.Value,0),1);
            s = f_PMMA + f_But; if s<=0, f_PMMA=0.3; f_But=0.7; s=1; end
            f_PMMA = f_PMMA/s; f_But = f_But/s;

            rho_PMMA   = 1.18;
            M_PMMA     = 5*refProps.C_nat.A + 2*refProps.O_nat.A + 8*refProps.H_nat.A;
            Nf_PMMA    = rho_PMMA * NA / M_PMMA;
            sig_PMMA_b = (5*refSigmaThermal_b.C_nat + 2*refSigmaThermal_b.O_nat + 8*refSigmaThermal_b.H_nat) * sE;
            S_PMMA     = Nf_PMMA * sig_PMMA_b * 1e-24;

            rho_B      = 0.58;
            M_B        = 4*refProps.C_nat.A + 10*refProps.H_nat.A;
            Nf_B       = rho_B * NA / M_B;
            sig_B_b    = (4*refSigmaThermal_b.C_nat + 10*refSigmaThermal_b.H_nat) * sE;
            S_But      = Nf_B * sig_B_b * 1e-24;

            Sigma_ref  = f_PMMA * S_PMMA + f_But * S_But;
            note = sprintf(' (thickness-mix: PMMA %.0f%% / butane %.0f%%, 1/sqrt(E))', 100*f_PMMA, 100*f_But);

    case 'Battery (Li anode only)'
            Sigma_ref = Sigma_element_cm1('Li_nat');
            note = sprintf(' (Li metal nat, %.1f meV, 1/sqrt(E))', E_meV);

    case 'Battery (entire)'
            f_shell = min(max(fB_Shell.Value,0),1);
            f_cath  = min(max(fB_Cath.Value,0),1);
            f_elec  = 1 - f_shell - f_cath;
            if f_elec < 0, f_elec = 0; end
            s = f_shell + f_cath + f_elec; if s<=0, f_shell=0.1; f_cath=0.5; f_elec=0.4; s=1; end
            f_shell = f_shell/s; f_cath=f_cath/s; f_elec=f_elec/s;

    S_Fe = Sigma_element_cm1('Fe_nat');

    rho_MnO2   = 5.0;  % g/cc
    M_MnO2     = refProps.Mn_nat.A + 2*refProps.O_nat.A;
    Nf_MnO2    = rho_MnO2 * NA / M_MnO2;
    sig_MnO2_b = (refSigmaThermal_b.Mn_nat + 2*refSigmaThermal_b.O_nat) * sE;
    S_MnO2     = Nf_MnO2 * sig_MnO2_b * 1e-24;

    rho_PMMA   = 1.18;
    M_PMMA     = 5*refProps.C_nat.A + 2*refProps.O_nat.A + 8*refProps.H_nat.A;
    Nf_PMMA    = rho_PMMA * NA / M_PMMA;
    sig_PMMA_b = (5*refSigmaThermal_b.C_nat + 2*refSigmaThermal_b.O_nat + 8*refSigmaThermal_b.H_nat) * sE;
    S_PMMA     = Nf_PMMA * sig_PMMA_b * 1e-24;

    Sigma_ref = f_shell*S_Fe + f_cath*S_MnO2 + f_elec*S_PMMA;
    note = sprintf(' (thickness-mix: Fe %.0f%% / MnO2 %.0f%% / polymer %.0f%%, 1/sqrt(E))', ...
                   100*f_shell, 100*f_cath, 100*f_elec);

     
        otherwise
            Sigma_ref = NaN; note = '';
    end

    function S = Sigma_water(sEfac)
        M  = 2*refProps.H_nat.A + refProps.O_nat.A;
        Nf = 1.0 * NA / M;
        sigmaMol_b = (2*refSigmaThermal_b.H_nat + refSigmaThermal_b.O_nat) * sEfac;
        S = Nf * sigmaMol_b * 1e-24;
    end
end

%% ---------- Reference Model & Comparison (Tab 9) ----------
function computeRefCompare()
    % Compare experimental ROI Sigma from the Transmission tab with a reference model.
    if ~isfield(data,'results') || ~isfield(data.results,'has_roi') || ~data.results.has_roi
        xsResultLabel.Text = 'Result: No ROI Sigma(exp) yet. Go to "Transmission & Sigma" and draw an ROI on the Sigma map.';
        cla(axXSPlot); return;
    end
    sigma_exp = data.results.sigma_exp_cm1;

    [Sigma_known, infoStr] = computeReferenceSigma();
    cla(axXSPlot);
    if isfinite(Sigma_known)
        bar(axXSPlot, categorical({'Experimental','Reference'}), [sigma_exp, Sigma_known]);
        ttl = 'Sigma comparison';
        xsResultLabel.Text = sprintf('Experimental Sigma(exp, ROI) = %.4f cm^-1 | Reference Sigma(ref) = %.4f cm^-1%s', ...
                                     sigma_exp, Sigma_known, infoStr);
    else
        bar(axXSPlot, categorical({'Experimental'}), sigma_exp);
        ttl = 'Sigma (no reference)';
        xsResultLabel.Text = sprintf('Experimental Sigma(exp, ROI) = %.4f cm^-1%s', sigma_exp, infoStr);
    end
    ylabel(axXSPlot,'Macroscopic Cross Section (cm^-1)','Interpreter','none');
    title(axXSPlot, ttl,'Interpreter','none'); grid(axXSPlot,'on');
end

function [Sigma_known, infoStr] = computeReferenceSigma()
    % 1) Manual override
    tmp = str2double(strtrim(knownXSField.Value));
    if isfinite(tmp)
        Sigma_known = tmp; infoStr = ' | (Manual Sigma provided)'; return;
    end

    % 2) ENDF route
    if useENDF.Value
        [Sigma_known, infoStr] = computeFromENDF();
        if isfinite(Sigma_known), return; end
    end

    % 3) Preset (NIST 1/v) route
    [Sigma_known, infoStr] = computeFromNIST1v();
    if isfinite(Sigma_known), return; end

    % 4) FINAL FALLBACK: use Tab-8's reference (already computed there)
    if isfield(data,'results') && isfinite(data.results.sigma_ref_cm1)
        Sigma_known = data.results.sigma_ref_cm1;
        infoStr = ' | Using Transmission-tab reference';
        return;
    end

    % Nothing available
    Sigma_known = NaN;
    infoStr = ' | (No reference available for this selection)';
end


function [Sigma_known, infoStr] = computeFromNIST1v()
    mat     = ddMaterial.Value;
    E_meV   = max(enField.Value, eps);
    gfact   = max(gField.Value, 0);
    applyG  = useWestcott.Value;

    switch mat
        case 'Boron (nat)';    key='B_nat';
        case 'Boron-10';       key='B10';
        case 'Cadmium (nat)';  key='Cd_nat';
        case 'Cadmium-113';    key='Cd113';
        case 'Gadolinium (nat)'; key='Gd_nat';
        case 'CdZnTe (CZT)';   key='CZT';
        otherwise;             key='';
    end

    if strcmp(key,'CZT')
        x   = czt_x_Field.Value; rho = czt_rho_Field.Value;
        M   = (1-x)*refProps.Cd_nat.A + x*refProps.Zn_nat.A + refProps.Te_nat.A;
        Nform = rho * NA / M;
        s = sqrt(25.3 / E_meV); if applyG, s = s * gfact; end
        sigCd = refSigmaThermal_b.Cd_nat * s;
        sigZn = refSigmaThermal_b.Zn_nat * s;
        sigTe = refSigmaThermal_b.Te_nat * s;
        Sigma_known = ((1-x)*Nform*sigCd + x*Nform*sigZn + Nform*sigTe) * 1e-24;
        infoStr = sprintf(' | Preset sigma0 + 1/sqrt(E)%s (CZT, x=%.2f, rho=%.2f, E=%.1f meV)', ...
                          tern(applyG,' (g applied)',''), x, rho, E_meV);
        return;
    end

    if isempty(key) || ~isfield(refSigmaThermal_b,key)
        Sigma_known = NaN; infoStr = ' | (No sigma0 entry for this selection)'; return;
    end

    sigma0 = refSigmaThermal_b.(key);            % barns at 25.3 meV
    s  = sqrt(25.3 / E_meV); if applyG, s = s * gfact; end
    sigma_ref_b = sigma0 * s;

    rho = rhoField.Value; A = awField.Value;
    if ~(isfinite(rho) && rho>0 && isfinite(A) && A>0)
        Sigma_known = NaN;
        infoStr = sprintf(' | sigma(E)=%.3g b (no rho/A -> Sigma not computed)', sigma_ref_b);
    else
        N = rho * NA / A;
        Sigma_known = N * sigma_ref_b * 1e-24;
        infoStr = sprintf(' | Preset sigma0 + 1/sqrt(E)%s: sigma(E)=%.3g b; Sigma_ref=%.4f cm^-1 (E=%.1f meV)', ...
                          tern(applyG,' (g)',''), sigma_ref_b, Sigma_known, E_meV);
    end
end

function [Sigma_known, infoStr] = computeFromENDF()
    if isempty(data.spectrum) || ~isfield(data.spectrum,'E_meV')
        infoStr=' | ENDF mode: load a Spectrum phi(E) CSV first.'; Sigma_known=NaN; return;
    end
    Emin = EminField.Value; Emax = EmaxField.Value;
    if Emin >= Emax
        infoStr=' | ENDF mode: Emin<Emax required.'; Sigma_known=NaN; return;
    end
    specE = data.spectrum.E_meV(:); specF = data.spectrum.phi(:);
    inwin = specE>=Emin & specE<=Emax;
    if ~any(inwin), infoStr=' | ENDF mode: spectrum has no points in the window.'; Sigma_known=NaN; return; end
    specE = specE(inwin); specF = specF(inwin);

    mat = ddMaterial.Value;
    switch mat
        case 'CdZnTe (CZT)'
            need = {'Cd_nat','Zn_nat','Te_nat'};
            for k=1:numel(need)
                if isfield(data.endf,need{k})
                    infoStr = sprintf(' | ENDF (CZT): load sigma(E) for Cd, Zn, and Te. Missing: %s', need{k});
                    Sigma_known = NaN; return;
                end
            end
            sCd = interp1(data.endf.Cd_nat.E_meV, data.endf.Cd_nat.sigma_b, specE, 'linear','extrap');
            sZn = interp1(data.endf.Zn_nat.E_meV, data.endf.Zn_nat.sigma_b, specE, 'linear','extrap');
            sTe = interp1(data.endf.Te_nat.E_meV, data.endf.Te_nat.sigma_b, specE, 'linear','extrap');
            x=czt_x_Field.Value; rho=czt_rho_Field.Value;
            M=(1-x)*refProps.Cd_nat.A + x*refProps.Zn_nat.A + refProps.Te_nat.A;
            Nform=rho*NA / M;
            SigmaE=((1-x)*Nform.*sCd + x*Nform.*sZn + Nform.*sTe) * 1e-24;
            num=trapz(specE, SigmaE .* specF);
            den=trapz(specE, specF);
            Sigma_known = num / max(den, eps);
            infoStr = sprintf(' | ENDF mix (CZT): Sigma_ref=%.4f cm^-1 over [%.1f, %.1f] meV', Sigma_known, Emin, Emax);

        otherwise
            key = materialKeyForENDF(mat);
            if isempty(key) || ~isfield(data.endf,key)
                infoStr=' | ENDF mode: load sigma(E) CSV for the selected material.'; Sigma_known=NaN; return;
            end
            sE = interp1(data.endf.(key).E_meV, data.endf.(key).sigma_b, specE, 'linear','extrap');
            sigma_eff_b = trapz(specE, sE .* specF) / max(trapz(specE, specF), eps);
            rho=rhoField.Value; A=awField.Value;
            if ~(isfinite(rho)&&rho>0&&isfinite(A)&&A>0)
                Sigma_known=NaN;
                infoStr=sprintf(' | ENDF: sigma_eff=%.3g b (no rho/A -> Sigma not computed)', sigma_eff_b);
            else
                N=rho * NA / A;
                Sigma_known = N * sigma_eff_b * 1e-24;
                infoStr=sprintf(' | ENDF: sigma_eff=%.3g b; Sigma_ref=%.4f cm^-1 over [%.1f, %.1f] meV', ...
                                 sigma_eff_b, Sigma_known, Emin, Emax);
            end
    end
end

function key = materialKeyForENDF(mat)
    switch mat
        case 'Boron (nat)';     key='B_nat';
        case 'Boron-10';        key='B10';
        case 'Cadmium (nat)';   key='Cd_nat';
        case 'Cadmium-113';     key='Cd113';
        case 'Gadolinium (nat)';key='Gd_nat';
        case 'CdZnTe (CZT)';    key='CZT';
        otherwise;              key='';
    end
end

function loadENDFcsvForMaterial()
    mat=ddMaterial.Value;
    key=materialKeyForENDF(mat);
    if strcmp(key,'CZT')
        uialert(fig,'Select sigma(E) CSV for Cd (nat). Then Zn (nat). Then Te (nat).','CZT sigma(E)');
        loadOneENDF('Cd_nat'); loadOneENDF('Zn_nat'); loadOneENDF('Te_nat');
    else
        if isempty(key), uialert(fig,'Select a known material first.','Info'); return; end
        loadOneENDF(key);
    end
    refInfo.Value = {'ENDF sigma(E) loaded. If using ENDF mode, also load Spectrum phi(E) CSV.'};
end

function loadOneENDF(storeKey)
    [file,path]=uigetfile({'*.csv;*.txt','CSV/TXT files'}, 'Select sigma(E) CSV', data.lastCsvPath);
    if isequal(file,0), return; end
    data.lastCsvPath=path;
    T=readmatrix(fullfile(path,file));
    if isempty(T) || size(T,2)<2
        uialert(fig,'sigma(E) CSV must have at least 2 numeric columns.','Error'); return;
    end
    T=T(~any(isnan(T),2),:);
    E=T(:,1); sig=T(:,2);
    if strcmp(ddUnits.Value,'eV'), E=E*1000; end
    data.endf.(storeKey)=struct('E_meV',E,'sigma_b',sig);
end

function loadSpectrumCsv()
    [file,path]=uigetfile({'*.csv;*.txt','CSV/TXT files'}, 'Select Spectrum phi(E) CSV', data.lastCsvPath);
    if isequal(file,0), return; end
    data.lastCsvPath=path;
    T=readmatrix(fullfile(path,file));
    if isempty(T) || size(T,2)<2
        uialert(fig,'Spectrum CSV must have at least 2 numeric columns.','Error'); return;
    end
    T=T(~any(isnan(T),2),:);
    E=T(:,1); phi=T(:,2);
    if strcmp(ddUnits.Value,'eV'), E=E*1000; end
    data.spectrum=struct('E_meV',E,'phi',phi);
    refInfo.Value={sprintf('Loaded Spectrum with %d points (%s).', numel(E), ddUnits.Value)};
end

function layoutCompositePanel()
    try
        if ~exist('fig','var') || ~isvalid(fig), return; end
        if ~exist('panelComp','var') || ~isvalid(panelComp), return; end
        if ~exist('tabTrans','var') || ~isvalid(tabTrans), return; end

        % Figure size: [left bottom width height]
        fb = fig.Position;
        panW = 700; 
        panH = 60;

        % Center horizontally, keep a safe Y inside Tab 8
        panX = max(20, (fb(3) - panW)/2);
        panY = 760;   % your original Y; adjust if needed so it never overlaps

        panelComp.Position = [panX panY panW panH];
    catch ME
        warning('layoutCompositePanel failed: %s', ME.message);
    end
end



%% ------------------ Inspection: range helpers used elsewhere ------------------
function onRangeModeChanged()
    isManual = strcmp(ddRangeMode.Value,'Manual');
    eMinField.Enable = iff(isManual,'on','off');
    eMaxField.Enable = iff(isManual,'on','off');

    if ~isManual
        % Auto defaults
        if strcmp(ddRangeMode.Value,'LEM (auto)')
            eMinField.Value = 65;  eMaxField.Value = 195;
        else
            eMinField.Value = 180; eMaxField.Value = 450;  % HEM default updated
        end
    end
end

function [eMin, eMax, modeLabel] = getInspectRange()
    isManual = strcmp(ddRangeMode.Value,'Manual');
    if isManual
        eMin = eMinField.Value;  eMax = eMaxField.Value;  modeLabel = 'Manual';
    else
        if strcmp(ddRangeMode.Value,'LEM (auto)')
            eMin = 65;  eMax = 195; modeLabel = 'LEM(auto)';
        else
            eMin = 180; eMax = 450; modeLabel = 'HEM(auto)';  % HEM default updated
        end
        eMinField.Value = eMin; eMaxField.Value = eMax;
    end
    if ~(eMax > eMin), eMax = eMin + eps; end
end

function showInspectImage()
    if ~isfield(data,'lastUniformity') || isempty(fieldnames(data.lastUniformity))
        uialert(fig,'Generate images on Tab 4 first.','Info'); return;
    end
    switch ddInspectSrc.Value
        case 'Uniformity'
            img = data.lastUniformity.uniform; ttl = sprintf('%s - Uniformity',data.prefix); cm = 'jet';
        case 'Attenuation'
            img = data.lastUniformity.atten;   ttl = sprintf('%s - Attenuation',data.prefix); cm = 'jet';
        otherwise
            img = data.lastUniformity.uniform; ttl = sprintf('%s - Grayscale Uniformity',data.prefix); cm = 'gray';
    end
    cla(axInspectImg);
    data.inspect.hIm = imagesc(axInspectImg, flipud(img));
    axis(axInspectImg,'image'); colormap(axInspectImg,cm); colorbar(axInspectImg); title(axInspectImg,ttl);
    data.inspect.hIm.HitTest = 'on'; data.inspect.hIm.PickableParts = 'all';
    data.inspect.hIm.ButtonDownFcn = @onInspectImageDown;
    clearInspectPoints();
end

function armInspectPicking()
    data.inspect.armed = ~data.inspect.armed;
    if data.inspect.armed, btnArm.Text = 'Stop Picking';
    else,                  btnArm.Text = 'Start Picking (3 pts)'; end
end

function onInspectImageDown(~,~)
    if ~data.inspect.armed, return; end
    if size(data.inspect.pts,1) >= data.inspect.maxPts, return; end
    if isempty(data.inspect.hIm) || ~isgraphics(data.inspect.hIm), return; end
    if ~isfield(data,'object_Ecorrect')
        uialert(fig,'Load object file first.','Error'); return;
    end
    cp = axInspectImg.CurrentPoint;
    cx = round(cp(1,1)); cyDisp = round(cp(1,2));
    if ~(isfinite(cx) && isfinite(cyDisp)), return; end
    cx = max(1, min(NX, cx));
    cyDisp = max(1, min(NY, cyDisp));
    cyData = NY - cyDisp + 1;
    linIdx = cyData + (cx-1)*NX;

    data.inspect.pts = [data.inspect.pts; [cx cyData linIdx]]; 
    updatePickOverlay();

    if size(data.inspect.pts,1) >= data.inspect.maxPts
        data.inspect.armed = false;
        btnArm.Text = 'Start Picking (3 pts)';
    end
end

function updatePickOverlay()
    if ~isempty(data.inspect.markObjs)
        try, delete(data.inspect.markObjs(isvalid(data.inspect.markObjs))); end
    end
    if ~isempty(data.inspect.textObjs)
        try, delete(data.inspect.textObjs(isvalid(data.inspect.textObjs))); end
    end
    data.inspect.markObjs = gobjects(0);
    data.inspect.textObjs = gobjects(0);

    n = size(data.inspect.pts,1);
    for k=1:n
        cx = data.inspect.pts(k,1);
        cyData = data.inspect.pts(k,2);
        cyDisp = NY - cyData + 1;
        hold(axInspectImg,'on');
        mk = plot(axInspectImg, cx, cyDisp, 'wo','MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.5);
        tx = text(axInspectImg, cx+2, cyDisp, sprintf('P%d (%d,%d)',k,cyData,cx), ...
            'Color','w','FontWeight','bold','BackgroundColor',[0 0 0 0.25],'Margin',1);
        hold(axInspectImg,'off');
        data.inspect.markObjs(end+1) = mk;
        data.inspect.textObjs(end+1) = tx; 
    end
    lblPickStatus.Text = sprintf('Points: %d/3', n);
end

function clearInspectPoints()
    if isfield(data.inspect,'markObjs') && ~isempty(data.inspect.markObjs)
        try, delete(data.inspect.markObjs(isvalid(data.inspect.markObjs))); end
    end
    if isfield(data.inspect,'textObjs') && ~isempty(data.inspect.textObjs)
        try, delete(data.inspect.textObjs(isvalid(data.inspect.textObjs))); end
    end
    data.inspect.markObjs = gobjects(0);
    data.inspect.textObjs = gobjects(0);
    data.inspect.pts = zeros(0,3);
    localResetInspectAxes();
    try, lblPickStatus.Text = 'Points: 0/3'; catch, end
end

    function plotAllThreePoints()
    % Validate required data
    if ~isfield(data,'object_Ecorrect') || isempty(data.object_Ecorrect)
        uialert(fig,'Load the OBJECT file first (contains LMmap_Ecorrect).','Error'); return;
    end
    if ~isfield(data,'keV') || isempty(data.keV)
        uialert(fig,'keV axis not found. Load data first.','Error'); return;
    end
    if ~isfield(data,'inspect') || size(data.inspect.pts,1)==0
        uialert(fig,'Pick up to 3 points, then press "Plot All Three Points".','Info'); return;
    end

   
    keV_all = data.keV(:);                          
    [eMin, eMax, ~] = getInspectRange();
    inwin = (keV_all >= eMin) & (keV_all <= eMax);
    if ~any(inwin)
        uialert(fig,'Selected energy window has no channels.','Info'); return;
    end
    x = keV_all(inwin);                              

    
    axs = [axInspectSpec1, axInspectSpec2, axInspectSpec3];
    for a = axs
        cla(a); grid(a,'on'); xlabel(a,'Energy (keV)'); ylabel(a,'Counts');
        xlim(a,[eMin eMax]);
    end

    % Plot up to 3 picked pixels
    nPick = min(size(data.inspect.pts,1), 3);
    for k = 1:nPick
        pidx = data.inspect.pts(k,3);                
        y_raw = data.object_Ecorrect(pidx, inwin).'; 

        a = axs(k);
        y_use = y_raw;                              
        plot(a, x, y_use, 'LineWidth',1.2, 'DisplayName','Raw/Scaled'); hold(a,'on');

        % Shade the selected energy range
        yl = ylim(a);
        patch(a, [eMin eMax eMax eMin],[yl(1) yl(1) yl(2) yl(2)], ...
              [0.9 0.9 0.95], 'EdgeColor','none','FaceAlpha',0.25, ...
              'DisplayName','Range');

        hold(a,'off'); xlim(a,[eMin eMax]);
        title(a, sprintf('Pixel %d Spectra', k), 'Interpreter','none');
        legend(a,'Location','northeast');
    end

    % If fewer than 3 points, set titles/limits cleanly for the rest
    for k = (nPick+1):3
        a = axs(k);
        cla(a); grid(a,'on'); xlabel(a,'Energy (keV)'); ylabel(a,'Counts');
        xlim(a,[eMin eMax]);
        title(a, sprintf('Pixel %d Spectra', k), 'Interpreter','none');
    end
end


function localResetInspectAxes()
    try
        [eMin, eMax, ~] = getInspectRange();
        cla(axInspectSpec1); cla(axInspectSpec2); cla(axInspectSpec3);
        for a = [axInspectSpec1, axInspectSpec2, axInspectSpec3]
            grid(a,'on'); xlabel(a,'Energy (keV)'); ylabel(a,'Counts');
            xlim(a,[eMin eMax]);
        end
        title(axInspectSpec1, 'Pixel 1 Spectra', 'Interpreter','none');
        title(axInspectSpec2, 'Pixel 2 Spectra', 'Interpreter','none');
        title(axInspectSpec3, 'Pixel 3 Spectra', 'Interpreter','none');
    catch
    end
end

%% ------------------ Utility (percentiles etc.) ------------------
function p = robust_percentile(v, pct)
    v = v(isfinite(v));
    if isempty(v), p = NaN; return; end
    v = sort(v);
    if pct >= 100, p = v(end); return; end
    k = max(1, min(numel(v), round(pct/100 * numel(v))));
    p = v(k);
end

function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end

function out = iff(cond,a,b)
    if cond, out = a; else, out = b; end
end

function p = percentile_no_tb(x, pct)
    x = x(isfinite(x)); if isempty(x), p = NaN; return; end
    x = sort(x);
    r = pct/100*(numel(x)-1)+1;
    k = floor(r); d = r-k;
    if k >= numel(x), p = x(end); return; end
    p = x(k) + d*(x(k+1)-x(k));
end

function save_png(img, fname, ttl, cmap, doAxisImage)
    f = figure('Visible','off');
    imagesc(img); if doAxisImage, axis image off; else, axis off; end
    colormap(gca,cmap); colorbar; title(ttl,'Interpreter','none');
    exportgraphics(gca, fname, 'Resolution', 300);
    close(f);
end

%% ---- Windowing helpers ----
function initWinRange()
    if ~isfield(data,'lastUniformity') || ~isfield(data.lastUniformity,'uniform'), return; end
    U = data.lastUniformity.uniform;
    Um = U(data.mask); Um = Um(isfinite(Um));
    H = 1;
    if ~isempty(Um)
        H = max(Um); if ~isfinite(H) || H<=0, H = 1; end
    end
    setWinLimits(0, H);
    setWinValues([0 H]);
    onWinRangeChanged();
    try
        lbl0.Text   = sprintf('%.3g', 0);
        lblMid.Text = sprintf('%.3g', H/2);
        lblMax.Text = sprintf('%.3g', H);
    end
end

function setWinLimits(L, H)
    sMin.Limits = [L H];
    sMax.Limits = [L H];
    data.winLimits = [L H];
end

function setWinValues(v)
    v = sort(v(:).');
    v(1) = max(v(1), sMin.Limits(1));
    v(2) = min(v(2), sMax.Limits(2));
    if ~(v(2) > v(1)), v(2) = v(1) + eps; end
    sMin.Value = v(1);
    sMax.Value = v(2);
end

function onMinChanging(v)
    rng = diff(data.winLimits);
    db  = max(1e-6, 0.005*rng);
    if v >= sMax.Value - db, sMax.Value = min(v + db, sMax.Limits(2)); end
    sMin.Value = max(v, sMin.Limits(1));
    onWinRangeChanged();
end

function onMaxChanging(v)
    rng = diff(data.winLimits);
    db  = max(1e-6, 0.005*rng);
    if v <= sMin.Value + db, sMin.Value = max(v - db, sMin.Limits(1)); end
    sMax.Value = min(v, sMax.Limits(2));
    onWinRangeChanged();
end

function onWinRangeChanged()
    if ~isfield(data,'lastUniformity') || isempty(fieldnames(data.lastUniformity)), return; end
    v = sort([sMin.Value sMax.Value]); vmin=v(1); vmax=v(2);
    if ~(vmax > vmin), vmax = vmin + eps; setWinValues([vmin vmax]); end

    if isgraphics(axUniform) && ~isempty(data.hImg.uniform)
        caxis(axUniform,[vmin vmax]);
    end

    U  = data.lastUniformity.uniform;
    Um = U(data.mask); Um = Um(isfinite(Um)); if isempty(Um), Um = [0 1]; end
    uL = 0; uH = max(Um); if ~isfinite(uH) || uH<=0, uH = 1; end
    alpha = @(vv) (vv - uL) / max(uH - uL, eps);
    av = [alpha(vmin) alpha(vmax)]; av = min(max(av,0),1);

    if isfield(data.lastUniformity,'atten') && ~isempty(data.hImg.atten)
        A  = data.lastUniformity.atten;
        Am = A(data.mask); Am = Am(isfinite(Am)); if isempty(Am), Am = [0 1]; end
        Amin=min(Am); Amax=max(Am);
        aMinWin=Amin + av(1)*(Amax - Amin);
        aMaxWin=Amin + av(2)*(Amax - Amin);
        if ~(aMaxWin > aMinWin), aMaxWin=aMinWin + eps; end
        caxis(axAtten,[aMinWin aMaxWin]);
    end

    if ~isempty(data.hImg.gray)
        mode = ddGrayMode.Value;
        switch mode
            case 'Grayscale (Uniformity)'
                caxis(axGray,[vmin vmax]);
            case 'Grayscale of Amap'
                A  = data.lastUniformity.atten;
                Am = A(data.mask); Am = Am(isfinite(Am)); if isempty(Am), Am = [0 1]; end
                Amin=min(Am); Amax=max(Am);
                aMinWin=Amin + av(1)*(Amax - Amin);
                aMaxWin=Amin + av(2)*(Amax - Amin);
                if ~(aMaxWin > aMinWin), aMaxWin=aMinWin + eps; end
                caxis(axGray,[aMinWin aMaxWin]);
            otherwise % Transmission (T)
                T  = data.lastUniformity.Tmap;
                Tm = T(data.mask); Tm = Tm(isfinite(Tm)); if isempty(Tm), Tm = [0 1]; end
                Tmin=min(Tm); Tmax=max(Tm);
                tMinWin=Tmin + av(1)*(Tmax - Tmin);
                tMaxWin=Tmin + av(2)*(Tmax - Tmin);
                tMinWin=max(min(tMinWin,1),0);
                tMaxWin=max(min(tMaxWin,1),tMinWin+eps);
                caxis(axGray,[tMinWin tMaxWin]);
        end
    end
    drawnow limitrate
end

%% ---------- Tiny formatting helpers used above ----------
function s = numOrNaN(x, fmt)
    if nargin < 2, fmt = '%.4g'; end
    if ~isfinite(x), s = 'n/a'; else, s = sprintf(fmt, x); end
end

function p = pctDiff(meas, refv)
    if ~isfinite(meas) || ~isfinite(refv) || refv == 0
        p = NaN;
    else
        p = 100 * (meas - refv) / refv;
    end
end

function lines = formatTransReport(E_meV, t_cm, sigma_exp, T_mean, T_pred, Sigma_ref, T_ref, dSigmaPct, dTpredPct, dTmeanPct, refNote)
    if nargin < 11 || isempty(refNote), refNote = ''; end
    hasRef = isfinite(Sigma_ref);

    hdr  = sprintf('Energy: %.3g meV    Thickness: %.3f cm', max(E_meV,eps), max(t_cm,eps));
    exp1 = sprintf('ROI: Sigma_exp = %s cm^{-1}   |   T(mean) = %s   |   T(pred) = %s', ...
               numOrNaN(sigma_exp), numOrNaN(T_mean), numOrNaN(T_pred));

    if hasRef
       ref1   = sprintf('Reference: Sigma_ref = %s cm^{-1}; T_ref = %s', ...
                 numOrNaN(Sigma_ref), numOrNaN(T_ref));
        deltas = sprintf('Delta%% vs ref:  Sigma: %s%%   |   T_pred: %s%%   |   T_mean: %s%%', ...
                         numOrNaN(dSigmaPct), numOrNaN(dTpredPct), numOrNaN(dTmeanPct));
        lines = {hdr, exp1, ref1, deltas};
    else
        lines = {hdr, exp1, '(No reference model available for this selection.)'};
    end
end
end
