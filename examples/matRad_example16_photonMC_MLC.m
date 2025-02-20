%% Example: Photon Treatment Plan using VMC++ dose calculation
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this example we will show 
% (i) how to load patient data into matRad
% (ii) how to setup a photon dose calculation based on the VMC++ Monte Carlo algorithm 
% (iii) how to inversely optimize the beamlet intensities directly from command window in MATLAB. 
% (iv) how to visualize the result

%% set matRad runtime configuration
matRad_rc %If this throws an error, run it from the parent directory first to set the paths

%% Patient Data Import
% Let's begin with a clear Matlab environment and import the boxphantom
% into your workspace. 
load('TG119.mat')
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.

pln.radiationMode           = 'photons';  
pln.machine                 = 'Generic';
pln.numOfFractions          = 1;
pln.propStf.gantryAngles    = [0:72:359];
pln.propStf.couchAngles     = [0 0 0 0 0];
%pln.propStf.gantryAngles    = [0];
%pln.propStf.couchAngles     = [0];
pln.propStf.bixelWidth      = 10;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
% Enable sequencing and direct aperture optimization (DAO).
pln.propOpt.runSequencing   = 1;
pln.propOpt.runDAO          = 1;

quantityOpt    = 'physicalDose';                                     
modelName      = 'none';  

% retrieve bio model parameters
pln.bioModel = matRad_bioModel(pln.radiationMode,modelName);

% retrieve scenarios for dose calculation and optimziation
pln.multScen = matRad_NominalScenario(ct);
% dose calculation settings
pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation

dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for IMRT
pln.propOpt.quantityOpt = quantityOpt;
resultGUI = matRad_fluenceOptimization(dij,cst,pln);
%% Sequencing
% This is a multileaf collimator leaf sequencing algorithm that is used in 
% order to modulate the intensity of the beams with multiple static 
% segments, so that translates each intensity map into a set of deliverable 
% aperture shapes.
resultGUI = matRad_siochiLeafSequencing(resultGUI,stf,dij,5,1);
[pln,stf] = matRad_aperture2collimation(pln,stf,resultGUI.sequencing,resultGUI.apertureInfo);
%% Aperture visualization
% Use a matrad function to visualize the resulting aperture shapes
matRad_visApertureInfo(resultGUI.apertureInfo)

%% Dose Calculation
%resultGUI_MC = matRad_calcDoseInfluence(ct,cst,stf,pln);
pln.propDoseCalc.engine = 'TOPAS';
pln.propDoseCalc.beamProfile = 'mlc';
pln.propDoseCalc.externalCalculation = 'write';
resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% readout
%waitforbuttonpress; %We will wait since we do need to do the external calculation first
pln.propDoseCalc.externalCalculation = resultGUI_MC.meta.TOPASworkingDir;
resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% Plot and compare the Resulting Dose Slices
plane      = 3;
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);
doseWindow = [0 max(resultGUI.physicalDose(:))];
doseWindow_MC = [0 max(resultGUI_MC.physicalDose(:))];

figure,title('original plan - matRad pencil beam')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('recalculated plan - TOPAS MC')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_MC.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow_MC,[]);