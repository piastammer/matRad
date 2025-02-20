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
load('BOXPHANTOM.mat');

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% structure requires input from the treatment planner and defines the most
% important cornerstones of your treatment plan.

pln.radiationMode           = 'photons';  
pln.machine                 = 'Generic';
pln.numOfFractions          = 30;
pln.bioModel = 'none'; 
pln.multScen = 'nomScen';

pln.propStf.gantryAngles    = [0];
pln.propStf.couchAngles     = [0];
pln.propStf.bixelWidth      = 10;
pln.propStf.numOfBeams      = numel(pln.propStf.gantryAngles);
pln.propStf.isoCenter       = ones(pln.propStf.numOfBeams,1) * matRad_getIsoCenter(cst,ct,0);
pln.propSeq.runSequencing   = 0;
pln.propOpt.runDAO          = 0;

% dose calculation settings
% We can choose a different dose calculation engine, here "ompMC", by
% setting the engine parameter in propDoseCalc. Further settings for the
% engine can be written into the propDoseCalc struct and will be passed on
% to the dose algorithm, if such a setting can be found
pln.propDoseCalc.engine = 'ompMC';

pln.propDoseCalc.doseGrid.resolution.x = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.y = 3; % [mm]
pln.propDoseCalc.doseGrid.resolution.z = 3; % [mm]

%% Generate Beam Geometry STF
stf = matRad_generateStf(ct,cst,pln);

%% Dose Calculation
% Calculate dose influence matrix for unit pencil beam intensities using 
% a Monte Carlo algorithm
dij = matRad_calcDoseInfluence(ct,cst,stf,pln);

%% Inverse Optimization for IMRT
resultGUI = matRad_fluenceOptimization(dij,cst,pln);

%% Plot the Resulting Dose Slice
% Just let's plot the transversal iso-center dose slice
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);
figure,
imagesc(resultGUI.physicalDose(:,:,slice)),colorbar, colormap(jet)

%%
% Exemplary, we show how to obtain the dose in the target and plot the histogram
ixTarget     = cst{2,4}{1};
doseInTarget = resultGUI.physicalDose(ixTarget);
figure
hist(doseInTarget);

 % use hist for compatibility with GNU Octave
title('dose in target'),xlabel('[Gy]'),ylabel('#');

%% compute integral energy
matRad_calcIntEnergy(resultGUI.physicalDose,ct,pln);
