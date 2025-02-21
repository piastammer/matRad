%% Example: Photon Treatment Plan with Pencil Beam vs TOPAS MC
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
% (i) how to load patient data into matRad and setup a photon dose calculation 
% (ii) how to inversely optimize beamlet intensities
% (iii) how to write out files for TOPAS MC to recalculate the dose for the optimized plan 
% (iv) how to read the dose back in from the TOPAS result files
% (v) how to visually and quantitatively compare the results of TOPAS and matRad pencil beam

%% Patient Data Import
% Let's begin with a clear Matlab environment. Then, import the TG119
% phantom into your workspace. The phantom is comprised of a 'ct' and 'cst'
% structure defining the CT images and the structure set. Make sure the 
% matRad root directory with all its subdirectories is added to the Matlab 
% search path.

matRad_cfg = matRad_rc; %If this throws an error, run it from the parent directory first to set the paths
load('TG119.mat');

ixTarget = 3;

%% Treatment Plan
% The next step is to define your treatment plan labeled as 'pln'. This 
% matlab structure requires input from the treatment planner and defines 
% the most important cornerstones of your treatment plan.

pln.radiationMode   = 'photons';  
pln.machine         = 'Generic';
pln.numOfFractions  = 1; %one fraction for simplicity

%%
% Define the biological model used for modeling biological dose (esp. for
% particles).
% Possible biological models are:
% none:        use no specific biological model
% constRBE:    use a constant RBE
% MCN:         use the variable RBE McNamara model for protons
% WED:         use the variable RBE Wedenberg model for protons
% LEM:         use the biophysical variable RBE Local Effect model for carbons
% As we are  using photons, we simply set the parameter to 'none'
pln.bioModel = 'none';

%% 
% It is possible to request multiple error scenarios for robustness
% analysis and optimization. Here, we just use the "nominal scenario"
% (nomScen)
pln.multScen = 'nomScen';

% use a 50 degree geantry beam spacing
pln.propStf.gantryAngles = [0:50:359];
pln.propStf.couchAngles  = zeros(1,numel(pln.propStf.gantryAngles));
pln.propStf.numOfBeams   = numel(pln.propStf.gantryAngles);

stf                      = matRad_generateStf(ct,cst,pln);
pln.propStf.isoCenter    = vertcat(stf.isoCenter);
dij                      = matRad_calcDoseInfluence(ct,cst,stf,pln);
resultGUI         = matRad_fluenceOptimization(dij,cst,pln);

%% Dose Calculation
%resultGUI_MC = matRad_calcDoseInfluence(ct,cst,stf,pln);
pln.propDoseCalc.engine = 'TOPAS';
pln.propDoseCalc.beamProfile = 'uniform';
pln.propDoseCalc.numHistoriesDirect = 1e7; % probably dose has to be scaled somehow
pln.propDoseCalc.externalCalculation = 'write';
resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%% readout
%waitforbuttonpress; %We will wait since we do need to do the external calculation first
pln.propDoseCalc.externalCalculation = resultGUI_MC.meta.TOPASworkingDir;
resultGUI_MC = matRad_calcDoseForward(ct,cst,stf,pln,resultGUI.w);

%%  Visual Comparison of results
% Let's compare the new recalculation against the optimization result.
% Check if you have added all subdirectories to the Matlab search path,
% otherwise it will not find the plotting function
plane      = 3;
slice = matRad_world2cubeIndex(pln.propStf.isoCenter(1,:),ct);
slice = slice(3);
doseWindow = [0 max(resultGUI.physicalDose(:))];
doseWindow_MC = [0 max(resultGUI_MC.physicalDose(:))];

figure,title('original plan - matRad pencil beam')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow,[]);
figure,title('recalculated plan - TOPAS MC')
matRad_plotSliceWrapper(gca,ct,cst,1,resultGUI_MC.physicalDose,plane,slice,[],0.75,colorcube,[],doseWindow_MC,[]);

%% 
% At this point we would like to see the absolute difference of the first 
% optimization (pencil beam algorithm) and the recalculation (TOPAS MC)
absDiffCube = resultGUI.physicalDose-resultGUI_MC.physicalDose;
figure,title( 'difference pencil beam - MC plans')
matRad_plotSliceWrapper(gca,ct,cst,1,absDiffCube,plane,slice,[],[],colorcube);

%% Obtain dose statistics
% Two more columns will be added to the cst structure depicting the DVH and
% standard dose statistics such as D95,D98, mean dose, max dose etc.
resultGUI = matRad_planAnalysis(resultGUI,ct,cst,stf,pln);
resultGUI_MC = matRad_planAnalysis(resultGUI_MC,ct,cst,stf,pln);

%%
% The treatment plan using more beams should in principle result in a
% better OAR sparing. Therefore lets have a look at the D95 of the OAR of 
% both plans
ixOAR = 2;
disp(resultGUI.qi(ixOAR).D_95);
disp(resultGUI_MC.qi(ixOAR).D_95);
