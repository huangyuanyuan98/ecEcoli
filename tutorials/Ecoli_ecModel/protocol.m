% The function of this script is to demonstrate the reconstruction and
% analysis of a *full* ecModel. The tutorial uses the iML1515 model of
% E.coli as starting point. However, this script does not claim to 
% construct a "production-ready" ecEcoliGEM model:
% dependent on how you intend to use the ecModel it may require additional
% curation and evaluation.

%% Installation
% - Install RAVEN & GECKO
%   - Install RAVEN by following the installation instructions:
%     https://github.com/SysBioChalmers/RAVEN/wiki/Installation
%   - Simplest, RAVEN can be installed as MATLAB Add-On:
%     https://se.mathworks.com/help/matlab/matlab_env/get-add-ons.html
%   - The installation of Gurobi as LP solver is highly recommended
checkInstallation; % Confirm that RAVEN is functional, should be 2.9.2 or later.

%   - Install GECKO by following the installation instructions:
%     https://github.com/SysBioChalmers/GECKO/wiki/Installation-and-upgrade
%   - As RAVEN, GECKO can also be installed as MATLAB Add-On (link above)
%   - Add the appropriate GECKO (sub)folders to MATLAB path:
GECKOInstaller.install

%% STAGE 1: Preparation stage for ecModel reconstruction
% STEP 1: Preparation of project file structure and data
% - Initiate a basic structure of files and folders for your intended
%   project. This includes a copy of the template adapter.
%   The next line is commented out as the project structure is already
%   available in GECKO/tutorials/Ecoli_ecGEM.
% startGECKOproject()

% STEP 2: Store the starting GEM
% - Find a high-quality GEM of E. coli. ecEcoliGEM is
%   based on a high-quality E. coli GEM such as:
%   - iML1515 (http://bigg.ucsd.edu/static/models/iML1515.xml)
%   - Or other recent E. coli genome-scale models

% STEP 3: Modify the model adapter
% Model-specific parameters in the model adapter, which for this tutorial is
% located at tutorials/Ecoli_ecGEM/ecEcoliGEMadapter.m.


%% STAGE 2: Expansion from a starting metabolic model to an ecModel structure

% STEP 4: Set up model-specific adapter
% Specify the location of the EcoliGEMAdapter and set it as the default adapter.
adapterLocation = fullfile(findGECKOroot,'tutorials','Ecoli_ecModel','EcoliGEMAdapter.m');
ModelAdapterManager.setDefault(adapterLocation);

% Load the default ModelAdapter into the Workspace to access its parameters.
ModelAdapter = ModelAdapterManager.getDefault();
params = ModelAdapter.getParameters();
% Note: Any modifications to model parameters should be made in the
% adapter file itself (EcoliGEMAdapter.m), not directly in the Workspace.

% STEP 5: Load conventional E.coli-GEM
% Load the high-quality E. coli genome-scale metabolic model (iML1515)
% If the model path is already specified in the ModelAdapter (param.convGEM),
% simply use loadConventionalGEM().
model = loadConventionalGEM();
% Alternatively, the model can be imported directly from an SBML file:
% model = importModel(fullfile(findGECKOroot,'tutorials','Ecoli_ecGEM','models','iML1515.xml'));

% STEP 6: Prepare ecModel
% Convert the conventional GEM into an enzyme-constrained model (ecModel)
[ecModel, noUniprot] = makeEcModel(model,false);
% Verify UniProt annotation coverage
% If noUniprot is empty, all genes were successfully matched to UniProt dataset.

% STEP 7: Annotate ecModel with protein complex data
% Escherichia coli complexes are available in Complex Portal 
% and included in ecEcoliGEM.
[ecModel, foundComplex, proposedComplex] = applyComplexData(ecModel);

% STEP 8: Save the partially constructed ecModel in YAML format
% In this script there is code at the end of each stage to store a copy of
% the ecModel. 
saveEcModel(ecModel, 'ecEcoliModel_stage1.yml');

%% STAGE 3: integration of kcat into the ecModel structure
% ecModel=loadEcModel('ecEcoliModel_stage1.yml'); 

% STEP 9: Gather EC numbers
% First, extract existing EC numbers from the starting model. 
% Missing EC numbers are then retrieved from the UniProt and KEGG databases.
ecModel = getECfromGEM(ecModel);
noEC = cellfun(@isempty, ecModel.ec.eccodes);   % Identify reactions with missing EC numbers
ecModel = getECfromDatabase(ecModel, noEC);     % Update only the missing EC numbers

% STEP 10: Gather kcat values from BRENDA
% Retrieve kcat based on the EC numbers and substrate information of reactions.
% If no exact match is found, the function automatically relaxes 
% the matching criteria
kcatList_fuzzy = fuzzyKcatMatching(ecModel);
% Now we have a kcatList, which will be used to update ecModel in a later step.

% STEP 11-14: Gather kcat values from DLKcat
% STEP 11: Gather metabolite SMILES, required for DLKcat as kcat source
% Metabolite SMILES are gathered from PubChem.
[ecModel, noSmiles] = findMetSmiles(ecModel);

% STEP 12: Prepare DLKcat input file
% Generate the input file for DLKcat using enzyme sequences and metabolite SMILES.
% Note: Uncommenting the line below will overwrite existing DLKcat input files.
% writeDLKcatInput(ecModel, [], [], [], [], true);

% STEP 13: Run DLKcat
% Run the DLKcat prediction algorithm via Docker.
% Existing kcat values in DLKcat.tsv will be overwritten if the file exists.
% runDLKcat();

% STEP 14: Load DLKcat output
% Load the predicted kcat values and integrate them into the ecModel.
kcatList_DLKcat = readDLKcatOutput(ecModel);

% STEP 15-16: Combine kcat values from BRENDA and DLKcat
% Combine enzyme kinetic parameters from different sources based on priority rules.
% Current setting (1,1,0) retains only high-confidence BRENDA measurements 
% (origin level ≤1) and excludes EC wildcard matches. DLKcat predictions 
% supplement missing data to increase coverage.
% To increase the coverage of BRENDA data, the matching criteria can be relaxed.
kcatList_merged = mergeDLKcatAndFuzzyKcats(kcatList_DLKcat, kcatList_fuzzy, 1, 1, 0);
ecModel  = selectKcatValue(ecModel, kcatList_merged);

% STEP 17: Apply custom kcat values
% During the construction of the ecModel, high-confidence intracellular
% enzyme turnover numbers (kapp) were collected from literature to ensure
% that the model can accurately reproduce experimental phenotypes.
% These curated values are stored in the following files:
%   - customKcats.tsv: experimentally measured intracellular kapp from literature
%   - customKcats_ML.tsv: machine learning-predicted intracellular kapp based on literature data
% Apply kapp values calculated from literature
ecModel = applyCustomKcats(ecModel,fullfile(params.path,'data','customKcats.tsv'));

% Apply kapp values predicted by machine learning from literature
[ecModel, rxnUpdated, notMatch] = applyCustomKcats(ecModel,fullfile(params.path,'data','customKcats_ML.tsv'));

% STEP 18: Get kcat values across isozymes
% Some reactions are catalyzed by multiple isozymes, and certain isozymes may
% lack experimentally determined kcat values. This can lead to underestimated
% protein costs, causing these isozymes to be preferentially used in simulations.
% To avoid this bias, missing kcat values are filled with the mean kcat of
% the corresponding isozymes, providing a more reasonable estimate of the
% enzymatic parameters for all isozyme variants.
ecModel = getKcatAcrossIsozymes(ecModel);

% STEP 19: Apply kcat constraints from ecModel.ec.kcat to ecModel.S
% While STEP 16-18 add or modify values in ecModel.ec.kcat, these contraints
% are not yet applied to the S-matrix: the enzyme is not yet set as pseudo-
% substrate with the -MW/kcat stoichiometry. Now, applyKcatConstraints will
% take the values from ecModel.ec.kcat, considers any complex/subunit data
% that is tracked in ecModel.ec.rxnEnzMat, together with the MW in
% ecModel.ec.mw, and uses this to modify the enzyme usage coefficients
% directly in ecModel.S. Any time a change is made to the .kcat, .rxnEnzMat
% or .mw fields, the applyKcatConstraints function should be run again to
% reapply the new constraints onto the metabolic model.
ecModel = applyKcatConstraints(ecModel);

% STEP 20-21: Set upper bound of protein pool
% The protein pool usage in the model is constrained by three factors:
% 1. f-factor: ratio of enzymes/proteins
% 2. sigma-factor: how saturated enzymes are on average
% 3. Ptot: total protein content
% The theoretical maximum protein pool is: Protein pool = f × sigma × Ptot
% Here, sigma is set to 0.65 because most kcat values are derived from intracellular enzyme activity (kapp),
% which better reflects the actual enzyme saturation in vivo.
Ptot  = params.Ptot;
f     = params.f;
sigma = params.sigma; 

% Sets the limit of the total protein usage in the model
ecModel = setProtPoolSize(ecModel,Ptot,f,sigma);

% STEP 22: Save the STAGE-2 ecGEM to a YAML file
saveEcModel(ecModel, 'ecEcoliModel_stage2.yml');

%% STAGE 4: model tuning
% ecModel=loadEcModel('ecEcoliModel_stage2.yml');

% STEP 23: Test maximum growth rate
% Test whether the model is able to reach maximum growth if glucose uptake
% is unlimited. First set glucose uptake unconstraint.
ecModel = setParam(ecModel,'lb','EX_glc__D_e',-1000);
% And set growth maximization as the objective function.
ecModel = setParam(ecModel,'obj','BIOMASS_Ec_iML1515_core_75p37M',1);
% Run FBA.
sol = solveLP(ecModel,1);
bioRxnIdx = getIndexes(ecModel,params.bioRxn,'rxns');
fprintf('Growth rate: %f /hour.\n', sol.x(bioRxnIdx)) % Growth rate: 0.898699 /hour.

% STEP 24: Simulate Crabtree effect considering protein pool constraints
% We will below run a custom plotCrabtree function that is kept in the code
% subfolder. To run this function we will need to navigate into the folder
% where it is stored, but we will navigate back to the current folder
% afterwards.
currentFolder = pwd;
cd(fullfile(params.path,'code'))
% Plot Crabtree effect to analyze metabolic flux distribution
plotCrabtree(ecModel);

% STEP 25-26: Address pyruvate overflow
% Issue: The model exhibits pyruvate overflow
% Approach: Increase the kcat of reaction ACt2rpp_REV (acetate transport) by 100-fold
% Rationale: Enhances flux through the acetate transport reaction, which indirectly
% promotes pyruvate consumption and mitigates intracellular pyruvate accumulation
ecModel = setKcatForReactions(ecModel,'ACt2rpp_REV',551.74);
ecModel = applyKcatConstraints(ecModel);

% Plot Crabtree effect to analyze metabolic flux distribution
plotCrabtree(ecModel);

% STEP 27-28: Address insufficient acetate production
% Issue: Acetate production is lower than expected
% Solution: Reduce kcat values of respiration-related reactions 
% to redirect metabolic flux towards acetate formation 
% Modified reactions:
% ATPS4rpp: ATP synthase reactions, kcat reduced to limit respiratory ATP production
% ACONTa/ACONTb: Aconitase reactions, kcat reduced to limit TCA cycle flux
% These adjustments redirect metabolic flux towards acetate formation by constraining
% respiratory energy metabolism and TCA cycle activity
ecModel = setKcatForReactions(ecModel,'ATPS4rpp_EXP_1',82);
ecModel = setKcatForReactions(ecModel,'ATPS4rpp_EXP_2',82);
ecModel = setKcatForReactions(ecModel,'ACONTa_EXP_1',2.961);
ecModel = setKcatForReactions(ecModel,'ACONTa_EXP_2',3.093);
ecModel = setKcatForReactions(ecModel,'ACONTb_EXP_1',2.961);
ecModel = setKcatForReactions(ecModel,'ACONTb_EXP_2',3.093);
ecModel = applyKcatConstraints(ecModel);

% Plot Crabtree effect to analyze metabolic flux distribution
% Validation of adjustments: Acetate overflow can now be observed, 
% and the simulation results are close to experimental measurements
plotCrabtree(ecModel);

% STEP 29: Save the Crabtree effect plot
% The plot will also be saved in the output subfolder.
savePDF(gcf, fullfile(params.path,'output','crabtree.pdf'));

% STEP 30: Save the final ecModel
% Save the adjusted ecModel in YAML format for subsequent analysis and applications
saveEcModel(ecModel,'ecEcoliModel.yml');

fprintf('Final ecEcoliModel successfully constructed and saved.\n');
