function [outV,gRate] = plotCrabtree(ecModel)

%% Define growth rate range
gRate_full = 0:0.01:1;
gRate = gRate_full;

%% Initialize output
outV = zeros(numel(ecModel.rxns), numel(gRate));
glcEx = getIndexes(ecModel,'EX_glc__D_e','rxns');       % Index of glucose exchange
ecModel = setParam(ecModel,'obj','EX_glc__D_e',1);       % Set glucose uptake as objective
totP  = -ecModel.lb(strcmp(ecModel.rxns,'prot_pool_exchange')); % Total protein pool

%% Loop over growth rates
for i = 1:numel(gRate)
    tmpModel = setParam(ecModel,'lb','BIOMASS_Ec_iML1515_core_75p37M',gRate(i));
    sol = solveLP(tmpModel);
    if ~isempty(sol.x)
        tmpModel = setParam(tmpModel,'lb','EX_glc__D_e',sol.x(glcEx)*1.01);
        tmpModel = setParam(tmpModel,'obj','prot_pool_exchange',1);
        sol = solveLP(tmpModel);
        outV(:,i) = sol.x;
    end
end

%% Truncate data to last valid solution
lastValidIndex = find(any(outV ~= 0, 1), 1, 'last');
gRate_plot = gRate_full(1:lastValidIndex);
outV_plot = outV(:, 1:lastValidIndex);

%% Load experimental data
fID = fopen(fullfile(findGECKOroot,'tutorials','Ecoli_ecModel','data','exp.tsv'),'r');
expData = textscan(fID,'%f %f %f %f','Delimiter',';','HeaderLines',2);
fclose(fID);
fluxToPlot = [expData{2} expData{3} expData{4}];

%% Reaction indices to plot
rxnsToPlot = getIndexes(ecModel,{'EX_o2_e','EX_glc__D_e','EX_ac_e','EX_pyr_e'},'rxns');

%% Colors and markers
fluxColors = [127,178,211; 141,211,201; 255,181,95; 252,127,113]/255; % Flux curves
poolColor = [110,160,180]/255;                                           % Protein pool
markerShapes = {'o','s','^','d'};

%% Plot
figure('Position',[100 100 900 400])
tiledlayout(1,2,'Padding','compact','TileSpacing','compact')

% Flux curves + experimental points
nexttile
hold on
hLines = gobjects(length(rxnsToPlot),1);
for r = 1:length(rxnsToPlot)
    hLines(r) = plot(gRate_plot, abs(outV_plot(rxnsToPlot(r),:)), ...
        'Color', fluxColors(r,:), 'LineWidth', 1.8);
    if r <= size(fluxToPlot,2)
        scatter(expData{1}, fluxToPlot(:,r), 60, fluxColors(r,:), markerShapes{r}, ...
            'filled', 'HandleVisibility','off')
    end
end
legend(hLines, ecModel.rxnNames(rxnsToPlot),'Location','northwest','FontSize',10)
xlabel('Growth rate (/hour)','FontSize',11)
ylabel('Absolute flux (mmol/gDCW·h)','FontSize',11)
xlim([0 1])
ylim([0 max(max(abs(outV_plot(rxnsToPlot,:))))*1.1])
ax = gca; ax.Box = 'off'; ax.TickDir='out'; ax.XColor=[0 0 0]; ax.YColor=[0 0 0];
ax.FontSize = 14;
hold off

% Protein pool usage
nexttile
poolRxn = getIndexes(ecModel,'prot_pool_exchange','rxns');
plot(gRate_plot, abs(outV_plot(poolRxn,:))/totP, 'LineWidth',1.8, 'Color', poolColor)
xlabel('Growth rate (/hour)','FontSize',11)
ylabel('Fraction of protein pool used','FontSize',11)
xlim([0 1])
ylim([0 1])
ax = gca; ax.Box='off'; ax.TickDir='out'; ax.XColor=[0 0 0]; ax.YColor=[0 0 0];
ax.FontSize = 14;

%% Output truncated data
outV = outV_plot;
gRate = gRate_plot;

end
