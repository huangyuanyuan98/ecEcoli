classdef EcoliGEMAdapter < ModelAdapter 
	methods
		function obj = EcoliGEMAdapter()
			obj.params.path = fullfile(findGECKOroot,'tutorials','Ecoli_ecModel');

			obj.params.convGEM = fullfile(obj.params.path,'models','iML1515.xml');

			obj.params.sigma = 0.65; % kapp

			obj.params.Ptot = 0.56; % literature

			obj.params.f = 0.5524;
			
			obj.params.gR_exp = 0.69; % literature

			obj.params.org_name = 'escherichia coli';
			
			obj.params.complex.taxonomicID = 83333;

			obj.params.kegg.ID = 'eco';

			obj.params.kegg.geneID = 'kegg';

			obj.params.uniprot.type = 'proteome';

			obj.params.uniprot.ID = 'UP000000625';

			obj.params.uniprot.geneIDfield = 'gene_oln';

			obj.params.uniprot.reviewed = true;

			obj.params.c_source = 'EX_glc__D_e'; 

			obj.params.bioRxn = 'BIOMASS_Ec_iML1515_core_75p37M';

			obj.params.enzyme_comp = 'cytosol';			
		end
	
		function [spont,spontRxnNames] = getSpontaneousReactions(obj,model)
			spont = contains(model.rxnNames,'spontaneous');
			spontRxnNames = model.rxnNames(spont);
		end
		
	end
end
