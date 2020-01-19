#This code is used to introduce enzyme concentration  constraint in GEMs by COBRApy and to calculate the parameters that need to be entered during the construction of the enzyme-constrained model.

import cobra
from cobra.io.dict import model_to_dict, model_from_dict,metabolite_from_dict,gene_from_dict,reaction_from_dict
import pandas as pd
import csv
#Dividing the reversible reaction in GEM into two irreversible reactions.
def irreversible_model(model):
    cobra.manipulation.modify.convert_to_irreversible(model)
    return model

#Retrieving genes and gene_reaction_rule from GEM.
def get_genes_and_GPR(model,ID_kcat_file):
    GPR = pd.DataFrame()
    model_news = model_to_dict(model,sort=False)
    df_ID_kcat = pd.read_csv(ID_kcat_file,index_col=0)
    genes = pd.DataFrame(model_news['genes'])
    genes.drop(['name','annotation'],axis=1,inplace=True)
    genes = genes[~genes['id'].isin(['s0001'])]
    genes.set_index('id',inplace=True)
    genes.to_csv("./genes_file.csv", sep=',')
    all_GPR = pd.DataFrame(model_news['reactions'])
    all_GPR.drop(['name','metabolites','lower_bound','upper_bound','annotation','objective_coefficient','notes'],axis=1,inplace=True)
    all_GPR.to_csv("./all_GPR.csv", sep=',',index=False)
    df_all_GPR = pd.read_csv("./all_GPR.csv",index_col=0)
    for reaction_id in df_ID_kcat.index:
        GPR.loc[reaction_id,'GPR'] = df_all_GPR.loc[reaction_id,'gene_reaction_rule']
    GPR.to_csv("./ID_GPR_file.csv", sep=',')
    return(genes)

#Calculate the molecular weight of the enzyme that catalyzes each reaction in GEM based on the number of subunits and molecular weight of each gene.
def calculate_MW(ID_GPR_file,gene_subunit_file,subunit_molecular_weight_file):
    df_ID_GPR = pd.read_csv(ID_GPR_file,index_col=0)
    df_gene_subunit = pd.read_csv(gene_subunit_file,index_col=0)
    df_subunit_molecular_weight = pd.read_csv(subunit_molecular_weight_file,index_col=0)
    ID_MW = pd.DataFrame()
    for reaction_id in df_ID_GPR.index:
        reaction_GPR = df_ID_GPR.loc[reaction_id,'GPR'].split('or')
        s = ''
        for i in range(0, len(reaction_GPR)):
            genes = reaction_GPR[i].replace('(', '').replace(")", '').replace(" ", '').split('and')
            sum_mw = 0
            for m in range(0,len(genes)):
                if genes[m] == 'b3991' and 'b3990' not in genes:
                    sum_mw += 1 * df_subunit_molecular_weight.loc[genes[m],'mw']
                else:
                    sum_mw += df_gene_subunit.loc[genes[m],'subunit'] * df_subunit_molecular_weight.loc[genes[m],'mw']
            s = s + str(round(sum_mw,4)) + ' or '
        s = s.rstrip(' or ')
        ID_MW.loc[reaction_id, 'MW'] = s
        ID_MW.to_csv("./ID_MW_file.csv", sep=',')
    return(ID_MW)

#Calculating kcat/MW (When the reaction is catalyzed by several isozymes, the maximum was retained).
def calculate_kcat_MW(ID_kcat_file,ID_MW):
    df_ID_MW = ID_MW
    df_ID_kcat = pd.read_csv(ID_kcat_file,index_col=0)
    ID_kcat_MW = pd.DataFrame()
    for reaction_id in df_ID_kcat.index:
        MW = df_ID_MW.loc[reaction_id,'MW'].split('or')
        minmw = min(map(float,MW))
        kcat_mw = df_ID_kcat.loc[reaction_id,'kcat'] / minmw
        ID_kcat_MW.loc[reaction_id, 'kcat'] = df_ID_kcat.loc[reaction_id,'kcat']
        ID_kcat_MW.loc[reaction_id, 'MW'] = minmw
        ID_kcat_MW.loc[reaction_id, 'kcat_MW'] = kcat_mw
    ID_kcat_MW.to_csv("./ID_kcat_MW_file.csv", sep=',')
    return(ID_kcat_MW)

#Calculating f (the mass fraction of enzymes that are accounted in the model out of all proteins) based on the protein abundance which can be obtained from PAXdb database.
def calculate_f(genes,gene_abundance_file,subunit_molecular_weight_file):
    df_gene = genes
    df_gene_abundance = pd.read_csv(gene_abundance_file,index_col=0)
    df_subunit_molecular_weight = pd.read_csv(subunit_molecular_weight_file,index_col=0)
    enzy_abundance = 0
    pro_abundance = 0
    for i in df_gene.index:
        enzy_abundance += df_gene_abundance.loc[i,'abundance'] * df_subunit_molecular_weight.loc[i,'mw']
    for i in df_gene_abundance.index:
        pro_abundance += df_gene_abundance.loc[i,'abundance'] * df_subunit_molecular_weight.loc[i,'mw']
    f = enzy_abundance/pro_abundance
    return f

#Introducing enzyme concentration constraint by COBRApy using the calculated parameters.
def set_enzyme_constraint(model,ID_kcat_MW,lowerbound,upperbound):
    df_kcat_MW = ID_kcat_MW
    coefficients = dict()
    for rxn in model.reactions:
        if rxn.id in df_kcat_MW.index:
            coefficients[rxn.forward_variable] = 1/float(df_kcat_MW.loc[rxn.id,'kcat_MW'])
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model