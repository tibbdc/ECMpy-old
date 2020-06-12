# This code is used to introduce enzyme concentration  constraint in GEMs
# by COBRApy and to calculate the parameters that need to be entered
# during the construction of the enzyme-constrained model.
from cobra.util.solver import set_objective
from cobra.core import Reaction
from cobra.io.dict import model_to_dict
import pandas as pd
from warnings import warn


def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    :param model: A Model object which will be modified in place.
    :return:
    """
    warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction._metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)


def get_genes_and_gpr(model, id_kcat_file):
    """Retrieving genes and gene_reaction_rule from GEM.

    :param model: A genome scale metabolic network model for
    constructing the enzyme-constrained model.
    :param id_kcat_file: A CSV file contains the kcat values for each
    reaction in the model.
    :return:All the genes in the model.
    """
    gpr = pd.DataFrame()
    model_news = model_to_dict(model, sort=False)
    df_id_kcat = pd.read_csv(id_kcat_file, index_col=0)
    genes = pd.DataFrame(model_news['genes'])
    genes.drop(['name', 'annotation'], axis=1, inplace=True)
    genes = genes[~genes['id'].isin(['s0001'])]
    genes.set_index('id', inplace=True)
    genes.to_csv("./genes_file.csv", sep=',')
    all_gpr = pd.DataFrame(model_news['reactions'])
    all_gpr.drop(['name', 'metabolites', 'lower_bound', 'upper_bound',
                  'annotation', 'objective_coefficient', 'notes'],
                 axis=1, inplace=True)
    all_gpr.to_csv("./all_GPR.csv", sep=',', index=False)
    df_all_gpr = pd.read_csv("./all_GPR.csv", index_col=0)
    for reaction_id in df_id_kcat.index:
        gpr.loc[reaction_id, 'GPR'] = \
            df_all_gpr.loc[reaction_id, 'gene_reaction_rule']
    gpr.to_csv("./ID_GPR_file.csv", sep=',')
    return genes


def calculate_mw(
        id_gpr_file, gene_subunit_file, subunit_molecular_weight_file):
    """Calculate the molecular weight of the enzyme that catalyzes each
    reaction in GEM based on the number of subunits and
    molecular weight of each gene.

    :param id_gpr_file: A CSV file contains the GPR relationship
     for each reaction in the GEM model.
    :param gene_subunit_file: A CSV file contains the number of subunit
     components of each gene expressed protein.
    :param subunit_molecular_weight_file: A CSV file contains the molecular
    weight of each gene expressed protein.
    :return: The molecular weight of the enzyme that catalyzes each reaction
     in the GEM model.
    """
    df_id_gpr = pd.read_csv(id_gpr_file, index_col=0)
    df_gene_subunit = pd.read_csv(gene_subunit_file, index_col=0)
    df_subunit_molecular_weight = \
        pd.read_csv(subunit_molecular_weight_file, index_col=0)
    id_mw = pd.DataFrame()
    for reaction_id in df_id_gpr.index:
        reaction_gpr = df_id_gpr.loc[reaction_id, 'GPR'].split('or')
        s = ''
        for i in range(0, len(reaction_gpr)):
            genes = reaction_gpr[i].replace('(', '').replace(")", '').\
                replace(" ", '').split('and')
            sum_mw = 0
            for m in range(0, len(genes)):
                # When the enzyme is composed of protein subunit expressed
                #  by only b3991, the number of subunit is 1.
                if genes[m] == 'b3991' and 'b3990' not in genes:
                    sum_mw += 1 * \
                              df_subunit_molecular_weight.loc[genes[m], 'mw']
                # When the enzyme is composed of protein subunits expressed by
                #  b3991 and b3990, the subunit composition ratio is 6:6.
                else:
                    sum_mw += df_gene_subunit.loc[genes[m], 'subunit'] * \
                              df_subunit_molecular_weight.loc[genes[m], 'mw']
            s = s + str(round(sum_mw, 4)) + ' or '
        s = s.rstrip(' or ')
        id_mw.loc[reaction_id, 'MW'] = s
        id_mw.to_csv("./ID_MW_file.csv", sep=',')
    return id_mw


def calculate_kcat_mw(id_kcat_file, id_mw):
    """Calculating kcat/MW

    When the reaction is catalyzed by several isozymes,
    the maximum was retained.

    :param id_kcat_file: A CSV file contains the kcat values for each
    reaction in the model.
    :param id_mw: The molecular weight of the enzyme that catalyzes
     each reaction in the GEM model.
    :return: The kcat/MW value of the enzyme catalyzing each reaction
     in the GEM model.
    """
    df_id_mw = id_mw
    df_id_kcat = pd.read_csv(id_kcat_file, index_col=0)
    id_kcat_mw = pd.DataFrame()
    for reaction_id in df_id_kcat.index:
        mw = df_id_mw.loc[reaction_id, 'MW'].split('or')
        min_mw = min(map(float, mw))
        kcat_mw = df_id_kcat.loc[reaction_id, 'kcat'] / min_mw
        id_kcat_mw.loc[reaction_id, 'kcat'] = \
            df_id_kcat.loc[reaction_id, 'kcat']
        id_kcat_mw.loc[reaction_id, 'MW'] = min_mw
        id_kcat_mw.loc[reaction_id, 'kcat_MW'] = kcat_mw
    id_kcat_mw.to_csv("./ID_kcat_MW_file.csv", sep=',')
    return id_kcat_mw


def calculate_f(genes, gene_abundance_file, subunit_molecular_weight_file):
    """Calculating f (the mass fraction of enzymes that are accounted
    in the model out of all proteins) based on the protein abundance
    which can be obtained from PAXdb database.

    :param genes: All the genes in the model.
    :param gene_abundance_file: The protein abundance of each gene
     in the E. coli genome.
    :param subunit_molecular_weight_file: The molecular weight of the
     protein subunit expressed by each gene.
    :return: The enzyme mass fraction f.
    """
    df_gene = genes
    df_gene_abundance = pd.read_csv(gene_abundance_file, index_col=0)
    df_subunit_molecular_weight = \
        pd.read_csv(subunit_molecular_weight_file, index_col=0)
    enzy_abundance = 0
    pro_abundance = 0
    for i in df_gene.index:
        enzy_abundance += df_gene_abundance.loc[i, 'abundance'] * \
                          df_subunit_molecular_weight.loc[i, 'mw']
    for i in df_gene_abundance.index:
        pro_abundance += df_gene_abundance.loc[i, 'abundance'] * \
                         df_subunit_molecular_weight.loc[i, 'mw']
    f = enzy_abundance/pro_abundance
    return f


def set_enzyme_constraint(model, id_kcat_mw, lowerbound, upperbound):
    """Introducing enzyme concentration constraint
    by COBRApy using the calculated parameters.

    :param model:
    :param id_kcat_mw: The kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model.
    :param lowerbound: The lower bound of enzyme concentration constraint in
     the enzyme-constrained model.
    :param upperbound: The upper bound of enzyme concentration constraint in
     the enzyme-constrained model.
    :return: Construct an enzyme-constrained model.
    """
    df_kcat_mw = id_kcat_mw
    coefficients = dict()
    for rxn in model.reactions:
        if rxn.id in df_kcat_mw.index:
            coefficients[rxn.forward_variable] = \
                1/float(df_kcat_mw.loc[rxn.id, 'kcat_MW'])
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model
