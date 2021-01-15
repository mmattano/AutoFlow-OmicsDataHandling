import numpy as np
import pandas as pd
from pyopenms import FeatureMap, FeatureXMLFile
from cobra.core.formula import Formula

__version__ = "0.0.1"


INVALID_FORMULA_STR = ['(', 'Generic', 'R', 'X']


def extractNamesAndIntensities(feature_dir, sample_names, database):
    metabolites_unique = database[0].unique()
    extracted_data_dict = {}
    cnt = 0
    for name in sample_names:
        features = FeatureMap()
        FeatureXMLFile().load(
            feature_dir + "/" + name + ".featureXML", features
        )
        for f in features:
            try:
                peptideRef = f.getMetaValue("PeptideRef").decode("utf-8")
            except AttributeError:
                peptideRef = f.getMetaValue("PeptideRef")
            if (
                peptideRef
                in metabolites_unique
            ):
                formula = database[
                    database[0] == peptideRef
                ][1]
                extracted_data_dict[cnt] = {
                    "sample_group_name": name,
                    "Metabolite": peptideRef,
                    "Formula": list(formula)[0],
                    "Intensity": f.getMetaValue("peak_apex_int"),
                }
                cnt = cnt + 1
    return pd.DataFrame.from_dict(extracted_data_dict, "index")


def calculateMeanVarRSD(
    extracted_data_all, sample_name_2_replicate_groups, min_reps=3
):
    stats_all_dict = {}
    cnt = 0
    for replicate_group in sample_name_2_replicate_groups[
        "replicate_group_name"
    ].unique():
        list_of_compound_intensities = pd.DataFrame(
            columns=["Metabolite", "Formula", "Intensities"]
        )
        for sample_name_index, sample_name in sample_name_2_replicate_groups[
            sample_name_2_replicate_groups["replicate_group_name"]
            == replicate_group
        ].iterrows():
            for met_index, met in extracted_data_all[
                extracted_data_all["sample_group_name"]
                == sample_name["sample_group_name"]
            ].iterrows():
                list_of_compound_intensities = (
                    list_of_compound_intensities.append(
                        pd.DataFrame(
                            [
                                {
                                    "Metabolite": met["Metabolite"],
                                    "Formula": met["Formula"],
                                    "Intensities": met["Intensity"],
                                }
                            ]
                        ),
                        ignore_index=True,
                    )
                )
        for met_name in list_of_compound_intensities["Metabolite"].unique():
            intensities = list_of_compound_intensities[
                list_of_compound_intensities["Metabolite"] == met_name
            ]
            if len(intensities) >= min_reps:
                mean = np.mean(intensities["Intensities"])
                var = np.var(intensities["Intensities"])
                rsd = np.sqrt(var) / mean
                stats_all_dict[cnt] = {
                    "replicate_group_name": replicate_group,
                    "Metabolite": met_name,
                    "Formula": list(
                        list_of_compound_intensities[
                            list_of_compound_intensities["Metabolite"]
                            == met_name
                        ]["Formula"]
                    )[0],
                    "Mean": mean,
                    "Variance": var,
                    "RSD": rsd,
                }
                cnt = cnt + 1
    return pd.DataFrame.from_dict(stats_all_dict, "index")


def print_formula(elements):
    return ''.join([f'{k}{int(v)}' for k, v in elements.items()])


def zero_charge(metabolite):
    formula = Formula(metabolite.formula)
    if metabolite.charge != 0:
        if 'desulfurated' in metabolite.name:
            formula.elements['S'] += metabolite.charge / 2
            formula = Formula(print_formula(formula.elements))
        else:
            if 'H' in formula.elements and not np.isnan(metabolite.charge):
                formula.elements['H'] += -metabolite.charge
                formula = Formula(print_formula(formula.elements))
    return formula


def is_valid(metabolite):
    if not metabolite.formula:
        return False
    for string in INVALID_FORMULA_STR:
        if string in metabolite.formula:
            return False


def make_struct(formulas, ids):
    masses = [0]*len(formulas)
    df_formulas = pd.DataFrame(columns=['mass', 'formula', 'id'])
    df_formulas['mass'] = masses
    df_formulas['formula'] = formulas
    df_formulas['id'] = ids
    return df_formulas


def make_mapping(df_formulas):
    map_id = {}
    for formula, df in df_formulas.groupby('formula'):
        map_id[formula] = '\t'.join(df['id'])
    df_mapping = df_formulas.drop_duplicates(subset=['formula'])
    df_mapping['id'] = df_mapping['formula'].map(map_id)
    return df_mapping


def store_struct(df_formulas, name, dirpath):
    df_struct_mapping = pd.DataFrame(columns=['id', 'formula', 'smiles', 'inchi'])
    df_struct_mapping['formula'] = df_formulas['formula']
    df_struct_mapping['id'] = df_formulas['id']
    df_struct_mapping['smiles'] = 'smiles'
    df_struct_mapping['inchi'] = 'inchi'
    df_struct_mapping.to_csv(f'{dirpath}/{name}_struct.tsv', sep='\t', header=None, index=None)


def store_mapping(df_mapping, name, dirpath):
    filename_mapping_csv = f'{dirpath}/{name}_mapping_csv.tsv'
    filename_mapping = f'{dirpath}/{name}_mapping.tsv'
    df_mapping.to_csv(filename_mapping_csv, sep='\t', header=None, index=None, quoting=0)
    with open(filename_mapping_csv) as f:
        with open(filename_mapping, 'w') as fm:
            fm.write(f'database_name	{name}\ndatabase_version	{name}\n')
            for l in f:
                fm.write(l.replace('"', ''))


def create_database(metabolites, name, dirpath):
    formulas = []
    ids = []
    for m in metabolites:
        if not is_valid(m):
            continue
        formula = zero_charge(m)
        formulas.append(str(formula))
        ids.append(m.id)

    df_formulas = make_struct(formulas, ids)
    df_mapping = make_mapping(df_formulas)

    store_struct(df_formulas, name, dirpath)
    store_mapping(df_mapping, name, dirpath)