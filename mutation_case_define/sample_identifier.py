import json
import pandas as pd
from collections import Counter

if __name__ == '__main__':

    cell_line_dbs = pd.read_csv("../dbs/DepMap-2018q3-celllines.csv")
    Counter(cell_line_dbs["Subtype_Disease"].to_list())
    chosen_subtype_keys = ['carcinoma',
                   'adenocarcinoma',
                   'lymphoid_neoplasm',
                   'malignant_melanoma',
                   'glioma',
                   'Basal',
                   'squamous_cell_carcinoma',
                   'Luminal',
                   'hepatocellular_carcinoma',
                   'haematopoietic_neoplasm',
                   'neuroblastoma']

    CCLE_mutation_dbs_non_filtered = pd.read_csv("../dbs/CCLE_mutations.csv")
    CCLE_mutation_dbs = CCLE_mutation_dbs_non_filtered[CCLE_mutation_dbs_non_filtered["Variant_Classification"] == "Missense_Mutation"]

    just_position_conversion_ls = [str(row["Protein_Change"])[0:-1] for index, row in CCLE_mutation_dbs.iterrows()]
    CCLE_mutation_dbs["Protein_Change_Position"] = just_position_conversion_ls

    #print(len(CCLE_mutation_dbs["Protein_Change_Position"].to_list()),len(CCLE_mutation_dbs["Protein_Change"].to_list()) )

    subtype_disease_case_defined_dict = dict()
    for subtype_disease in chosen_subtype_keys:
        cases_of_chosen_subtype = cell_line_dbs[cell_line_dbs["Subtype_Disease"] == subtype_disease]["Broad_ID"].to_list()

        subtype_disease_case_defined_dict[subtype_disease], subtype_interested_mutation_list, case_mutation_defined_dict = cases_of_chosen_subtype, list(), dict()
        subtype_mutation_case_defined_dict_to_json = dict()

        for case in cases_of_chosen_subtype:
            mutation_dbs_of_chosen_case = (CCLE_mutation_dbs[CCLE_mutation_dbs["DepMap_ID"] == case])

            Hugo_Symbols, Protein_Changes = mutation_dbs_of_chosen_case["Hugo_Symbol"].to_list(),\
                                            mutation_dbs_of_chosen_case["Protein_Change_Position"].to_list()

            for i in range(len(Hugo_Symbols)):
                try:
                    subtype_interested_mutation_list.append(str(Hugo_Symbols[i]) + "_" + str(Protein_Changes[i]))
                except:
                    pass

        for mutation_HugoSym_Position in list(set(subtype_interested_mutation_list)):
            Hugo_Sym, Position = mutation_HugoSym_Position.split("_")[0], mutation_HugoSym_Position.split("_")[1]

            objected_case_collection_list = CCLE_mutation_dbs[(CCLE_mutation_dbs["Hugo_Symbol"] == Hugo_Sym) & (CCLE_mutation_dbs["Protein_Change_Position"] == Position)]["DepMap_ID"].to_list()

            subtype_mutation_case_defined_dict_to_json[str(Hugo_Sym)+"_"+str(Position).split("p.")[1]] = objected_case_collection_list

        print(subtype_disease)
        subtype_disease_json_mutation_profile_file = open("../rsults/"+subtype_disease+".json","r")
        json.dump(subtype_mutation_case_defined_dict_to_json, subtype_disease_json_mutation_profile_file)
        subtype_disease_json_mutation_profile_file.close()

