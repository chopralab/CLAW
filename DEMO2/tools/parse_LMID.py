
import requests
# import requests
import json
import pandas as pd
import re





def process_data(json_text):
    data = json.loads(json_text)
    # print(len(data))
    results = []
    try:
        if "Row1" in data:


            for row_data_index in range(len(data)):

                string  = "Row"+str(row_data_index+1)

                row_data = data[string]

                result = {
                    "abbreviation": row_data.get("abbrev"),
                    "core_class": row_data.get("core"),
                    "main_class": row_data.get("main_class"),
                    "sub_class": row_data.get("sub_class"),
                    "lm_id": row_data.get("lm_id"),
                    "name": row_data.get("name"),
                    "sys_name": row_data.get("sys_name"),
                    "inchi_key": row_data.get("inchi_key"),
                    "formula": row_data.get("formula"),
                    "exact_mass": row_data.get("exactmass"),
                    "smiles": row_data.get("smiles"),
                    "pubchem_cid": row_data.get("pubchem_cid"),
                    "chebi_id": row_data.get("chebi_id"),
                }

                results.append(result)

        else:
            result = {
                "abbreviation": data.get("abbrev"),
                "core_class": data.get("core"),
                "main_class": data.get("main_class"),
                "sub_class": data.get("sub_class"),
                "lm_id": data.get("lm_id"),
                "name": data.get("name"),
                "sys_name": data.get("sys_name"),
                "inchi_key": data.get("inchi_key"),
                "formula": data.get("formula"),
                "exact_mass": data.get("exactmass"),
                "smiles": data.get("smiles"),
                "pubchem_cid": data.get("pubchem_cid"),
                "chebi_id": data.get("chebi_id"),
            }

            results.append(result)

    #         results.append(result)

        return results
    except:
            return "Empty List"
        
def get_compound_info(compound_name):
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    compound_url = f"{base_url}/compound/name/{compound_name}/property/MolecularFormula,MolecularWeight,IUPACName,InChIKey/JSON"
    response = requests.get(compound_url)
    
    if response.status_code == 200:
        data = json.loads(response.text)
        compound_info = data["PropertyTable"]["Properties"][0]
        return compound_info
    else:
        return None

def get_lipid_info_from_inchi(inchikey):
    base_url = "http://www.lipidmaps.org/rest"
    lipid_url = f"{base_url}/compound/inchi_key/{inchikey}/all/json"

    response = requests.get(lipid_url)
    


    if response.status_code == 200:
        return response.text.strip()  # remove leading/trailing whitespace
    else:
        return None

#Use the function:


    
def get_lipid_id_from_name(lipid_name):
    url = f"https://www.lipidmaps.org/rest/compound/abbrev/{lipid_name}/all/json"

    response = requests.get(url)

    if response.status_code == 200:
        json_string=  response.text.strip()  # remove leading/trailing whitespace
        return process_data(json_string)
    else:
        info = get_compound_info(lipid_name)
        if info:
            inchi = info['InChIKey']
            json_string = get_lipid_info_from_inchi(inchi)

    
            return process_data(json_string)
    
        else:
            return "Empty List"
    
    



def parse_lipid_maps_info(Lipid_info):
    LM_IDs = []
    LM_names = []
    if Lipid_info == ["Empty List"]:
        LM_IDs.append("Not Found")
        LM_names.append("Not Found")
    else:
        for i in Lipid_info:
            LM_IDs.append(i["lm_id"])
            LM_names.append(i["name"])
            
    return LM_IDs, LM_names


def get_lipid_info(input_string):
    # Split the string into individual lipids
    lipid_names = input_string.split('|')
    info_list = []
    for i in lipid_names:
        i = i.strip()
        temp_info = get_lipid_id_from_name(i)
        if temp_info == "Empty List":
            continue
        else:
            
            info_list = info_list + temp_info
            
    if len(info_list) == 0:
            info_list.append("Empty List")
        
    

    return parse_lipid_maps_info(info_list)



# x,y = get_lipid_info('DG(26:0)')
try:
    df = pd.read_csv("lipid_database/Custom_MRMs.csv")

    # apply function to lipid column
    lipid_info = df['Lipid'].apply(get_lipid_info)

    # split the returned tuple into two separate Series
    IDs, Names = zip(*lipid_info)

    # If the IDs and Names are lists and you want to convert them to strings
    IDs = [', '.join(map(str, x)) for x in IDs]
    Names = [', '.join(map(str, x)) for x in Names]

    # Now, add these lists to the DataFrame using .assign()
    df = df.assign(IDs=IDs, Names=Names)
    df.fillna("N/A", inplace=True)

    df.to_csv("lipid_database/Custom_MRMs.csv")
except:
    pass

