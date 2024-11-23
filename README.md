# scaffold-structures-rdkit
Search  by SMILES of interest and highlight Scaffolds (Rs)

**Scaffold-Structures-RDKit** is a Python script designed for scaffold and substructure analysis in molecular datasets using **RDKit**. It processes an Excel sheet containing SMILES strings, converting them into molecular objects for visualization. The user can input a SMILES query to search for specific substructures, highlighting matches in the molecules and visualizing the results in a grid. Additionally, the script identifies and highlights dissimilar scaffold regions (Rs) for detailed analysis. The results, including highlighted images and related data, are saved as output files. This tool is ideal for cheminformatics workflows involving scaffold identification and visualization.


1. Import Libraries
```
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

from io import StringIO
```

2. Input Excel file
```
input_xlsx = input("Enter the directory of the Excel sheet: ")
data_df = pd.read_excel(input_xlsx)
index_flag = input("Does the Excel sheet have an index? Type 'True' or 'False': ").lower() == 'false'
```

3. Convert excel to CSV file
```
csv_data = data_df.to_csv(index=index_flag)
csv_file = StringIO(csv_data)
df = pd.read_csv(csv_file)
```

4. convert SMILES to MOLS for visualization
```
df['mols'] = [Chem.MolFromSmiles(smiles) for smiles in df['SMILES']]
```

5. Search by SMILES of interest and Highlight Rs
```
smilesofinterest = input('write SMILES of interest, ex: O=c1[nH]cnc2c1ccn2[C@H]1CCCO1')
molofinterest = Chem.MolFromSmiles(smilesofinterest)
mols_to_display = []
for mol_smi, code in zip(df["SMILES"], df["Code"]):
    mol = Chem.MolFromSmiles(mol_smi)
    if mol is None:
        continue
    hit_ats = list(mol.GetSubstructMatch(molofinterest))
    if hit_ats:
        mols_to_display.append((mol, hit_ats, code))

mols = [item[0] for item in mols_to_display]
hit_ats_list = [item[1] for item in mols_to_display]
codes = [item[2] for item in mols_to_display]
img1 = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(300, 300), 
                            legends=[name for name in df["Code"]], 
                            highlightAtomLists=hit_ats_list ,returnPNG=False)

print("Similarity search by SMILES of interest")
img1.save(input("To save agents of SMILES of interest: /directory/../name.png"))
display(img1)


df_R = pd.DataFrame({"Code": codes, "mols": mols})
df_R = df_R.merge(df[['Code', 'MRC-5', 'T. cruzi', 'L. inf', 'T. bruc', 'T. rhod','PMM']], 
                       on='Code', how='left')

# Higlight Rs (Disimilarity)
def get_atoms_to_highlight(mol_scaffold_of_interest, mol_compound_of_interest):
    atom_list = [i for i in range(mol_compound_of_interest.GetNumAtoms())]
    grouped_common_atoms = mol_compound_of_interest.GetSubstructMatches(mol_scaffold_of_interest)
    flattened_grouped_common_atoms = list(sum(grouped_common_atoms, ()))
    atoms_to_highlight = list(set(atom_list) - set(flattened_grouped_common_atoms))
    return atoms_to_highlight

hl_atom_lists = []

for mol in df_R["mols"]:
    highlight = get_atoms_to_highlight(molofinterest, mol)
    mol.__sssAtoms = highlight
    hl_atom_lists.append(highlight)

img = Draw.MolsToGridImage(df_R["mols"], 
                    legends=[name for name in df["Code"]], 
                           highlightAtomLists=hl_atom_lists, 
                           molsPerRow=4, maxMols= 1000,
                           subImgSize=(300, 300), useSVG=False,returnPNG=False)
print("Highlighting Rs")
img.save(input("To save Rs image: /directory/.../name.png"))
display(img)
display(df_R)
df_0.to_csv(input("To save csv file for the SMILES of interest: /directory/.../name.png"), index=False)
```
