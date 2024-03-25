#!/usr/bin/env python
# coding: utf-8

# In[13]:


import os
import pandas as pd
import sys
from optparse import OptionParser

from rdkit import Chem
from rdkit.Chem import Draw,AllChem
from rdkit.Chem.Draw import IPythonConsole
from multiprocessing import Pool


from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
import py3Dmol
import argparse
from ipywidgets import interact, interactive, fixed # For interactive display of conformers
from espsim import EmbedAlignConstrainedScore, EmbedAlignScore, ConstrainedEmbedMultipleConfs, GetEspSim


"The reference file is added in the code. "




# In[10]:


def Mol2MolSupplier (file=None,sanitize=True):
    mols=[]
    with open(file, 'r') as f:
        line =f.readline()
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol = []
                mol.append(line)
                line = f.readline()
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        mol.append(line)
                        break
                mol[-1] = mol[-1].rstrip() # removes blank line at file end
                block = ",".join(mol).replace(',','')
                m=Chem.MolFromMol2Block(block,sanitize=sanitize)
            mols.append(m)
    return(mols)


# In[14]:


filePath = "ref.mol2"
database=Mol2MolSupplier(filePath,sanitize=True)

prbMol=Chem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(database[0])))
AllChem.EmbedMolecule(prbMol,AllChem.ETKDG())
AllChem.UFFOptimizeMolecule(prbMol)


# In[ ]:


def run(query_smi):

    mol_id = query_smi[1]
    smiles = query_smi[0]
    
    output = []
    all_prbmols = []
    all_refmols = []
   
    try:

        refMols = [Chem.AddHs(Chem.MolFromSmiles(smiles))]        
        simShape,simEsp=EmbedAlignScore(prbMol,refMols, partialCharges='mmff',metric='tanimoto',renormalize=True)
        all_data = []
        for j in range(10):
            for k in range(10):
                all_data.append([j,k,GetEspSim(prbMol,refMols[0],j,k)])
                                
        i_ = int(pd.DataFrame(all_data, columns=['j','k','EspSim']).sort_values(by=['EspSim'],ascending=False).reset_index(drop=True).iloc[0]["j"])
        j_ = int(pd.DataFrame(all_data, columns=['j','k','EspSim']).sort_values(by=['EspSim'], ascending=False).reset_index(drop=True).iloc[0]["k"])
        
 
        prb = Chem.MolFromMolBlock(Chem.MolToMolBlock(prbMol, confId=i_))
#         prb.SetProp("_Name", mol_id)
        
        ref = Chem.MolFromMolBlock(Chem.MolToMolBlock(refMols[0], confId=j_))
#         ref.SetProp("_Name", mol_id)
#         ref.SetProp("simShape", str(simShape[0]))
#         ref.SetProp("simEsp", str(simEsp[0]))
#         ref.SetProp("combi",str(simShape[0]*simEsp[0]))
        
        
        output.append([prb,ref,simShape[0],simEsp[0],simShape[0]*simEsp[0],mol_id])
 
    except:
        output.append(["NA","NA"])

    
    return output
    
    


# In[12]:


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description= "python Parallel_shape_similarity referencesdf ncpus refout_sdf queryout_sdf")
    parser.add_argument("referencesdf", help="reference sdf")
    parser.add_argument("ncpus", help="Input the number of cpus")
    parser.add_argument("refout_sdf", help="This is the output file of all the reference sdf files. Espsim moves the reference sdf files")
    parser.add_argument("queryout_sdf", help="This is the output file of all the query sdf files")
    args = parser.parse_args()
    
#     query = args.querysdf
    reference = args.referencesdf
    ref_sdf = args.refout_sdf
    query_sdf = args.queryout_sdf
    out = []
    ncpus = int(args.ncpus)
    
#     refmols = []
#     querymols = []
    
    suppl = Chem.SDMolSupplier(reference)  
    refMols = [(Chem.MolToSmiles(x),x.GetProp("idnumber")) for x in suppl if x is not None]
#    refMols = [(Chem.MolToSmiles(x),x.GetProp("Catalog ID")) for x in suppl if x is not None]

    p = Pool(ncpus)
    
    for x in p.map(run, refMols):
        for y in x:
            if(y[0] != "NA"):
                out.append(y)

df = pd.DataFrame(out)
df.columns =['ref','query','shapesim', 'espsim', 'combi','id']
df.sort_values(by=['combi'], ascending=False, inplace=True)

refmols = []
querymols = []

for _,x in df.iterrows():
    x['ref'].SetProp("_Name", x['id'])
    x['query'].SetProp("_Name", x['id'])
    x['query'].SetProp("Shape_Tanimoto", "{:0.2f}".format(x['shapesim']))
    x['query'].SetProp("Electrostatic_Tanimoto", "{:0.2f}".format(x['espsim']))
    x['query'].SetProp("Shape_Tanim_X_Elec_Tanim", "{:0.2f}".format(x['combi']))
    refmols.append(x['ref'])
    querymols.append(x['query'])
    
 # Finally write sorted sdf files

w = Chem.SDWriter(ref_sdf)
for m in refmols:
    w.write(m)

w = Chem.SDWriter(query_sdf)
for m in querymols:
    w.write(m)
