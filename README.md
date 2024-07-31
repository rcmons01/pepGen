# pepGen
A stand-alone script for generation of combinatorial peptoid (NOT PEPTIDE) libraries for virtual screening.
![image](https://github.com/user-attachments/assets/0164a960-b94c-4e54-af79-b3e432535ad0)


pepGen is a python script for generating n-mer peptoid libraries of any composition. [Genscript](https://www.genscript.com/peptoid_synthesis.html) currently has 20 peptoid monomers listed and available for synthesis and so the current library ([genscript_SMILES.txt](https://github.com/rcmons01/pepGen/blob/main/genscript_SMILES.txt)) reflects these residues. However, this list can be expanded or trimmed down as needed. The requirements for creating and incorporating new residues is that there is a Bromine (Br) group attached at the terminal nitrogen (N) group and that there is an amide group at the opposite end. This script can generate either all combinations of an n-mer or a single defined peptoid molecule. It can generate as SDF or MOL2. 

This script heavily utilizes libraries from [stk](https://github.com/thestk/stk), [rdkit](https://github.com/rdkit/rdkit), [openbabel](https://github.com/openbabel/openbabel), and utilized ChatGPT for debugging. 

## Conda install

```sh
conda create -n peptoid_env python=3.9
conda activate peptoid_env
conda install -c conda-forge rdkit
conda install -c conda-forge openbabel
pip install stk
```

## PIP install

### Linux
```sh
python -m venv peptoid_env
source peptoid_env/bin/activate
pip install rdkit-pypi
pip install openbabel
pip install stk
```

### Windows
```sh
python -m venv peptoid_env
source peptoid_env\Scripts\activate
pip install rdkit-pypi
pip install openbabel
pip install stk
```

## Examples

The program will either create every combination of length n or can be used to generate a specific sequence. the "--format" flag can be enabled to switch to .sdf. 

### Create all combinations of 3-mer peptoids (20^3 = 8,000 molecules)
```py
python pepgen.py --smiles_file genscript_SMILES.txt --output_dir pepgen_output --mode combinations --length 3
```

### Create a distinct 10-mer peptoid
```py
python pepgen.py --smiles_file genscript_SMILES.txt --output_dir pepgen_output --mode single --sequence "Ndip_NVal_NVal_Nmba_NIle_NLeu_Nffa_Ntbu_Nffa_Nffa"
```

---

***KEEP IN MIND THAT IF ALL 20 RESIDUES ARE USED, THE PEPTOID LIBRARY WILL SCALE AS 20^n!!! So a 5-mer library will be = 3.2 million peptoids***
