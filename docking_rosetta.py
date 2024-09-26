import os
import sys
import pymol
import argparse
import subprocess

def get_bcc_charge(mol_path, mol_charge, antechamber_path="~/miniconda3/envs/ambertools/bin/antechamber"):
    """
    This function calculates the BCC charge of a molecule using antechamber.
    """
    if not os.path.exists(mol_path):
        sys.exit("The molecule file does not exist")
    
    mol_name = os.path.basename(mol_path)
    mol_name_noext = os.path.splitext(mol_name)[0]
    mol_dir = os.path.dirname(mol_path)
    mol_ext = os.path.splitext(mol_name)[1]
    
    if mol_ext != ".mol2":
        sys.exit("The molecule file extension must be .mol2")
    
    antechamber_cmd = f"{antechamber_path} -i {mol_path} -fi mol2 -o {mol_dir}/{mol_name_noext}.bcc.mol2 -fo mol2 -c bcc -nc {mol_charge}"
    subprocess.run(antechamber_cmd, shell=True)
    
    return None

def get_rosetta_param(mol_path, mol2genpy_path="$ROSETTA/main/source/scripts/python/public/generic_potential/mol2genparams.py"):
    """
    This function generates the Rosetta parameters for a molecule using mol2genpy.
    """
    if not os.path.exists(mol_path):
        sys.exit(f"The molecule file {mol_path} does not exist")
    
    mol_name = os.path.basename(mol_path)
    mol_name_noext = os.path.splitext(mol_name)[0]
    mol_dir = os.path.dirname(mol_path)
    mol_ext = os.path.splitext(mol_name)[1]
    
    if mol_ext != ".mol2":
        sys.exit("The molecule file extension must be .mol2")
    
    mol2genpy_cmd = f"python {mol2genpy_path} -s {mol_path} --outdir {mol_dir} --prefix LG1"

    subprocess.run(mol2genpy_cmd, shell=True)
    os.system("rm -rf *A") # remove the old files
    os.system("rm -rf *sam") # remove the old files
    return None

def make_complex(pro_pdb, lig_pdb):
    """
    Use PyMOL to make a protein-ligand complex.
    """
    pymol.cmd.load(pro_pdb, "protein")
    pymol.cmd.load(lig_pdb, "ligand")
    pymol.cmd.create("complex", "protein or ligand")
    save_name = os.path.splitext(pro_pdb)[0] + "_complex.pdb"
    pymol.cmd.save(save_name, "complex")
    
    return save_name

def GALigandDock(complex, flags_file, dock_xml, num_process=5, nstruct=1):
    """
    This function docks a ligand to a protein using Rosetta.
    """
    if not os.path.exists(complex):
        sys.exit(f"The complex file {complex} does not exist")
    
    complex_name = os.path.basename(complex)
    complex_name_noext = os.path.splitext(complex_name)[0]
    complex_dir = os.path.dirname(complex)
    complex_ext = os.path.splitext(complex_name)[1]
    
    if complex_ext != ".pdb":
        sys.exit("The complex file extension must be .pdb")
    
    with open(flags_file, "r") as f:
        lines = f.readlines()
    
    with open(flags_file, "w") as f:
        for line in lines:
            if line.startswith("-extra_res_fa "):
                f.write(f"-extra_res_fa {complex_dir}/LG1.params\n")
            elif line.startswith("-out:path:all "):
                f.write(f"-out:path:all {complex_dir}\n")
            elif line.startswith("-nstruct "):
                f.write(f"-nstruct {nstruct}\n")
            else:
                f.write(line)
        f.close()
    # save the flags file to the complex directory
    # flags_file = f"{complex_dir}/flags"

    
    docking_cmd = f"mpirun -np {num_process} rosetta_scripts.mpi.linuxgccrelease -s {complex_dir}/{complex_name_noext}.pdb -parser:protocol {dock_xml} @{flags_file}"
    print("Running the docking command:", docking_cmd)
    subprocess.run(docking_cmd, shell=True)
    
    return None

def main():
    description = (
        "This script docks a ligand to a protein using Rosetta GALigandDock.\n\n"
        "It calculates the BCC charge of the ligand, generates the Rosetta parameters,"
        "makes a protein-ligand complex, and docks the ligand to the protein.\n\n"
        "The ligand must be in mol2 format, and its coordinates must near the protein pocket.\n\n"
        "Author: DingLuo\n"
        "Email: luoding@sust.edu.cn\n"
        "version: 1.0\n"
        "Date: 2024-08-03"
    )
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-p", "--protein", required=True, help="The protein PDB file")
    parser.add_argument("-l", "--ligand", required=True, help="The ligand mol2 file")
    parser.add_argument("-f", "--flags", help="The flags file", default="./input/flags")
    parser.add_argument("-d", "--dock", help="The docking XML file", default="./input/dock.xml")
    parser.add_argument("-n", "--process", type=int, help="The number of processes", default=5)
    parser.add_argument("-s", "--nstruct", type=int, help="The number of structures", default=1)
    parser.add_argument("-v", "--version", action="version", version="%(prog)s 1.0")
    parser.add_argument("-c", "--charge", type=bool, help="Whether to calculate the BCC charge", default=True)
    parser.add_argument("-nc", "--charge_num", type=int, help="The charge number of the ligand", default=0)
    args = parser.parse_args()
    
    inputfile = args.ligand
    ligand_name = os.path.splitext(inputfile)[0]
    dir = os.path.dirname(inputfile)
    if args.charge:
        bcc_output = ligand_name + ".bcc.mol2"
        
        if not os.path.exists(bcc_output):
            print("Calculating BCC charges...")
            get_bcc_charge(args.ligand,args.charge_num)
    # bcc_output = ligand_name + ".bcc.mol2"
    
    # if not os.path.exists(bcc_output):
    #     print("Calculating BCC charges...")
    #     get_bcc_charge(args.ligand)
    
        if not os.path.exists(f"{dir}/LG1.params"):
            print("Generating Rosetta parameters...")
            get_rosetta_param(bcc_output)
    else:
        if not os.path.exists(f"{dir}/LG1.params"):
            print("Generating Rosetta parameters...")
            get_rosetta_param(args.ligand)
    
    save_name = make_complex(args.protein, f"{dir}/LG1_0001.pdb")
    
    GALigandDock(save_name, args.flags, args.dock, args.process)
    print("Docking is done")

if __name__ == "__main__":
    main()
