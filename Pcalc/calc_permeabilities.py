import os
import glob
import re
import numpy as np

def get_tprs(base_dir, zzmin, zzmax, zzstep, output_file):
    """Generates a list of .tpr files and saves to output_file."""
    with open(output_file, "w") as f:
        for u in range(int(zzmin * 10), int(zzmax * 10) + 1, int(zzstep * 10)):
            umb = f"{u / 10:.1f}"
            tpr_path = os.path.join(base_dir, f"umbrella_d_{umb}/01_prod_r1/001.tpr")
            gz_path = tpr_path + ".gz"
            
            if os.path.exists(gz_path):
                os.system(f"gunzip {gz_path}")
            
            if os.path.exists(tpr_path):
                f.write(f"{tpr_path}\n")

def get_pulls_all(base_dir, zzmin, zzmax, zzstep, tinit, tend, output_dir, filesf_path):
    """Extracts data from .xvg files and compiles file paths."""
    os.makedirs(output_dir, exist_ok=True)
    
    with open(filesf_path, "w") as filesf:
        for u in range(int(zzmin * 10), int(zzmax * 10) + 1, int(zzstep * 10)):
            umb = f"{u / 10:.1f}"
            for var in ["f", "x"]:
                umb_xvg_path = os.path.join(output_dir, f"umb_{umb}{var}.xvg")
                with open(umb_xvg_path, "w") as xvg_out:
                    for r in [1, 2, 3]:
                        xvg_files = glob.glob(os.path.join(base_dir, f"umbrella_d_{umb}/01_prod_r{r}", f"???{var}.xvg"))
                        for xvg in xvg_files:
                            with open(xvg, "r") as infile:
                                n = 0
                                for line in infile:
                                    parts = line.split()
                                    if len(parts) >= 2:
                                        try:
                                            time = float(parts[0])
                                            value = float(parts[1])
                                            if tinit < time <= tend:
                                                n += 1
                                                xvg_out.write(f"{n * 10} {value}\n")
                                        except ValueError:
                                            continue
                if var == "f":
                    filesf.write(f"{umb_xvg_path}\n")

def get_pulls_rep(base_dir, zzmin, zzmax, zzstep, tinit, tend, output_dir, filesf_path, rep):
    """Extracts data from .xvg files and compiles file paths for specified reps."""
    os.makedirs(output_dir, exist_ok=True)
    
    with open(filesf_path, "w") as filesf:
        for u in range(int(zzmin * 10), int(zzmax * 10) + 1, int(zzstep * 10)):
            umb = f"{u / 10:.1f}"
            for var in ["f", "x"]:
                umb_xvg_path = os.path.join(output_dir, f"umb_{umb}{var}.xvg")
                with open(umb_xvg_path, "w") as xvg_out:
                    xvg_files = glob.glob(os.path.join(base_dir, f"umbrella_d_{umb}/01_prod_r{rep}", f"???{var}.xvg"))
                        
                    for xvg in xvg_files:
                        with open(xvg, "r") as infile:
                            n = 0
                            for line in infile:
                                parts = line.split()
                                if len(parts) >= 2:
                                    try:
                                        time = float(parts[0])
                                        value = float(parts[1])
                                        if tinit < time <= tend:
                                            n += 1
                                            xvg_out.write(f"{n * 10} {value}\n")
                                    except ValueError:
                                        continue
                    if var == "f":
                        filesf.write(f"{umb_xvg_path}\n")

def g_wham(gromacs_path, output_dir, largebins, zprof0, bs):
    """Runs the WHAM analysis."""
    os.system(f"rm -f {output_dir}/bsResult.xvg")
    wham_command = (
        f"{gromacs_path} wham -it {output_dir}/tpr-files.dat -if {output_dir}/filesf.dat "
        f"-o {output_dir}/profilef.xvg -hist {output_dir}/histf.xvg "
        f"-unit kCal -tol 1e-6 -bins {largebins} -zprof0 {zprof0} "
        f"-min -0.01 -max 3.61 -nBootstrap {bs} -xvg none > {output_dir}/aux-log 2>&1"
    )
    os.system(wham_command)
    os.system(f"awk 'NR%10==1{{print $0}}' {output_dir}/bsResult.xvg > {output_dir}/aux-bsResult")

def isdm_all(script_path, output_dir, method):
    """Runs the ISDM permeability calculation."""
    os.system(f"{script_path} > {output_dir}/Permeability_{method}")
    #os.system(f"rm -f {output_dir}/aux-bsResult *#")

def isdm_rep(script_path, output_dir, method, rep):
    """Runs the ISDM permeability calculation."""
    os.system(f"{script_path} > {output_dir}/Perm_single_{rep}_{method}")
    #os.system(f"rm -f {output_dir}/aux-bsResult *#")

def extract_value(filename):
    """Extracts the numerical permeability value from the file."""
    with open(filename, 'r') as file:
        for line in file:
            match = re.search(r"#Peff = ([\d\.]+)", line)
            if match:
                return float(match.group(1))
    return None 

def jackknife_error(replicates):
    """Computes the Jackknife error and standard error (SE)."""
    n = len(replicates)
    mean_total = np.mean(replicates)
    pseudo_values = [n * mean_total - (n - 1) * x for x in replicates]
    
    jackknife_var = np.var(pseudo_values, ddof=1) * (n - 1) / n
    jackknife_se = np.sqrt(jackknife_var)
    
    # Traditional standard error (SE)
    standard_se = np.std(replicates, ddof=1) / np.sqrt(n)
    
    return jackknife_var, jackknife_se, standard_se

def extract_permeability(method, output_file):
    """Extracts permeability values for a given method and writes to file."""
    total_value = extract_value(f"Permeability_{method}")
    
    replicates = [extract_value(f"Perm_single_{i}_{method}") for i in range(1, 4)]
    
    jackknife_var, jackknife_se, standard_se = jackknife_error(replicates)
    
    output_file.write(f"{method} {total_value:.6f} {jackknife_se:.6f} {standard_se:.6f} "
                      f"{replicates[0]:.6f} {replicates[1]:.6f} {replicates[2]:.6f}\n")


if __name__ == "__main__":
    import argparse
    
    # Get compound name from user input
    parser = argparse.ArgumentParser(description="Run permeability analysis for a given compound.")
    parser.add_argument("compound_name", help="Name of the compound")
    args = parser.parse_args()

    # Define paths dynamically based on compound name
    compound_name = args.compound_name
    base_dir = f"/home/afortuna/permeability/01_sims/ligands/literature_compounds/{compound_name}"
    output_dir = f"/home/afortuna/permeability/02_analysis/{compound_name}/calc_perm"

    zzmin, zzmax, zzstep = 0.0, 3.6, 0.2
    tinit, tend = 50000, 150000
    tpr_output_file = os.path.join(output_dir, "tpr-files.dat")
    filesf_output_file = os.path.join(output_dir, "filesf.dat")
    gromacs_path = "/gromacs/gromacs-2021.2/bin/gmx"
    script_path = os.path.join(output_dir, "isdm.py")
    largebins = 191  # Example value
    zprof0 = 3.6
    bs = 100
    
    for method in ["no_EP", "EP1"]: 
        method_dir = os.path.join(base_dir, method, "05_umbrellas")
        get_tprs(method_dir, zzmin, zzmax, zzstep, tpr_output_file)
        get_pulls_all(method_dir, zzmin, zzmax, zzstep, tinit, tend, output_dir, filesf_output_file)
        g_wham(gromacs_path, output_dir, largebins, zprof0, bs)
        isdm_all(script_path, output_dir, method)
        
        with open("permeability_results.txt", "w") as output_file:
            output_file.write("Method Total_Value Jackknife_Error SE Rep1 Rep2 Rep3\n")
            extract_permeability("no_EP", output_file)
            extract_permeability("EP1", output_file)

        for rep in range(1, 4): 
            get_pulls_rep(method_dir, zzmin, zzmax, zzstep, tinit, tend, output_dir, filesf_output_file, rep)
            g_wham(gromacs_path, output_dir, largebins, zprof0, bs)
            isdm_rep(script_path, output_dir, method, rep)







        

