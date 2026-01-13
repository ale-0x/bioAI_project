<!-- # bioAI_project

Code written for the exam 'Bio-inspired AI' 

We present a genetic algorithm capable of evolving a small peptide sequence to bind a pre-defined pocket in a protein structure.

To show this, we employed the resolved crystal structure of F HIV protease complexed with TL-3 inhibitor as an example ([2P3D](https://www.rcsb.org/structure/2P3D)). 
Our effort was to replace the TL-3 inhibitor with a custom made peptide sequence evolved with a genetic algorithm and evaluate its match to the shape of the pocket.

To achieve this result, the genetic algorithm is designed as such:

## Population
The individuals that make up the population are candidate peptide sequences. Their generation is taken care of with `ga_problem.peptide_generator()`, that will create a random string to begin the computation.

## Evaluation
Every candidate will be passed as a molecule to Vina, a physics-based docking software that will predict the complex 3D structure and binding pose of the protein with the candidate sequence in a _specified binding box_.
The binding will be assessed in terms of binding energy of the final complex, the baseline minimization objective; further constrains on the solution space are imposed by a candidate hydrophobicity measure, with the algorithm penalizing hydrophobicity in candidates to favor drug absorption. This criterion is included in the final objective function that looks like:

$$\text{binding energy} + \alpha \times \text{hydrophobicity}$$

## Mutation and Selection
After the candidate population is evaluated, each individual will undergo crossover and/or point mutation to introduce variability. The crossover operator will produce a hybrid sequence between two candidates with a random breakpoint, while the point mutator will scan the candidate string and for each position employ blosum weights to produce another likely aminoacid to replace the original one. Both processes are governed by user-defined mutation likelihoods.

## Final assessment
The final candidate is returned along with a structural prediction and is available for further inspection. -->

# Evolutionary Peptide-Protein Docking

This project explores the application of **Evolutionary Algorithms (EA)** for the *de novo* design of peptide ligands targeting specific protein receptors. The problem of finding a peptide sequence with high binding affinity is combinatorial and computationally expensive. We propose a **Genetic Algorithm (GA)** implemented using the `inspyred` library, coupled with `AutoDock Vina` for fitness evaluation. The pipeline integrates `OpenBabel` for chemical file conversion and customized genetic operators (mutation, crossover) designed for amino acid sequences. Experimental results demonstrate the algorithm's ability to evolve peptides with progressively lower binding energies (higher affinity) over generations.

## 🧬 Project Overview

This tool automates the search for optimal peptide sequences that bind to a specific pocket of a target protein. Instead of screening a static library of compounds, it **evolves** the solution using biological principles.

### Key Features
* **Genetic Algorithm**: Uses `inspyred` to manage population, selection, and evolution.
* **Physics-Based Scoring**: Utilizes **AutoDock Vina** to evaluate the binding energy (kcal/mol) of generated peptides.
* **Bio-Inspired Operators**:
    * **BLOSUM62 Mutation**: Amino acid substitutions are weighted based on evolutionary likelihood (using the BLOSUM62 matrix), not just random chance.
    * **Single-Point Crossover**: Recombines parent sequences to explore new solution spaces.
* **Parallel Processing**: Evaluates multiple individuals concurrently using Python's `multiprocessing` to speed up docking.
* **Multi-Objective Optimization**: Optimizes for binding energy while penalizing high hydrophobicity to improve solubility and drug-likeness.

---

## ⚙️ Installation & Prerequisites

### 1. System Dependencies
You need the following tools installed and accessible in your system path:
* **AutoDock Vina** (executable named `vina`)
* **OpenBabel** (for chemical conversion)

### 2. Python Environment
We recommend using Conda to manage dependencies. A complete environment file is provided. To replicate the exact environment, run the following commands in order:

```bash
# Create the environment from file
conda env create -f environment.yml

# Activate the environment
conda activate bioai
```

Or manually:


```bash
# 1. Create and activate the environment
conda create -n bioai python=3.9
conda activate bioai

# 2. Install chemical and analysis tools from conda-forge
conda install -c conda-forge openbabel rdkit pymol-open-source biopython

# 3. Install AutoDock Vina from bioconda
conda install -c bioconda autodock-vina

# 4. Install evolutionary library via pip
pip install inspyred
```

If you prefer pip, the core requirements are:
* `inspyred`
* `numpy`
* `matplotlib`
* `biopython` (for BLOSUM matrices)
* `rdkit` (for advanced chemical handling)

## 🚀 Usage
The main entry point is `main.py`. You need to provide the target receptor name and the prepared PDBQT file.

### Basic Command
```bash
python main.py receptor_name path/to/receptor.pdbqt
```
### Example (using provided data)
To run the optimization on the HIV Protease (2P3D) structure:

```bash
python3 -u main.py                                              \
    2P3D_no_ligand                                              \
    $RECEPTOR                                                   \
    --job_id                $SLURM_JOB_ID                       \
    --cpus                  $SLURM_CPUS_PER_TASK                \
    --peptide_length        6                                   \
    --population_size       5                                   \
    --generations           8                                   \
    --initial_mutation_rate 0.30                                \
    --final_mutation_rate   0.05                                \
    --hydrophobicity_weight 0.10                                \
    --temp_dir_base         $TMP_BASE                           \
    --output                ../results/result                   \
    --center_x=8.084                                            \
    --center_y=-13.829                                          \
    --center_z=-0.140                                           \
    --size_x                32                                  \
    --size_y                32                                  \
    --size_z                32                                  \
    --exhaustiveness        4                                   \
    --vina_exe_path         vina                                \
    --no_delete
```

### Command Line Arguments

The script accepts several arguments to customize the simulation. They are categorized below by function.

#### 1. Required Arguments

| Argument | Description |
| --- | --- |
| `receptor_name` | Identifier code for the target protein (e.g., `2P3D`). Used for naming output files. |
| `receptor_file` | Relative path to the prepared receptor file in **.pdbqt** format (must include charges and hydrogens). |

#### 2. Genetic Algorithm & General Settings

| Argument | Flag | Default | Description |
| --- | --- | --- | --- |
| `--job_id` | `-j` | `local_test` | Unique identifier for the job (useful for HPC/Slurm clusters). |
| `--cpus` | `-c` | *(from constants)* | Number of CPU cores to use for parallel docking. |
| `--peptide_length` | `-l` | `10` | Length of the evolved peptide sequence (number of amino acids). |
| `--population_size` | `-n` | `40` | Number of individuals (peptides) per generation. |
| `--generations` | `-g` | `30` | Maximum number of generations to run. |
| `--initial_mutation_rate` |  | `0.3` | Probability of mutation at the start of the simulation. |
| `--final_mutation_rate` |  | `0.01` | Probability of mutation at the end (simulated annealing approach). |
| `--hydrophobicity_weight` |  | `0.1` | Weight factor () used to penalize hydrophobic peptides in the fitness function. |
| `--output` | `-o` | `result.txt` | Name of the final output folder/file. |
| `--temp_dir_base` |  | `../resources/tmp` | Base directory for temporary files generated during docking. |

#### 3. AutoDock Vina Grid Configuration

Define the search space (Grid Box) where Vina will look for binding.

| Argument | Flag | Default | Description |
| --- | --- | --- | --- |
| `--center_x` | `-x` | *(from constants)* | X coordinate of the grid box center. |
| `--center_y` | `-y` | *(from constants)* | Y coordinate of the grid box center. |
| `--center_z` | `-z` | *(from constants)* | Z coordinate of the grid box center. |
| `--size_x` | `-X` | `32` | Size of the search box along the X axis (Ångström). |
| `--size_y` | `-Y` | `32` | Size of the search box along the Y axis (Ångström). |
| `--size_z` | `-Z` | `32` | Size of the search box along the Z axis (Ångström). |
| `--exhaustiveness` | `-e` | `8` | Search exhaustiveness (higher = more precise but slower). |
| `--vina_exe_path` |  | `vina` | Path to the AutoDock Vina executable. |

#### 4. Utility Flags

| Flag | Description |
| --- | --- |
| `--verbose` (`-v`) | Enable verbose logging to the terminal. |
| `--no_delete` | If set, temporary files (PDB/PDBQT candidates) are **not** deleted after execution (useful for debugging). |
| `--help` (`-h`) | Show help for using the command |
| `--version` | Show version information and exit. |

## 🧪 Methodology

### The Evolutionary Cycle
1. **Initialization**: A population of random peptide sequences (strings of amino acids) is generated
2. **Evaluation**:
    - The peptide string is converted into a 3D PDB structure (linear chain).
    - Converted to PDBQT format (OpenBabel)
    - **Docking**: AutoDock Vina attempts to fit the peptide into the target pocket
    - **Fitness Calculation**: The best binding affinity (lowest energy) is returned
3. **Selection**: Tournament selection picks the best individuals
4. **Reproduction**:
    - **Crossover**: Mixes parts of two sequences
    - **Mutation**: Changes amino acids based on BLOSUM62 probability (preserving chemical properties)
5. **Termination**: Stops after `MAX_GENERATIONS`.
### Fitness Function
The algorithm minimizes the objective function:$$ F(x) = E_{binding} + (\alpha \times H_{index}) $$
Where:
- $E_{binding}$: Vina docking score (negative kcal/mol is better)
- $H_{index}$: Average hydrophobicity of the peptide (calculated at pH 7)
- $\alpha$: Weight factor (`HYDROPHOBICITY_WEIGHT`) to penalize overly hydrophobic peptides (which are often insoluble).

## 📊 Results & Visualization
At the end of the run, the `plots/` folder will contain:
1. **Generation Plot**: Trends of Best, Average, and Worst fitness over time.
2. **Landscape Plot**: Scatter plot of Binding Energy vs. Hydrophobicity, showing the population's evolution trajectory.

## 👥 Authors

| Name | Tag | Role |
| --- | --- | --- |
| [Alex Callegaro](https://github.com/ale-0x) | @ale-0x | Master's student in Quantitative and Computational Biology, University of Trento |
| [Clemente Calabrese](https://github.com/cl3mente) | @cl3mente | Master's student in Quantitative and Computational Biology, University of Trento |
| [Federico Cavasin](https://github.com/Clausola) | @Clausola | Master's student in Quantitative and Computational Biology, University of Trento |
| [Giada Silvaggi](https://github.com/giadasil) | @giadasil | Master's student in Quantitative and Computational Biology, University of Trento |
