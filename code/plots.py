# plots.py

import matplotlib
matplotlib.use("Agg")                   # Use non-interactive backend for file saving without display

import csv
import logomaker         as lm
import matplotlib.pyplot as plt
import numpy             as np
import pandas            as pd
import seaborn           as sns
import warnings

from inspyred.ec           import analysis, Individual
from matplotlib            import cm
from matplotlib.colors     import LinearSegmentedColormap
from matplotlib.patches    import Patch
from matplotlib.ticker     import MaxNLocator
from pathlib               import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import OneHotEncoder


import constants as C

from peptide_operators  import get_hydrophobicity
from utils              import print_error

gradient_map = LinearSegmentedColormap.from_list("custom_gradient", ['#FF0000', "#FFC800"])

def plot_observer_statistics(observer_file_path: str | Path, plot_folder_directory: str | Path) -> None:
    observer_file_path    = Path(observer_file_path)
    plot_folder_directory = Path(plot_folder_directory)

    plot_folder_directory.mkdir(parents = True, exist_ok = True)

    plt.rcParams.update(C.IEEE_PLOT_PARAMS)
    plt.figure(figsize = (C.IEEE_FIGURE_WIDTH, C.IEEE_FIGURE_HEIGHT))
    
    with warnings.catch_warnings(), observer_file_path.open('r') as f:
        warnings.filterwarnings("ignore", message = "FigureCanvasAgg is non-interactive")
        analysis.generation_plot(f)
    
    plt.savefig(plot_folder_directory / f"generation_plot_{observer_file_path.stem.split('_')[-1]}.png")  # generation_plot_<JobID>.png
    plt.savefig(plot_folder_directory / f"generation_plot_{observer_file_path.stem.split('_')[-1]}.pdf")  # generation_plot_<JobID>.pdf
    plt.close()

def plot_energy_vs_hydrophobicity(individuals_file: str | Path, plot_folder_directory: str | Path, hydrophobicity_weight: float = C.HYDROPHOBICITY_WEIGHT) -> None:
    individuals_file      = Path(individuals_file)
    plot_folder_directory = Path(plot_folder_directory)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    job_id = individuals_file.stem.split('_')[-1]

    generations_dict = {}
    with open(individuals_file, 'r') as f:
        for line in f.readlines():
            gen, ind_number, fitness, individual = line.strip().split(',')
            gen = int(gen)
            if gen not in generations_dict:
                generations_dict[gen] = []
            ind = Individual(candidate=individual)
            ind.fitness = float(fitness)
            ind.birthdate = gen
            ind.hydrophobicity = get_hydrophobicity(individual)*hydrophobicity_weight
            ind.energy = ind.fitness - ind.hydrophobicity
            generations_dict[gen].append(ind)

    # Get max generation for color normalization
    max_gen = max(gen for gen in generations_dict.keys()) if generations_dict else 1

    # Plot individuals with color based on age
    for ind in [ind for array in generations_dict.values() for ind in array]:
        energy = ind.energy  # x-axis: energy
        hydrophobicity = ind.hydrophobicity  # y-axis: hydrophobicity calculation
        age = ind.birthdate
        color_value = age / max_gen if max_gen > 0 else 0
        ax.scatter(
            energy, 
            hydrophobicity, 
            color=gradient_map.reversed()(color_value), 
            s=100, 
            alpha = 1,
            )
        ax.text(
            energy, 
            hydrophobicity, 
            ind.candidate, 
            fontsize=8, 
            ha='right', 
            va='bottom'
        )

    # Draw lines connecting individuals from the same generation
    for gen, individuals in sorted(generations_dict.items()):
        if len(individuals) > 1:
            energies = [ind.energy for ind in individuals]
            hydrophobicities = [ind.hydrophobicity for ind in individuals]
            color_value = gen / max_gen if max_gen > 0 else 0
            ax.plot(energies, hydrophobicities, linewidth=2, color=gradient_map.reversed()(color_value))

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=gradient_map.reversed(), 
                            norm=plt.Normalize(vmin=0, vmax=max_gen))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Generation', rotation=270, labelpad=20)

    ax.set_xlabel('Energy (kcal/mol)')
    ax.set_ylabel('Hydrophobicity index a pH 7 (media)')
    ax.set_title(f'Energy vs Hydrophobicity')
    ax.grid(True, alpha=0.3)

    scatter_plot_file = plot_folder_directory / f"energy_hydrophobicity_{job_id}.png"
    plt.savefig(scatter_plot_file, dpi=150, bbox_inches='tight')
    plt.savefig(scatter_plot_file.with_suffix(".pdf"), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Energy vs Hydrophobicity plot saved to {scatter_plot_file}")

def plot_observer_statistics_2(observer_file_path: str | Path, plot_folder_directory: str | Path) -> None:
    """
    Generates an alternative/refined plot for fitness statistics.

    Similar to `plot_observer_statistics`, but may use a different style 
    (e.g., Seaborn themes) or layout to visualize the Minimum, Average, 
    and Maximum fitness trends over generations.

    Parameters
    ----------
    observer_file_path : `str` or `Path`
        Path to the CSV statistics file.
    plot_folder_directory : `str` or `Path`
        Directory where the output plot will be saved.

    Returns
    -------
    `None`
    """
    observer_file_path    = Path(observer_file_path)
    plot_folder_directory = Path(plot_folder_directory)

    plot_folder_directory.mkdir(parents = True, exist_ok = True)

    job_id    = observer_file_path.stem.split('_')[-1]
    plot_file = plot_folder_directory / f"generation_plot_{job_id}.png"

    plt.rcParams.update(C.IEEE_PLOT_PARAMS)
    plt.figure(figsize = (C.IEEE_FIGURE_WIDTH, C.IEEE_FIGURE_HEIGHT))
    
    print(f"Generate statistics graphs for job '{job_id}'...")
    with warnings.catch_warnings(), observer_file_path.open('r') as f:
        warnings.filterwarnings("ignore", message = "FigureCanvasAgg is non-interactive")
        analysis.generation_plot(f)
    
    plt.title(f"Fitness Trend - Job {job_id}")
    plt.savefig(plot_file)
    plt.savefig(plot_file.with_suffix(".pdf"))
    plt.close()

    print(f"Graph saved in: {plot_file}")


def plot_energy_vs_hydrophobicity_2(individuals_file: str | Path, plot_folder_directory: str | Path, hydrophobicity_weight: float = C.HYDROPHOBICITY_WEIGHT) -> None:
    """
    Generates a refined scatter plot of Energy vs. Hydrophobicity.

    This version might include additional visual aids, such as:
    - A marginal histogram distribution.
    - A different color palette or marker style.
    - Highlighting of specific Pareto-optimal solutions.

    Parameters
    ----------
    individuals_file : `str` or `Path`
        Path to the file containing the history of individuals.
    plot_folder_directory : `str` or `Path`
        Directory where the output plot will be saved.
    hydrophobicity_weight : `float`
        The weight parameter used in the fitness function.

    Returns
    -------
    `None`
    """
    individuals_file      = Path(individuals_file)
    plot_folder_directory = Path(plot_folder_directory)

    plot_folder_directory.mkdir(parents = True, exist_ok = True)

    if not individuals_file.exists():
        print_error(f"File {individuals_file} not found.", code = 0)
        return
    
    job_id    = individuals_file.stem.split('_')[-1]
    plot_file = plot_folder_directory / f"energy_hydrophobicity_evol_{job_id}.png"
    
    generations_dict = {}
    with individuals_file.open('r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row: continue                # Skip blank lines

            # Parsing Management
            try:
                gen      = int(row[0])
                fitness  = float(row[2])
                sequence = row[3].strip()
            except (ValueError, IndexError):
                continue                        # Skip headers or malformed lines
            
            if gen not in generations_dict:
                generations_dict[gen] = []

            ind                = Individual(candidate = sequence)
            ind.fitness        = fitness
            ind.birthdate      = gen
            ind.hydrophobicity = get_hydrophobicity(sequence)
            ind.energy_val     = fitness - (ind.hydrophobicity * hydrophobicity_weight)

            generations_dict[gen].append(ind)
    
    if not generations_dict:
        print("No data found for the plot.")
        return
    
    plt.rcParams.update(C.IEEE_PLOT_PARAMS)
    fig, ax = plt.subplots(figsize = (12, 8))
    
    max_gen     = max(generations_dict.keys())                      # Get max generation for color normalization
    sorted_gens = sorted(generations_dict.keys())

    # Lists to draw the line of the best by generation
    best_gen_energy = []
    best_gen_hydro  = []

    for gen in sorted_gens:
        individuals = generations_dict[gen]
        if len(individuals) > 1:
            assert isinstance(individuals, list)
            individuals.sort(key = lambda x: x.energy_val)                  # We sort by energy to have a clean line
            
            gen_energies = [ind.energy_val     for ind in individuals]
            gen_hydros   = [ind.hydrophobicity for ind in individuals]
            color_val    = gen / max_gen if max_gen > 0 else 0              # Generation-based color
            
            ax.plot(
                gen_energies, gen_hydros, 
                color     = gradient_map.reversed()(color_val), 
                alpha     = 0.4,                                            # Slight transparency so as not to cover the stitches
                linewidth = 0.8,
                zorder    = 1,                                              # Lines under the dots
            )

    all_individuals = [ind for gen_list in generations_dict.values() for ind in gen_list]
    sc = ax.scatter(
        [ind.energy_val     for ind in all_individuals],
        [ind.hydrophobicity for ind in all_individuals], 
        c          = [ind.birthdate      for ind in all_individuals], 
        cmap       = gradient_map.reversed(), 
        s          = 80, 
        alpha      = 1,
        edgecolors = 'none',
        linewidth  = 0.5,
        zorder     = 2                                                      # Make sure the dots are ABOVE the lines
    )

    for gen in sorted_gens:
        best_in_gen = min(generations_dict[gen], key = lambda x: x.energy_val)
        best_gen_energy.append(best_in_gen.energy_val)
        best_gen_hydro .append(best_in_gen.hydrophobicity)

    ax.plot(best_gen_energy, best_gen_hydro, 
        color     = '#2ECC71',
        linestyle = '--',
        alpha     = 0.8,
        linewidth = 1.5,
        label     = "Best Trajectory",
        zorder    = 3
    )

    global_best = min(all_individuals, key = lambda x: x.fitness)
    ax.annotate(
        f" BEST: {global_best.candidate}\n({global_best.energy_val:.2f})",
        xy         = (global_best.energy_val, global_best.hydrophobicity),
        xytext     = (-20, 15),
        textcoords = 'offset points',
        fontsize   = 9, 
        fontweight = 'bold',
        arrowprops = dict(arrowstyle = '->', color = 'black'),
        ha         = 'right',
        va         = 'bottom', 
        color      = 'black',
        zorder     = 4
    )

    cbar = plt.colorbar(sc, ax = ax)
    cbar.set_label('Generation', rotation = 270, labelpad = 15)

    ax.set_xlabel('Binding Energy (kcal/mol)')
    ax.set_ylabel('Hydrophobicity Index (pH 7)')
    ax.set_title(f'Evolution Landscape (Job {job_id})')
    ax.grid(True, linestyle = '--', alpha = 0.3)

    plt.savefig(plot_file)
    plt.savefig(plot_file.with_suffix(".pdf"))
    plt.close()

    print(f"Energy vs Hydrophobicity plot saved to {plot_file}")


def get_best_sequences_data(csv_path: str, top_n: int = 50) -> pd.DataFrame:
    """
    Parses the individuals CSV to extract the top N unique sequences based on fitness/energy.
    
    Args:
        csv_path (str): Path to `ga_individuals_JOBID.csv`.
        top_n (int): Number of top sequences to retrieve.
        
    Returns:
        pd.DataFrame: DataFrame containing 'sequence' and 'energy' columns.
    """
    try:
        # Load data (assuming standard Inspyred CSV format or custom format)
        # We try to infer column names if headers are missing or specific
        df         = pd.read_csv(csv_path, sep = ", ", header = None)
        df.columns = ["generation", "individual", "fitness", "sequence"]
        
        # Normalize column names for safety
        df.columns = [c.lower().strip() for c in df.columns]
        
        # Identify key columns based on keywords
        seq_col = next((c for c in df.columns if 'candidate' in c or 'seq' in c), None)
        fit_col = next((c for c in df.columns if 'fitness' in c or 'energy' in c), None)
        
        if not seq_col or not fit_col:
            raise ValueError(f"Could not find sequence/fitness columns in {csv_path}")

        # Sort by fitness (ascending because binding energy is negative: lower is better)
        # We drop duplicates to avoid analyzing the same clone multiple times
        best_df = df.sort_values(by = fit_col, ascending = True).drop_duplicates(subset = [seq_col])
        
        return best_df.head(top_n)[[seq_col, fit_col]].rename(columns = {seq_col: 'sequence', fit_col: 'energy'})

    except Exception as e:
        print(f"[PLOT ERROR] Could not extract data from CSV: {e}")
        return pd.DataFrame()


def plot_sequence_consensus(csv_path: str, output_path: str, top_n: int = 50) -> None:
    """
    Generates a Sequence Logo to visualize motif conservation.

    Calculates the information content (in bits) for each position in the 
    peptide sequence based on the top surviving candidates. 
    Uses the `logomaker` library to create the visualization.

    Parameters
    ----------
    individuals_file : `str` or `Path`
        Path to the file containing the individuals.
    plot_folder_directory : `str` or `Path`
        Directory where the output image will be saved.
    top_n : `int`, optional
        Number of top individuals to consider for the consensus (default 75).

    Returns
    -------
    `None`
    """
    df = get_best_sequences_data(csv_path, top_n)
    if df.empty: return

    sequences = df['sequence'].tolist()
    
    # Create count matrix (Position x AminoAcid)
    matrix = lm.transform_matrix(
        lm.alignment_to_matrix(sequences, to_type = 'counts'),
        from_type = 'counts',
        to_type   = 'information'
    )
    
    # Plotting
    plt.figure(figsize = (10, 4))
    
    # 'skylign_protein' creates a standard information-content logo
    logo = lm.Logo(matrix, color_scheme = 'weblogo_protein', stack_order = 'big_on_top', shade_below = .5, fade_below = .5) #skylign_protein
    
    legend_dict = {
        'Hydrophobic': "#000000",
        'Polar'      : "#008000",
        'Basic'      : "#0000FF",
        'Acid'       : "#FF0000",        
    }

    legend_elements = [Patch(facecolor = color, label = label) for label, color in legend_dict.items()]

    logo.ax.legend(
        handles        = legend_elements,
        loc            = 'center left',
        bbox_to_anchor = (1.01, 0.5),
        borderaxespad  = 0,
        frameon        = False
    )

    logo.ax.set_title(f"Sequence Consensus (Top {top_n} Unique Individuals)", fontsize = 12, fontweight = 'bold')
    logo.ax.set_xlabel("Peptide Position", fontsize = 10)
    logo.ax.set_ylabel("Information (Bits)", fontsize = 10)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi = 300)
    plt.close()
    print(f"[PLOTS] Sequence Logo saved to: {output_path}")


def plot_diversity_heatmap(csv_path: str, output_path: str, top_n: int = 50) -> None:
    """
    Visualizes population diversity using a Pairwise Similarity Heatmap.

    1. Selects the top unique candidates.
    2. Computes the similarity score between every pair using the BLOSUM62 matrix.
    3. Plots a heatmap where:
    - High scores (Diagonal/Blocks) indicate conserved motifs.
    - Low scores indicate diversity.

    Parameters
    ----------
    individuals_file : `str` or `Path`
        Path to the file containing the individuals.
    plot_folder_directory : `str` or `Path`
        Directory where the output plot will be saved.
    top_n : `int`, optional
        Number of individuals to compare (default 75).

    Returns
    -------
    `None`
    """
    df = get_best_sequences_data(csv_path, top_n)
    if df.empty: return

    # df          = df.sort_values(by = ["sequence"])
    sequences   = df['sequence'].tolist()
    n           = len(sequences)
    dist_matrix = np.zeros((n, n))

    # Calculate pairwise Hamming distances
    # (Distance = number of mismatches between two strings)
    for i in range(n):
        for j in range(n):
            score = 0
            for a, b in zip(sequences[i], sequences[j]):
                try:
                    score += C.BLOSUM62[a, b]
                except (KeyError, IndexError):
                    score += -4  # Penalità di fallback se il carattere non è in matrice
            dist_matrix[i, j] = score

    # Plotting
    plt.figure(figsize = (8, 7))
    ax  = sns.heatmap(
        dist_matrix, 
        cmap        = "viridis", 
        cbar_kws    = {'label': 'BLOSUM62 Similarity Score'},
        xticklabels = df["sequence"],  # [seq if i == 0 else "" for i, seq in enumerate(df["sequence"])], 
        yticklabels = df["sequence"]   # [seq if i == 0 else "" for i, seq in enumerate(df["sequence"])]
    )

    ax.set_xticklabels(ax.get_xticklabels(), family = 'monospace', fontsize = 6)
    ax.set_yticklabels(ax.get_yticklabels(), family = 'monospace', fontsize = 6)
    
    plt.title(f"Similarity Heatmap (BLOSUM62 - Top {top_n})", fontsize = 12, fontweight = 'bold')
    plt.xlabel("Individual ID (Ranked by Fitness)")
    plt.ylabel("Individual ID (Ranked by Fitness)")
    
    plt.tight_layout()
    plt.savefig(output_path, dpi = 300)
    plt.close()
    print(f"[PLOTS] Diversity Heatmap saved to: {output_path}")


def plot_pca_landscape(csv_path: str, output_path: str, top_n: int = 100) -> None:
    """
    Visualizes the search space exploration using PCA.

    1. Encodes peptide sequences using One-Hot Encoding.
    2. Reduces dimensionality to 2D using Principal Component Analysis (PCA).
    3. Plots the candidates, coloring them by fitness or generation.
    
    Helps to understand if the algorithm is exploring diverse regions 
    or converging to a specific cluster.

    Parameters
    ----------
    individuals_file : `str` or `Path`
        Path to the file containing the individuals.
    plot_folder_directory : `str` or `Path`
        Directory where the output plot will be saved.

    Returns
    -------
    `None`
    """
    # We take more points (top_n=100) to better visualize the landscape
    df = get_best_sequences_data(csv_path, top_n)
    if df.empty: return

    sequences = df['sequence'].tolist()
    energies  = df['energy'].tolist()

    # Convert sequences to character lists for encoding
    # e.g. "MWW" -> [['M', 'W', 'W']]
    seq_matrix = [list(s) for s in sequences]

    # One-Hot Encoding
    # Converts categorical amino acids into a binary vector
    encoder = OneHotEncoder(sparse_output = False, handle_unknown = 'ignore')
    encoded_data = encoder.fit_transform(seq_matrix)

    # PCA Projection (2 Components)
    pca = PCA(n_components = 2)
    pca_result = pca.fit_transform(encoded_data)

    # Plotting
    plt.figure(figsize = (8, 6))
    
    sc = plt.scatter(
        pca_result[:, 0], 
        pca_result[:, 1], 
        c         = energies, 
        cmap      = 'viridis', # Reversed: Dark/Purple = Better Energy (Lower)
        edgecolor = 'k', 
        s         = 60, 
        alpha     = 0.8
    )
    
    cbar = plt.colorbar(sc)
    cbar.set_label("Binding Energy (kcal/mol)", rotation = 270, labelpad = 15)
    
    plt.title(f"PCA Landscape of Top {top_n} Sequences", fontsize = 12, fontweight = 'bold')
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)")
    plt.grid(True, linestyle = '--', alpha = 0.5)

    plt.tight_layout()
    plt.savefig(output_path, dpi = 300)
    plt.close()
    print(f"[PLOTS] PCA Landscape saved to: {output_path}")







def plot_fitness_statistics(observer_file_path: str | Path, plot_folder_directory: str | Path) -> None:
    """
    Plots the progression of fitness metrics (Best, Average, Worst) over time.

    A wrapper or standard implementation to visualize the convergence 
    of the evolutionary algorithm based on the logged statistics.

    Parameters
    ----------
    statistics_file : `str` or `Path`
        Path to the CSV file generated by the file observer.
    output_dir : `str` or `Path`
        Directory where the plot will be saved.

    Returns
    -------
    `None`
    """
    observer_file_path    = Path(observer_file_path)
    plot_folder_directory = Path(plot_folder_directory)
    plot_folder_directory.mkdir(parents=True, exist_ok=True)

    job_id = observer_file_path.stem.split('_')[-1]
    plot_file = plot_folder_directory / f"generation_plot_{job_id}.png"

    print(f"Generate statistics graphs for job '{job_id}'...")

    # 1. Load Data
    try:
        df = pd.read_csv(observer_file_path, header=None)
        # Mapping columns: 0:Gen, 1:PopSize, 2:Worst, 3:Best, 4:Median, 5:Average, 6:StdDev
        df.columns = ['Generation', 'PopSize', 'Worst', 'Best', 'Median', 'Average', 'StdDev']
    except Exception as e:
        print(f"[PLOT ERROR] Could not read statistics file: {e}")
        return

    # 2. Setup Style
    sns.set_theme(style="whitegrid", rc=C.IEEE_PLOT_PARAMS)
    plt.figure(figsize=(6, 4))

    # 3. Plotting Lines (Order matters for layering)
    
    # A) Worst -> RED
    sns.lineplot(
        data=df, x='Generation', y='Worst', 
        color='red', linewidth=1, alpha=0.6, label='Worst'
    )

    # B) Median -> BLUE
    sns.lineplot(
        data=df, x='Generation', y='Median', 
        color='blue', linewidth=1, alpha=0.8, label='Median'
    )

    # C) Average -> PURPLE (Dashed)
    sns.lineplot(
        data=df, x='Generation', y='Average', 
        color='purple', linestyle='--', linewidth=1.5, label='Average'
    )

    # D) Best -> GREEN (Thickest, on top)
    sns.lineplot(
        data=df, x='Generation', y='Best', 
        color='green', linewidth=2, label='Best'
    )

    plt.fill_between(
        df['Generation'], df['Best'], df['Worst'], 
        color='gray', alpha=0.1
    )

    # 4. Formatting
    ax = plt.gca()
    
    # Force X-axis to be Integers only
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    plt.title(f"Fitness Trend - Job {job_id}", fontweight='bold')
    plt.xlabel("Generation")
    plt.ylabel("Fitness (Energy)")
    plt.legend(frameon=True, loc='best')

    # 5. Save
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.savefig(plot_file.with_suffix(".pdf"))
    plt.close()

    print(f"Graph saved in: {plot_file}")

def plot_binding_energy_vs_hydrophobicity(individuals_file: str | Path, plot_folder_directory: str | Path, hydrophobicity_weight: float = C.HYDROPHOBICITY_WEIGHT) -> None:
    """
    Scatter plot correlating Binding Energy (Vina) and Hydrophobicity.

    Displays the population distribution in the objective space.
    - X-axis: Binding Energy (lower is better).
    - Y-axis: Hydrophobicity Index.
    - Color: Usually represents the Generation (Time).

    Parameters
    ----------
    data_file : `str` or `Path`
        Path to the dataset containing individual metrics.
    output_dir : `str` or `Path`
        Directory where the plot will be saved.

    Returns
    -------
    `None`
    """
    individuals_file      = Path(individuals_file)
    plot_folder_directory = Path(plot_folder_directory)
    plot_folder_directory.mkdir(parents=True, exist_ok=True)

    if not individuals_file.exists():
        print(f"Error: File {individuals_file} not found.")
        return
    
    job_id    = individuals_file.stem.split('_')[-1]
    plot_file = plot_folder_directory / f"energy_hydrophobicity_evol_{job_id}.png"

    print(f"Generating Energy/Hydrophobicity plot for Job {job_id}...")

    try:
        df = pd.read_csv(individuals_file, sep=", ", header=None, engine='python')
        df = df.iloc[:, [0, 2, 3]].copy() 
        df.columns = ['generation', 'fitness', 'sequence']
        df['sequence'] = df['sequence'].astype(str).str.strip()
        
        df['hydrophobicity'] = df['sequence'].apply(get_hydrophobicity)
        df['energy_val']     = df['fitness'] - (df['hydrophobicity'] * hydrophobicity_weight)
        
    except Exception as e:
        print(f"Error reading/processing CSV: {e}")
        return

    if df.empty: return

    sns.set_theme(style="whitegrid", rc=C.IEEE_PLOT_PARAMS)
    fig, ax = plt.subplots(figsize=(6, 4))

    # Colormap
    max_gen = df['generation'].max()
    norm = plt.Normalize(df['generation'].min(), max_gen)
    try:
        cmap = gradient_map.reversed()
    except NameError:
        cmap = cm.get_cmap('viridis_r')

    for gen, group in df.groupby('generation'):
        if len(group) > 1:
            sorted_group = group.sort_values('energy_val')
            ax.plot(
                sorted_group['energy_val'], 
                sorted_group['hydrophobicity'],
                color=cmap(norm(gen)), alpha=0.3, linewidth=0.8, zorder=1
            )

    sns.scatterplot(
        data=df, x='energy_val', y='hydrophobicity', hue='generation',
        palette=cmap, s=30, alpha=1.0, edgecolor='none', ax=ax, zorder=2, legend=False
    )

    best_trajectory = df.loc[df.groupby('generation')['energy_val'].idxmin()]
    sns.lineplot(
        data=best_trajectory, x='energy_val', y='hydrophobicity',
        color='#2ECC71', linestyle='--', linewidth=1.5, label='Best Trajectory',
        zorder=3, ax=ax, sort=False
    )

    global_best = df.loc[df['energy_val'].idxmin()]
    best_x = global_best['energy_val']
    best_y = global_best['hydrophobicity']
    x_mid = (df['energy_val'].max() + df['energy_val'].min()) / 2
    
    global_best = df.loc[df['energy_val'].idxmin()]
    best_x = global_best['energy_val']
    best_y = global_best['hydrophobicity']
    x_mid = (df['energy_val'].max() + df['energy_val'].min()) / 2
    
    if best_x < x_mid:
        text_pos_x, text_pos_y = 0.95, 0.90
        ha_align = 'right'
        connection_arc = "arc3,rad=0.15"   
    else:
        text_pos_x, text_pos_y = 0.05, 0.90
        ha_align = 'left'
        connection_arc = "arc3,rad=-0.15"

    ax.annotate(
        f" BEST: {global_best['sequence']}\n({global_best['energy_val']:.2f})",
        xy=(best_x, best_y),
        xytext=(text_pos_x, text_pos_y),
        textcoords='axes fraction',
        fontsize=9, fontweight='bold',
        arrowprops=dict(
            arrowstyle='->', 
            color='black', 
            shrinkB=5, 
            connectionstyle=connection_arc
        ),
        ha=ha_align, va='top', color='black',
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8, ec="none"),
        zorder=10
    )

    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Generation', rotation=270, labelpad=15)

    ax.set_xlabel('Binding Energy (kcal/mol)')
    ax.set_ylabel('Hydrophobicity Index (pH 7)')
    ax.set_title(f'Evolution Landscape (Job {job_id})')
    
    handles, labels = ax.get_legend_handles_labels()
    if handles: 
        ax.legend(
            handles=[handles[0]], 
            labels=[labels[0]], 
            loc='lower left',
            frameon=True,
            framealpha=0.9,
            edgecolor='gray'
        )

    plt.tight_layout()
    plt.savefig(plot_file)
    plt.savefig(plot_file.with_suffix(".pdf"))
    plt.close()

    print(f"Energy vs Hydrophobicity plot saved to {plot_file}")