# plots.py

import matplotlib
matplotlib.use("Agg")                   # Use non-interactive backend for file saving without display

import csv
import matplotlib.pyplot as plt
import warnings

from inspyred.ec        import analysis, Individual
from matplotlib.colors  import LinearSegmentedColormap
from pathlib            import Path


import constants as C

from peptide_operators  import get_hydrophobicity
from utils              import print_error

gradient_map = LinearSegmentedColormap.from_list("custom_gradient", ['#FF0000', '#FFFF00'])

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
    """
    scatterplot: energia vs idrofobicità media di ciascun individuo attraverso le generazioni.
    puoi chiamarlo con un individuals_file.csv
    """
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
    scatterplot: energia vs idrofobicità media di ciascun individuo attraverso le generazioni.
    puoi chiamarlo con un individuals_file.csv
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
        alpha      = 0.8,
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