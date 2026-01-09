import matplotlib.pyplot as plt
import inspyred.ec.analysis
from inspyred.ec import Individual
import constants as C
import random
import os
from matplotlib.colors import LinearSegmentedColormap

gradient_map = LinearSegmentedColormap.from_list("custom_gradient", ['#FF0000', '#FFFF00'])

def plot_observer_statistics(observer_file_path: str) -> None:

    job_id = os.path.basename(individuals_file).split('_')[-1].split('.')[0]
    plot_file = f"plots/generation_plot_{job_id}.png"
    inspyred.ec.analysis.generation_plot(open(observer_file_path, 'r'))
    plt.savefig(plot_file)
    plt.close()

def plot_energy_vs_hydrophobicity(individuals_file: str) -> None:
    """
    scatterplot: energia vs idrofobicità media di ciascun individuo attraverso le generazioni.
    puoi chiamarlo con un individuals_file.csv
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    job_id = os.path.basename(individuals_file).split('_')[-1].split('.')[0]

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
            ind.hydrophobicity = C.get_hydrophobicity(individual)
            generations_dict[gen].append(ind)

    # Get max generation for color normalization
    max_gen = max(gen for gen in generations_dict.keys()) if generations_dict else 1

    # Plot individuals with color based on age
    for ind in [ind for array in generations_dict.values() for ind in array]:
        energy = ind.fitness - ind.hydrophobicity*C.HYDROPHOBICITY_WEIGHT  # x-axis: energy
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
            energies = [ind.fitness for ind in individuals]
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

    scatter_plot_file = f"plots/energy_hydrophobicity_{job_id}.png"
    plt.savefig(scatter_plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Energy vs Hydrophobicity plot saved to {scatter_plot_file}")