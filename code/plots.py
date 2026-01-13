# plots.py

import matplotlib
matplotlib.use("Agg")                   # Usa backend non interattivo per salvataggio file senza display

import csv
import inspyred.ec.analysis
import matplotlib.pyplot as plt
import os
import random

from inspyred.ec        import Individual
from matplotlib.colors  import LinearSegmentedColormap


import constants as C

from peptide_operators  import get_hydrophobicity
from utils              import print_error

gradient_map = LinearSegmentedColormap.from_list("custom_gradient", ['#FF0000', '#FFFF00'])

def plot_observer_statistics(observer_file_path: str) -> None:

    job_id = os.path.basename(observer_file_path).split('_')[-1].split('.')[0]
    plot_file = f"../plots/generation_plot_{job_id}.png"
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
            ind.hydrophobicity = get_hydrophobicity(individual)*C.HYDROPHOBICITY_WEIGHT
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

    scatter_plot_file = f"../plots/energy_hydrophobicity_{job_id}.png"
    plt.savefig(scatter_plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Energy vs Hydrophobicity plot saved to {scatter_plot_file}")

def plot_observer_statistics_2(observer_file_path: str) -> None:
    if not os.path.exists(observer_file_path):
        print_error(f"File {observer_file_path} non trovato.", code = 0)
    else:
        job_id = os.path.basename(observer_file_path).split('_')[-1].split('.')[0]
        plot_file = f"../plots/generation_plot_{job_id}.png"

        print(f"Generazione grafico statistiche per Job '{job_id}'...")

        with open(observer_file_path, 'r') as f:
            data = f.readlines()
        
        plt.figure(figsize = (10, 6))
        inspyred.ec.analysis.generation_plot(data)
        plt.title(f"Andamento Fitness - Job {job_id}")
        plt.savefig(plot_file, dpi = 300, bbox_inches = 'tight')
        plt.close()
        
        print(f"Grafico salvato in: {plot_file}")

def plot_energy_vs_hydrophobicity_2(individuals_file: str) -> None:
    """
    scatterplot: energia vs idrofobicità media di ciascun individuo attraverso le generazioni.
    puoi chiamarlo con un individuals_file.csv
    """
    if not os.path.exists(individuals_file):
        print_error(f"File {individuals_file} non trovato.", code = 0)
        return
    
    job_id = os.path.basename(individuals_file).split('_')[-1].split('.')[0]
    
    generations_dict = {}
    with open(individuals_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row: continue                # Salta righe vuote

            # Gestione del parsing
            try:
                gen      = int(row[0])
                fitness  = float(row[2])
                sequence = row[3].strip()
            except ValueError:
                continue                        # Salta header o righe malformate
            
            if gen not in generations_dict:
                generations_dict[gen] = []

            ind                = Individual(candidate = sequence)
            ind.fitness        = fitness
            ind.birthdate      = gen
            ind.hydrophobicity = get_hydrophobicity(sequence)
            ind.energy_val     = fitness - (ind.hydrophobicity * C.HYDROPHOBICITY_WEIGHT)

            generations_dict[gen].append(ind)
    
    if not generations_dict:
        print("Nessun dato trovato per il plot.")
        return
    
    fig, ax = plt.subplots(figsize = (12, 8))
    max_gen = max(generations_dict.keys())                      # Get max generation for color normalization

    # Liste per tracciare la linea del migliore per generazione
    best_gen_energy = []
    best_gen_hydro  = []

    all_individuals = [ind for gen_list in generations_dict.values() for ind in gen_list]
    sorted_gens     = sorted(generations_dict.keys())

    for gen in sorted_gens:
        individuals = generations_dict[gen]
        if len(individuals) > 1:
            individuals.sort(key=lambda x: x.energy_val)                # Ordiniamo per energia per avere una linea pulita
            
            gen_energies = [ind.energy_val     for ind in individuals]
            gen_hydros   = [ind.hydrophobicity for ind in individuals]
            
            color_val = gen / max_gen if max_gen > 0 else 0             # Colore basato sulla generazione
            
            ax.plot(
                gen_energies, gen_hydros, 
                color     = gradient_map.reversed()(color_val), 
                alpha     = 0.4,                                        # Leggera trasparenza per non coprire i punti
                linewidth = 1
            )

    energies = [ind.energy_val     for ind in all_individuals]
    hydros   = [ind.hydrophobicity for ind in all_individuals]
    colors   = [ind.birthdate      for ind in all_individuals]

    sc = ax.scatter(
        energies, hydros, 
        c          = colors, 
        cmap       = gradient_map.reversed(), 
        s          = 80, 
        alpha      = 0.8,
        edgecolors = 'grey',
        linewidth  = 0.5,
        zorder     = 10                                                 # Assicura che i pallini siano SOPRA le linee
    )

    for gen in sorted_gens:
        best_in_gen = min(generations_dict[gen], key = lambda x: x.energy_val)
        best_gen_energy.append(best_in_gen.energy_val)
        best_gen_hydro.append(best_in_gen.hydrophobicity)

    ax.plot(best_gen_energy, best_gen_hydro, 
        color     = 'green',
        linestyle = '--',
        alpha     = 0.6,
        linewidth = 1.5,
        label     = "Best Trajectory",
        zorder    = 11
    )

    global_best = min(all_individuals, key = lambda x: x.fitness)
    ax.text(
        global_best.energy_val, 
        global_best.hydrophobicity, 
        f" BEST: {global_best.candidate}\n({global_best.energy_val:.2f})", 
        fontsize   = 9, 
        fontweight = 'bold',
        ha         = 'right', 
        va         = 'bottom',
        color      = 'black',
        zorder     = 12
    )

    cbar = plt.colorbar(sc, ax = ax)
    cbar.set_label('Generation', rotation = 270, labelpad = 20)

    ax.set_xlabel('Binding Energy (kcal/mol) [Lower is Better]')
    ax.set_ylabel('Hydrophobicity Index (pH 7)')
    ax.set_title(f'Evolution Landscape: Energy vs Hydrophobicity (Job {job_id})')
    ax.grid(True, linestyle = '--', alpha = 0.3)

    scatter_plot_file = f"../plots/energy_hydrophobicity_{job_id}.png"
    plt.savefig(scatter_plot_file, dpi = 300, bbox_inches = 'tight')
    plt.close()

    print(f"Energy vs Hydrophobicity plot saved to {scatter_plot_file}")