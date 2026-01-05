# bioAI_project

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
The final candidate is returned along with a structural prediction and is available for further inspection.
