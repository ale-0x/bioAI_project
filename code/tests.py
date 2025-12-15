# tests.py

# Puoi aggiungere questo blocco temporaneamente a ga_problem.py o main.py per testare
import random
from ga_problem import peptide_generator
from constants import PEPTIDE_LENGTH, POPULATION_SIZE

rand = random.Random()
rand.seed(99) # Per la riproducibilità del test

# Genera una popolazione di test:
pop_test = [peptide_generator(rand, {'peptide_length': PEPTIDE_LENGTH}) 
            for _ in range(POPULATION_SIZE)]

print("Popolazione iniziale di test:")
for i, seq in enumerate(pop_test):
    print(f"{i+1}: {seq} (Length: {len(seq)})")
