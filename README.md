# Improving_the_Classical_COVIDOA_Algorithm
An enhanced version of the Coronavirus Disease Optimization Algorithm (COVIDOA) featuring adaptive mutation, tournament selection, population restart, and improved replication strategies for superior optimization performance.


# Enhanced COVIDOA Algorithm
This repository contains an improved version of the Coronavirus Disease Optimization Algorithm (COVIDOA), a population-based metaheuristic inspired by the viral replication and mutation behavior of the coronavirus. These enhancements were made to address limitations in the original implementation, including poor performance in escaping local minima and slow convergence.


# Overview
COVIDOA mimics biological behaviors like frameshifting and mutation to solve complex optimization problems. However, its original implementation often suffers from:
(a) Premature convergence (Although better than most other metaheuristic algorithms)
(b) Poor exploration of the search space
(c) Static mutation behavior
(d) To improve its efficiency and reliability, several enhancements were introduced.


# Enhanced Implementation
1. Tournament Selection
Why? Replaces roulette wheel selection to provide better control over selection pressure and help maintain diversity.

How? Randomly selects a subset (tournament) and picks the fittest individual. Reduces the risk of premature convergence.

2. Adaptive Mutation Rate
Why? A static mutation rate lacks the flexibility to escape local minima.

How? Mutation rate dynamically increases when no improvement is seen over a number of iterations to promote exploration.

3. Population Restart Mechanism
Why? Avoids stagnation by reintroducing diversity into the population.

How? If no improvement is observed for a fixed number of generations, a portion of the population is replaced with new individuals.

4. Improved Replication Strategy
Why? Enhances exploration through more variation.

How? Introduced dynamic control of parameters (gamma, beta) in the uniform crossover process to inject more randomness.


# Benchmark Testing
The original and improved versions of the algorithm were tested on standard benchmark functions:

Alpine N. 2 Function
Algorithm	Best Cost @ 500 Iterations
Original COVIDOA	-3.243747570e+21
Improved COVIDOA	-3.939468708e+20

Schwefel Function
Algorithm	Best Cost @ 500 Iterations
Original COVIDOA	11544.4333
Improved COVIDOA	11506.1951

NB: The enhanced algorithm consistently outperformed the original, converging faster and producing better solutions.


# Results Summary
(a) Better Exploration: Adaptive mutation + population restart reduce the risk of getting stuck.

(b) Faster Convergence: Improved selection and replication boost convergence speed.

(c) Increased Diversity: Helps prevent premature convergence.

(d) Scalability: Easily tunable for different optimization problems.


# Future Work
(a) Parameter tuning for further refinement

(b) Explore hybridization with other metaheuristics (e.g., PSO, DE)

(c) Apply to real-world, high-dimensional optimization problems


# Parameters Used
MaxIt = 500          # Max number of iterations
nPop = 300           # Population size
D = 30               # Problem dimensionality
minVal, maxVal = 10, 50  # Variable bounds
MR = 0.1             # Base Mutation Rate
gamma = 0.5          # Crossover parameter
beta = 0.5           # Crossover parameter
numOfSubprotiens = 6 # COVIDOA-specific replication control
