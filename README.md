# RCGA

A Real-Coded Genetic Algorithm (RCGA) implementation in Ruby for optimizing linear classifier parameters (weights and bias).

## Overview

This program uses a real-coded genetic algorithm to search for optimal parameters of a linear discriminant function. Given a dataset of labeled feature vectors, RCGA evolves a set of weight parameters that maximize classification performance (precision, recall, or F-measure).

The discriminant function takes the form:

```
y = w1*x1 + w2*x2 + ... + wN*xN + bias
```

If `y > 0`, the sample is classified as label A; otherwise, label B.

## Requirements

- Ruby (or Docker)

## Usage

```
ruby RCGA.rb -c <config.yml> -i <input.tsv> -o <output.txt>
```

| Flag | Required | Description |
|------|----------|-------------|
| `-c` | Yes | YAML configuration file |
| `-i` | Yes | Input file (TSV with pre-computed feature values) |
| `-o` | Yes | Output file for optimized parameters |
| `-h` | No | Show help message |

### Input File Format

Tab-separated values where:
- Column 0: ground-truth label (must match labels defined in config)
- Columns 1..N: feature values

### Output File Format

The output file contains the best fitness score and the optimized parameter values.

## Configuration (YAML)

See `RCGA.yml` for a full example. Key settings include:

| Section | Parameter | Description |
|---------|-----------|-------------|
| `input.parameter_num` | int | Number of parameters (features) |
| `input.labels` | list | Classification labels (e.g., `["A", "B"]`) |
| `parameter_ranges` | list of [min, max] | Search range for each weight parameter |
| `bias_range` | [min, max] | Search range for the bias term |
| `ga.max_generation` | int | Maximum number of generations |
| `ga.initialization.population_num` | int | Population size |
| `ga.evaluation.opt_metric` | string | Metric to optimize: `precision`, `recall`, or `f_mesure` |
| `ga.evaluation.tolerance` | float | Convergence threshold |
| `ga.selection.strategy` | string | Selection strategy: `elite`, `roulette`, or `hybrid` |
| `ga.crossover.type` | string | Crossover method: `BLX-alpha` or `SPX` (Simplex) |
| `ga.crossover.prob` | float | Crossover probability |
| `ga.mutation.type` | string | Mutation type: `uniform` or `boundary` |
| `ga.mutation.prob` | float | Mutation probability |

## Algorithm

1. **Initialization** — Generate a random initial population within the specified parameter ranges.
2. **Evaluation** — Compute the fitness of each individual using the selected metric (precision, recall, or F-measure).
3. **Selection** — Select parents using one of three strategies:
   - **Elite**: Select the top-ranked individuals.
   - **Roulette**: Fitness-proportionate random selection.
   - **Hybrid**: Elite selection for a portion, roulette for the rest.
4. **Crossover** — Generate offspring using one of two methods:
   - **BLX-alpha**: Blend crossover with parameter `alpha`.
   - **SPX**: Simplex crossover using `n+1` parents to produce a child in n-dimensional space.
5. **Mutation** — Apply mutation with a given probability:
   - **Uniform**: Replace a gene with a random value within its range.
   - **Boundary**: Set a gene to either its upper or lower bound.
6. **Convergence Check** — If the relative change in average fitness of the top individuals falls below the tolerance, stop. Otherwise, return to step 2.

## Running with Docker

Build the image:

```bash
docker build -t rcga .
```

Run with the included test data:

```bash
docker run --rm -v $(pwd)/output:/output rcga \
  -c test_config.yml -i test_input.tsv -o /output/result.txt
```

Example output:

```
I, [...]  INFO -- : *** start RCGA ***
I, [...]  INFO -- : first generation
I, [...]  INFO -- : convergence!!
I, [...]  INFO -- : max f_mesure: 1.0
I, [...]  INFO -- : *** stop RCGA ***
```

## Test Data

A sample configuration (`test_config.yml`) and input file (`test_input.tsv`) are included for quick testing. The test data contains 20 linearly separable samples (10 per class) with 3 features.
