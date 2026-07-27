# MutationShelteringTheoryV2

Scripts used to perform the analyses in the manuscript:

> *Stepwise expansion of recombination suppression on sex chromosomes and other supergenes through lower load advantage and deleterious mutation sheltering*

For any question regarding the analyses or the scripts, please contact: **paul.yann.jay [at] gmail.com**

---

## Overview

Simulations were performed with [SLiM](https://messerlab.org/slim/), and fall into two categories:

1. **Analyses of the fate of many inversions** (Figures 3, S2–3, S10–19, S22–23, S25)
2. **Analyses of the formation of sex chromosomes** (Figures 5, 6, S24)

Deterministic simulations and figures were produced with R.

The original datasets from Jay et al. are available on figshare ([doi:10.6084/m9.figshare.19961033](https://doi.org/10.6084/m9.figshare.19961033)). The new datasets will be published upon publication of the manuscript; in the meantime, they are available upon request.

---

## Requirements

- [SLiM](https://messerlab.org/slim/) (tested with version 3/4 — see individual scripts for details)
- Perl (for the output-parsing scripts)
- R (for deterministic simulations and figures)
- [GNU parallel](https://www.gnu.org/software/parallel/) (optional, for running many simulations concurrently)
- [Snakemake](https://snakemake.readthedocs.io/) (for the pipeline introducing inversions into initial populations)

---

## 1. Analyses of the fate of many inversions

*(Figures 2c,3,4, S2–8,10–12,14)*

The workflow has two steps: first, "burn-in" simulations create initial populations under various parameter combinations; then, many inversions (or other recombination modifiers) are introduced into these initial populations.

### Step 1 — Burn-in: creating initial populations

Initial populations are created with:

- `ScriptNeutralInversion_DefineInitialAutosome_AllDistrib.slim` — basic simulation

**Basic usage:**

```bash
# Single simulation
slim -d N=1000 -d s=-0.01 -d h=0.1 -d r=1e-08 -d mu=1.45e-8 \
     ScriptNeutralInversion_DefineInitialStateAutosome_AllDistrib.slim

# Multiple simulations in parallel (GNU parallel)
parallel -j 2 slim -d s={3} -d mu={1} -d N=10000 '-d h="{2}"' -d Rep=1 -d r=1.0e-07 \
     ScriptNeutralInversion_DefineInitialStateAutosome_AllDistrib.slim \
     ::: 1.0e-08 ::: '"Distrib2"' ::: -0.0001 -0.001
```

**Parameters:**

| Parameter | Description |
|-----------|-------------|
| `mu`      | Mutation rate |
| `r`       | Recombination rate |
| `s`       | Selection coefficient (fixed value, or mean of a gamma distribution) |
| `h`       | Dominance coefficient (see below) |
| `N`       | Population size |

`h` accepts three kinds of values:
- A **float** — the dominance coefficient is fixed for all mutations.
- `"Distrib1"` — dominance coefficients are drawn with equal probability from {0, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5}.
- `"Distrib2"` — dominance coefficients are drawn with equal probability from {0.1, 0.2, 0.3, 0.4, 0.5}.

This writes the initial populations to a `../InitialState` directory, **which must be created beforehand**.

Other burn-in scripts are used for other types of simulations:

| Script | Description |
|--------|-------------|
| `ScriptNeutralInversion_DefineInitialState_HaploDiplo_SingleChrom.slim` | Haplodiplontic life cycle; from Jay et al. 2022 |
| `ScriptNeutralInversion_DefineInitialState_TwoChromFus.slim` | Chromosome fusion (Figure S18); from Jay et al. 2022 |
| `ScriptNeutralInversion_DefineInitialState_XY_Variable_sh.slim` | Variable *h* and *s*, drawn from a distribution; from Jay et al. 2022 |
| `ScriptNeutralInversion_DefineInitialState.slim` | Non-XY system with two mating types, all individuals always heterozygous (Figure S14) |

### Step 2 — Introducing inversions

Thousands of inversions (or other recombination modifiers) are then introduced into the initial populations, using the Snakemake pipeline:

```bash
SnakeMakeSlimSimulationV5 --jobs 20 --configfile ParamFile_Distrib.yaml \
    --config outputDir="Output/" initDir="InitialState/"
```

`ParamFile_Distrib.yaml` specifies the parameter ranges to be simulated and the number of inversions to simulate (see `ParamFile_Distrib.yaml` for an example). This pipeline runs `Script_IntroInv_XY_AllChrom_AllH_SingleRepV2.slim`, which covers most simulation scenarios.

Older scripts from Jay et al. 2022 handle other simulation types:

| Script | Description |
|--------|-------------|
| `Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_SingleChrom_HaploDiplo.slim` | Haplodiplontic life cycle (Figures S22–23) |
| `Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_Variablesh.slim` | Variable *h* and *s*, drawn from a distribution (Figure S17) |
| `Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY_ChromFus.slim` | Chromosome fusion (Figure S18) |
| `Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_HomoHeteroRecombMod.slim` | Non-XY system with two mating types; introduces recombination modifiers (rather than inversions) that suppress recombination regardless of genotype (Figure S14) |

**Basic usage of the "IntroduceInversionFromInitStat" scripts (Jay et al. 2022):**

```bash
slim -d N=1000 -d mu=1e-08 -d h=0.01 -d s=-0.01 -d r=1e-5 -d rep=1 \
     -d start=1 -d end=2000000 \
     -d "Init='../InitialState/slim_g15000_N=1000_r=1e-05_u=1e-08_s=-0.01_h=0.01_1635828730194.txt'" \
     -d "SexChrom='Y'" \
     Script_IntroduceInversionFromInitStat_IndivSimulPlot_BigChrom_NMut_XY.slim
# Introduces a single 2,000,000 bp inversion on the Y chromosome in an initial
# population. `rep` is only used to label the simulation (simulation number "Rep").
```

The inversion's position is defined by `start` and `end`. To mimic a chromosome fusion, set `start` and `end` to span the gap between the two chromosomes — e.g., for a gap at position 10 Mb, use `start=7000000, end=11000000`.

These scripts can be run with a modified Snakemake file (not provided) or with a job scheduler such as [GNU parallel](https://www.gnu.org/software/parallel/).

### Outputs

Simulations write their output to `../Output` (**must be created beforehand**), for example:

```
IntroInv_XYsyst_5M_SexChrom_N10000_r1.0e-07_u1.0e-09_s-0.0001_h0_sInv0_invSize1000000_All50000Rep_Stat.txt
```

This file contains the trajectories of 50,000 1-Mb sex-linked inversions simulated with `N=10000`, `s=-0.0001`, `h=0`, etc., produced by `Script_IntroInv_XY_AllChrom_AllH_SingleRepV2.slim`. Other scripts produce output files with slightly different naming conventions.

Before plotting, concatenate the output files:

```bash
cat IntroInv_XYsyst_5M_*N1* > IntroInv_XYsyst_5M_AllParam_All50000Rep_Stat_Full.txt
```

---

## 2. Analyses of the formation of sex chromosomes

*(Figures 5, 6, S24)*

Run `ScriptFormationXYChromosome_VarGamma_OnlyXY_Revers_Optimized.slim` directly with the desired parameters, for example:

```bash
slim -d N=10000 -d mu=1e-09 -d s=-0.03 -d r=1e-06 -d rep=5 -d N_BP=10000 \
     -d InvR=0.00000001 -d MaxSizeInv=20000000 \
     ScriptFormationXYChromosome_VarGamma_OnlyXY_Revers_Optimized.slim
```

| Parameter    | Description |
|--------------|-------------|
| `N_BP`       | Number of breakpoints in the genome |
| `InvR`       | Inversion rate |
| `MaxSizeInv` | Maximum inversion size |

### Parsing outputs

Simulation outputs need to be parsed before plotting, using the `Parse*` scripts:

- `ParseInvFreqOutput_Rev.pl`
- `ParseRecombinationOutput_SepSex.pl`

```bash
./ParseInvFreqOutput.pl \
    -i N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_InvFreq_IndivSimulation_XY.txt \
    -o N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_InvFreq_IndivSimulation_XY.Parsed.txt

./ParseRecombinationOutput.pl \
    -i N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.txt \
    -o N=1000_u=1e-08_r=1e-05_MaxSizeInv=50000000_Rep_6_Nrecomb_IndivSimulation_XY.Parsed.txt
```

---

## Deterministic simulations and figures

Deterministic simulations and the corresponding figures can be generated with the R scripts provided in this repository (see comments within each script for details).

---

## Data availability

The original datasets from Jay et al. are available on figshare: [doi:10.6084/m9.figshare.19961033](https://doi.org/10.6084/m9.figshare.19961033). The new datasets described in this manuscript will be made publicly available upon publication; until then, they are available upon reasonable request.

---

## Citation

If you use these scripts, please cite:

> *Stepwise expansion of recombination suppression on sex chromosomes and other supergenes through lower load advantage and deleterious mutation sheltering* (full citation to be added upon publication)

## Contact

For any question regarding the analyses or the scripts, please contact **paul.yann.jay [at] gmail.com**.
