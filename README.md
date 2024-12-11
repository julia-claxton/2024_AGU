# Steady-State Pitch Angle Distributions of Radiation Belt Electrons Subject to Multiple Bounces
### Julia Claxton (she/they), University of Colorado Boulder, *julia.claxton@colorado.edu*
---

This repository contains the code necessary to replicate the results of my work presented at the AGU 2024 meeting, as well as a copy of the poster and all figures used therein.

To view the poster and figures from this project, please see `./results` and the files & directories within. To see the code used to generate this data, see `./code` To perform a recalculation of the data used in these results, please email me for detailed instructions.

Feel free to email me with any questions, comments, or collaboration requests!

**Abstract: SM11D-2752**
> Pitch angle distributions (PADs) are an important metric in understanding magnetospheric wave-particle interactions. How a PAD varies with energy can provide key insight into energy-dependent scattering processes responsible for causing energetic particle precipitation.
> 
> When particles scattered into the loss cone collide with the atmosphere, a non-trivial fraction is scattered back along the field line in a new pitch angle and energy configuration. These backscattered particles form a new pitch angle distribution largely seen in the anti-loss cone. Traditionally, backscattered populations have been neglected in studies of wave-particle interactions. However, backscattered particles may play a role in shaping the PADs measured in orbit, and therefore untangling their dynamics is important for future studies of radiation belt wave-particle interactions. Since backscattered populations tend to reside in the anti-loss cone, they will precipitate at the magnetically conjugate hemisphere. This creates a second backscattered population, which will mirror back into the loss cone of the initial hemisphere. Therefore, loss and anti-loss cone measurements from orbit are not “clean” pictures of wave-particle interactions, but rather a composite sum of an initial input PAD and a series of atmospheric backscatters. The degree to which these backscatters impact measured PADs has not yet been explored.
> 
> In this work, we conduct multibounce particle simulations for radiation belt electrons using Geant4, creating composite steady-state distributions for a variety of typical input distributions. The impact of backscattered populations on the overall PAD is evaluated, as well as their atmosphere-driven lifetime before fully precipitating or reaching a steady state. The simulated composite distributions are compared to in-situ PADs measured by ELFIN, allowing for an evaluation and validation of the Geant4-based precipitation modeling commonly used to simulate the interaction of precipitating electrons with the atmosphere. Finally, the viability of inverse methods to recover an initial PAD from in-situ PADs is explored.