# Pebble Accretion and the Density Evolution of Kuiper Belt Objects
This repository is used to study the effects of pebble accretion on small, low density, Kuiper belt objects (KBOs).

## The Kuiper belt

The Kuiper belt is a band of icy planetesimals and dwarf planets that lies just beyond the orbit of Neptune. It is a favorable place to study planet and planetesimal formation because many of these objects -- particularly the Cold Classical KBOs are thought to have been undisturbed since its formation, serving as celestial ***fossils*** from which we can extract information about the young solar nebula. Owing to the distance between us and the Kuiper belt, it can be difficult to measure the **mass** and **size** of these objects to any reasonable accuracy. However, if a Kuiper belt object happens to be part of a binary (i.e., it *has* or *is* a moon), we can measure get the mass of the objects by taking a look at their interactions. The size of a KBO can be calculated in the event of stellar occultation, i.e., the KBO passes in front of a background star. 

Kuiper belt objects in which we have been lucky enough to obtain both the size and mass, allows us to calculate the density of Kuiper belt objects. With size, distance and density, astronomers have revelead in interesting trend; small KBOs have low densities (0.5 g cm<sup>-3</sup>) and large KBOs have high densities (2.5 g cm<sup>-3</sup>). If planets/planetesimals are consequence of the accumulation of pebbles in proximity to one another, their is no reason to assume that there would be such a drastic change in density from one size to another. The only exception is compaction through gravity; the gravitational pressure from larger masses is better able to *squeeze* away any 'empty space' in the planetesimal. However, this process would, at best, increase the density by a factor of 2, definitely not matching the whopping 5x difference seen in the Kuiper belt.


We explain the trend in the following way: planetesimals beyond the ice belt are formed primarily icy, and then become rocky through the process of ***pebble accretion***. **Calculating the effects of pebble accretion on icy planetesimals is the purpose of this repository.**

## How it works

Between 2020 and 2021, I ran a series of simulations study planetesimal formation via the **streaming instability**. Our simulations included two separate dust species, one representing ices and the other representing silicates. The result, was that planetesimals formed via the streaming instability were mostly icy. In this repository, you will find a csv file containing the ice and rock fraction of each planetesimal formed. These are the planetesimals I use to study the impact pebble accretion has on the density of KBOs. 


**For more information on the streaming instability, check out these articles:**    
[Youdin & Goodman (2005)](https://ui.adsabs.harvard.edu/link_gateway/2005ApJ...620..459Y/EPRINT_PDF)  
[Johansen et al. (2007)](https://ui.adsabs.harvard.edu/link_gateway/2007Natur.448.1022J/EPRINT_PDF)   
[Carrera et al. (2015)](https://ui.adsabs.harvard.edu/link_gateway/2015A%26A...579A..43C/EPRINT_PDF) 
[Johansen et al. (2015)](https://ui.adsabs.harvard.edu/link_gateway/2015SciA....1E0109J/EPRINT_PDF)    
[Nesvorný &  Vokrouhlický (2019)](https://ui.adsabs.harvard.edu/link_gateway/2019Icar..331...49N/EPRINT_PDF)    
[Li & Youdin (2021)](https://ui.adsabs.harvard.edu/link_gateway/2021ApJ...919..107L/EPRINT_PDF)

**For more information on pebble accretion, checkout this review by Anders Johansen and Michiel Lambrechts:**     
[Johansen & Lambrechts (2017)](https://www.annualreviews.org/doi/10.1146/annurev-earth-063016-020226)

