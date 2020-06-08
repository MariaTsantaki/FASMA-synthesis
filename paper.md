---
title: 'FASMA 2.0: A Python package for stellar parameters and chemical abundances'
tags:
  - Python
  - astronomy
  - spectroscopy
  - stellar parameters
  - chemical abundances
authors:
  - name: Maria Tsantaki
    orcid: 0000-0002-0552-2313
    affiliation: "1, 2"
  - name: Daniel Andreasen
    affiliation: "2, 3"
  - name: Guilherme Teixeira
    affiliation: 2
affiliations:
 - name: INAF -- Osservatorio Astrofisico di Arcetri, Largo E. Fermi 5, 50125 Firenze, Italy
   index: 1
 - name: Instituto de Astrofísica e Ciências do Espa\c{c}o, Universidade do Porto, CAUP, Rua das Estrelas, Porto, 4150-762, Portugal
   index: 2
 - name: Department of Molecular Medicine, Aarhus University Hospital, Aarhus, Denmark
   index: 3
date: 13 January 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Effective temperature (Teff), surface gravity (logg), and metallicity ([M/H]) are basic stellar atmospheric parameters necessary to characterize a star. Once these parameters are obtained, we can in turn, infer their chemical abundances of various elements and in conjunction with evolutionary models to estimate their evolution, i.e., mass and radius. In this work, we use spectroscopy as a powerful tool to extract this information from stellar atmospheres applied to stars with spectral type FGK both dwarfs and giants. The growing number of spectroscopic surveys dedicated to the study of the Galactic stellar populations has inflated the number of high quality spectra to several hundreds of thousands. This amount is expected to multiply with the forthcoming surveys, such as WEAVE [@deJong:2019] and 4MOST [@Dalton:2014]. The success of these surveys highly depends on the analysis tools to exploit sufficiently all spectral information. Moreover, it is a well-known axiom in exoplanetary studies that one can only determine the planetary properties once the ones of the host star are known. The planetary properties such as mass, radius, composition, are directly linked to their hosts and therefore, robust tools for the derivation of these parameters are necessary.

There are spectral packages available in the literature based on different methods to derive the atmospheric parameters. A standard method for solar-like stars relies on the measurements of equivalent widths of iron lines and by imposing excitation and ionization equilibrium [e.g., @Magrini:2013; @Andreasen:2017; @Tabernero:2019]. Other methods are based on fitting synthetic spectra created on-the-fly with observations under a minimization procedure [e.g., (SME) @Valenti:1996; (iSpec) @Blanco:2014]. Each of the above methods has different limitations depending on the accuracy of the provided atomic data of the line lists, the resolution of the spectrographs, on the quality of the spectra (e.g., signal-to-noise) but also due to the star itself (e.g., rotation, spectral type).

Our goal is to provide a tool to determine accurately and precisely the stellar parameters and chemical abundances to fulfill the above purposes. ``FASMA`` is a Python package to derive the main stellar atmospheric parameters based on the spectral synthesis technique. The principle of this technique relies on the comparison of synthetic spectra with observations to yield the best-fit parameters under a $\chi^{2}$ minimization process:

 $$\chi^{2} = \Sigma \frac{(obs_{i} - synth_{i})^{2}}{\sigma_{i}}$$

where obs is the observed spectrum, synth the synthetic, and $\sigma$ the error on the observed flux for each spectral point, i. The best parameters are the ones that minimize the non-linear least-squares based on the Levenberg-Marquardt algorithm [@Marquardt:1963; @Markwardt:2009].

The synthetic spectra, i.e., the fluxes at the top of the photosphere, are created on-the-fly with the radiative transfer code, [``MOOG``](https://www.as.utexas.edu/~chris/moog.html) [@Sneden:1973] in Fortran for a set of stellar parameters. The synthetic spectra are later convolved with different rotational profiles, such as macroturbulence (vmac), projected rotational velocity (vsini), and instrumental broadening (resolution) to match the observations. To synthesize a spectrum, the information on how the physical properties (such as temperature, electron pressure, gas pressure, etc.) behave at each layer (optical depth) of the atmosphere is necessary. This information is provided in a tabulated form by interpolating from a grid of pre-computed stellar atmospheres for a set of stellar parameters (Teff, logg, [M/H]). ``FASMA`` includes two sets of grids for this analysis: Kurucz [@Kurucz:1993] or marcs models [@Gustafsson:2008].

A key component for spectral synthesis is the list of atoms and molecules which are included in the wavelength intervals. The spectral line list required for the derivation of stellar parameters is described by @Tsantaki:2018 and for the chemical abundances is taken from @Adibekyan:2015 and @Delgado:2015.

Before we compare with observations, we perform a local normalization for the regions of the synthesis and filter out cosmic rays. The procedure to derive stellar parameters is presented in detail by [@Tsantaki:2018] and has been tested for medium and high resolution spectrographs. In the new version 2.0, we now include a new feature to derive chemical abundances for the following elements: Li, Na, Mg, Al, Si, Ca, Sc, Ti, V, Cr, Mn, and Ni. For this process, stellar parameters have to be set and the only free parameter is the abundance of each element. The abundances are calculated in a region of &pm;2m &angst; around each spectral line. The line list of the neighboring lines is taken from [VALD](http://vald.astro.uu.se/~vald/php/vald.php) [@Ryabchikova:2015].

``FASMA`` is run via terminal by setting the user options in a configuration file. ``FASMA`` includes all the standard inputs for spectral synthesis along with a manual for the derivation of stellar parameters and chemical abundances. The user has to provide solely the stellar spectrum for the analysis. ``FASMA`` can be used directly for most optical surveys, such as the Gaia-ESO [@Gilmore:2012] but also for the characterization planet hosts (see the [SWEET-cat](https://www.astro.up.pt/resources/sweet-cat/)).

# Figures

![Synthetic spectrum with solar values generated with ``FASMA`` (blue) compared with observed spectrum (orange).
](img/Sun_fasma.png)

# Acknowledgements

We acknowledge contributions from V. Adibekyan, H. Tabernero, and E. Delgado-Mena. We also thank the referees for improving this work. This research made use of the Vienna Atomic Line Database operated at Uppsala University, the Institute of Astronomy RAS in Moscow, and the University of Vienna. We thank the PyAstronomy and Astropy communities.

# References
