Code repository for Mousley, V. L., Contini, L., Re, R., Pinti, P. & Mareschal, D. (under review). Analytical pipeline optimisation in developmental fNIRS hyperscanning data: Neural coherence between preschoolers collaborating with their mothers. 

Please contact Victoria Mousley with any questions (v.mousley@bbk.ac.uk).

**Requirements**

_Software_

  - MATLAB 2020B (installation: https://uk.mathworks.com/help/install/install-products.html)
  - RStudio 4.2.2 (installation: https://rstudio.com/products/rstudio/download/)

_Toolboxes/Functions_

  - Homer2 Toolbox (https://homer-fnirs.org/download/). Huppert, T. J., Diamond, S. G., Franceschini, M. A., & Boas, D. A. (2009). HomER: A review of time-series analysis methods for near-infrared spectroscopy of the brain. Applied Optics, 48(10), D280–D298. https://doi.org/10.1364/ao.48.00d280
  - Cross Wavelet and Wavelet Coherence Toolbox (https://grinsted.github.io/wavelet-coherence/). Grinsted, A., Moore, J. C., & Jevrejeva, S. (2004). Application of the cross wavelet transform and wavelet coherence to geophysical time series. Nonlinear Processes in Geophysics, 11, 561–566. https://doi.org/1607-7946/npg/2004-11-561
  - Short-separation regression functions: Abdalmalak, A. et al. (2022). Effects of systemic physiology on mapping resting-state networks using functional near-infrared spectroscopy. Frontiers in Neuroscience, 16, 803297. https://doi.org/10.3389/fnins.2022.803297. Questions about these functions ('AdjustTemporalShift_fnirs_course.m', 'ApplySSC', and 'PhysiologyRegression_GLM_fnirs_course.m') should be directed to the original authors: Professor Rickson Mesquita (r.c.mesquita@bham.ac.uk) and Dr Sergio Luiz Novi Jr (novisl@ifi.unicamp.br).
  
**Example data**

After cloning the repository, you can run all scripts with the example data provided (see /example_data/) by adding the path at the beginning of the scripts. Example fNIRS dataset forthcoming. The WTC pipeline and analysis scripts also incorporate input from Excel files, such as excluded channels identified manually, background information about participants (e.g., children's ages), and behavioural data (e.g., collaborative task performance, parent-report questionnaire data). Templates for these file structures are also provided.

**Structure**

***WTC Pipeline***
1. config.m = study-specific parameters to be set before running pipeline
2. runWTC.m = main script
  - option to run with or without phase-scrambled pseudodyad analysis (line 37);
  - option to run with or without short-separation channel regression (line 39);
    - calculates channel-wise WTC;
    - calculates ROI-wise WTC;
    - averages across conditions of interest;
    - saves final data as labelled .struct in analysis path
4. exportData.m = collates outputted .m structs from pipeline for analysis in R

***Statistical Analysis***
1. analysis.R = statistical analysis script
  - unpacks MATLAB structs for manipulation in R;
  - runs analyses for each aim presented in Mousley et al. (under review)
