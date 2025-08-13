# Structural connectome architecture and biological vulnerability shape cortical atrophy in cocaine use disorder



Authors: Ziteng Han, Guoyuan Yang, Tianyi Yan



### Overview

Contents include the source codes and supporting data for our manuscript.



### Codes

The main analyses were carried out step-by-step using the codes as follows.

1. **Calculate the cocaine use disorder(CUD)-related cortical thickness deviation patterns and conduct network-level statistical analyses**: S1\_Wscore\_calculation.m; S1\_Network\_atrophy\_stats.ipynb;
2. **Calculate the spatial correlation between regional cortical atrophy and the mean atrophy of structurally connected neighbours**: S2\_S3\_Region\_to\_neighbour\_atrophy\_association.ipynb; S2\_System\_level\_stats.ipynb;
3. **Calculate the spatial correlation between regional cortical atrophy and the weighted mean atrophy of structurally connected neighbours,with weights derived from structural connecome(SC)-masked biological similarity matrices**: S2\_S3\_Region\_to\_neighbour\_atrophy\_association.ipynb;
4. **Identify the CUD epicenters with data-driven ranking-based method**: S4\_Epicenter\_ranking\_identification.ipynb;
5. **Perform epicenter-transcription analysis**: S5\_AHBA\_analyses.m;
6. **Assign PLS1 genes to seven canonical cell classes**: S6\_Cell\_type\_analysis.ipynb;
7. **Explore the association between the epicenter likelihood map and the clinical sympotom-response map derived from repeated transcranial magnetic stimulation (rTMS)**: S7\_TMS\_response\_epicenter\_regression.ipynb;
8. **Generate the individualized CUD epicenter likelihood map with data-driven ranking-based method**: S8\_Individualized\_epicenter\_identification.ipynb;
9. **Explore the multivariate associations between individual epicenters and behaviors**: S9\_PLS\_input\_preparation.m; S9\_Behavioral\_PLS\_analysis.ipynb;
10. **Susceptible-Infected-Removed (SIR) agent-based model**: S10\_SIR\_modelling.ipynb.



### Open source tools

1. Regional microarray expression data were preprocessed using a recommended pipeline with the abagen toolbox (https://abagen.readthedocs.io/en/stable/);
2. The spin test was conducted with the neuromaps toolbox (nulls.alexander\_bloch, https://netneurolab.github.io/neuromaps/);
3. The two rewired null models were generated seperately using a Matlab code package (https://www.brainnetworkslab.com/coderesources) and the randmio\_und\_connected function from the Brain Connectivity Toolbox ([https://brainconn.readthedocs.io/en/latest/generated/brainconn.reference.randmio\_und\_connected.html#brainconn.reference.randmio\_und\_connected](https://brainconn.readthedocs.io/en/latest/generated/brainconn.reference.randmio_und_connected.html#brainconn.reference.randmio_und_connected));
4. Gene enrichment analysis was performed using Metascape (https://metascape.org/), an online platform for gene function annotation and enrichment analysis;
5. Visualization was carried out with Connectome Workbench (https://www.humanconnectome.org/software/get-connectome-workbench).

### Data

1. Subject lists;
2. Cortical atrophy maps;
3. Epicenter likelihood maps:
4. PLS1 genes;
5. Cell type-specific enrichment results;
6. rTMS response maps;
7. Behabioural PLS1 yloadings.



### Contact information

For more questions, please email: ziteng.han07@gmail.com or ziteng.han@bit.edu.cn

