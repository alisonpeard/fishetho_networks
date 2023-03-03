# Network clustering on a dataset of fish welfare metrics
This repo contains code for [the DeMoS Institute](https://www.demos-institute.org) Carefish/Catch project which you can read about [here](https://www.demos-institute.org/carefish-catch). This code explores the relationships between 10 fish welfare criteria from the [FishEthoBase](https://fair-fish-database.net). For 55 fish, a weighted correlation matrix calculates the certainty-adjusted Spearman correlation between different criteria [[2]](#2). This correlation matrix is interpreted as a weighted graph and the SPONGE network clustering method of Cucuringu [[1]](#1) is applied to identify groups of criteria which appear to have positive relationships.

### Example
<img src="https://user-images.githubusercontent.com/41169293/222765606-b9224bb3-a8b6-4ce9-a4f5-8d989a418947.png" height="500" />
This is a network of Potential of a criteria being met for 55 fish. The strength of the relationships between each pair of criteria is given by the Spearman correlation between the Potential scores for each criteria across the fish. Nodes are colour coded by the communities they were assigned by SPONGE clustering. Yellow edges denote positive correlations and blue edges denote negative correlations, the thickness of an edge is proportional to the absolute value of the correlation.                                                               
                                                                                                                                      
## References
<a id="1">[1]</a> 
Cucuringu, M., Davies, P., Glielmo, A., & Tyagi, H. (2019, April). SPONGE: A generalized eigenproblem for clustering signed networks. In The 22nd International Conference on Artificial Intelligence and Statistics (pp. 1088-1098). PMLR.

<a id="2">[2]</a> 
Paul Bailey and Ahmad Emad. wCorr: Weighted Correlations, 2023. R
package version 1.9.6.


