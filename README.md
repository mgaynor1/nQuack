# nQuack

<img align="right" src="man/figures/nQuack.png" width=200> 


**Michelle L. Gaynor, Jacob B. Landis, Tim K. O’Connor, Robert G. Laport, Jeff J. Doyle, Douglas E. Soltis, José Miguel Ponciano, and Pamela S. Soltis**  


[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) 

## Overview
nQuack is a modified statistical framework to predict ploidy level based on sequence data. We build upon [Weib et al., 2018](https://doi.org/10.1186/s12859-018-2128-z) Gaussian Mixture Model approach to estimate ploidy level, which was originally written as [a C executable](https://github.com/clwgg/nQuire). 

## More on nQuack
Here we provided expanded tools and implementations to improve site-based heterozygosity inferences of ploidal level. 

nQuack provides data preparation guidance and tools to decrease noise in input data. These include a maximum sequence coverage quantile filter and sequence error-based filter, to remove biallelic sites that are likely not representative of copy number variance in the nuclear genome. We also consider only the frequency of allele A or B at each site, instead of both, as found in other methods. To learn more about best practices, see our [Data Preparation](https://mlgaynor.com/nQuack/articles/DataPreparation.html) guide.

Our model improves upon the nQuire framework by extending it to higher ploidal levels (pentaploid and hexaploid), correcting the augmented likelihood calculation, implementing more suitable distribution, and allowing additional ‘fixed’ models. We also decrease model selection errors by relying on BIC rather than likelihood ratio tests. To learn more about these methods, see our [Model Options](https://mlgaynor.com/nQuack/articles/ModelOptions.html) guide.

We provide 32 ways to estimates likelihood of a mixture of models with the expectation maximization algorithm ([see more here](https://mlgaynor.com/nQuack/articles/ModelOptions.html)) - 8 expectation maximization implementations with 4 model types each. In total, nQuack offers 128 models.

## Evaluation of nQuack  

To examine the utility of this method, we examined 513,792 models based on both simulated and real samples. Before using this method, we suggest you read our manuscript and consider the many limitations to a pattern-based approach for determining ploidal level.  


## Installation  

```
install.packages("devtools")
devtools::install_github("mgaynor1/nQuack")
```

### Warning: samtools must be local!   

If you are working on your personal computer, make sure samtools is installed and callable as "samtools" via terminal. If you are working on a cluster, you may need to symbolically-link samtools locally. Though the location of install may differ, here is how I make samtools callable locally on UF's amazing [HiPerGator](https://www.rc.ufl.edu/about/hipergator/) slurm cluster:    


```
mkdir bin
cd bin
ln -s /apps/samtools/1.15/bin/samtools samtools
```


For implementation, see our [Basic Example](https://mlgaynor.com/nQuack/articles/BasicExample.html) article.

## References 

Gaynor ML, Landis JB, O’Connor TK, Laport RG, Doyle JJ, Soltis DE, Ponciano JM, and Soltis PS. 2024. nQuack: An R package for predicting ploidy level from sequence data using site-based heterozygosity. *Applications in Plant Sciences* 12(4):e11606. [doi: 10.1002/aps3.11606](https://www.doi.org/10.1002/aps3.11606)


## Up Next:
 
- If you have sequence data with known plodial level for a mixed-ploidy system, let us know. We would love to collaborate with you. To be included in v2.0, please send me an email at shellyleegaynor at gmail.  



