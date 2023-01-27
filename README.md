# ZOAdaQN 

An adaptive quasi-Newton algorithm for zeroth-order stochastic optimization

![GitHub](https://img.shields.io/github/license/POptUS/ZOAdaQN)

This code was first introduced in the [paper](https://arxiv.org/abs/2109.12213) "Adaptive Sampling Quasi-Newton Methods for Zeroth-Order Stochastic 
Optimization" by Raghu Bollapragada and Stefan M. Wild.


## Getting started

Executing `Run_Experiments.m` will run the experiments and plot the 
results in the aforementioned manuscript. One has to choose 
* Name of the dataset (e.g., '15-absnormal')
* Noise variance (10^-3 or 10^-5)
* Number of random runs (e.g., 5)

When running `Run_Experiments.m`, the options given in the aforementioned manuscript are used. 
The following directories will be created:
* `Results`: Folder in which all the results are stored as .mat files
* `Plots`: Consists of all the plots showing the results

## Contents

The package consists of different folders:

* `PlotFunctions`: Consists of functions used to generate plots
                        This includes a version of Mark Schmidt's [prettyPlot](https://www.cs.ubc.ca/~schmidtm/Software/prettyPlot.html)
* `ZOAdaQNFunctions`: Consists of all the main functions used in the algorithms
                        e.g., line search, quasi-Newton techniques, variance functions, etc.
* `TestFunctions`: Contains the SG algorithm, and files given in the
    			[BenDFO](https://github.com/POptUS/BenDFO)
    			and
    			[YATSOp](https://github.com/POptUS/YATSOp)
    			repos
    			
The information about datasets used in the experiments is given in `TestFunctions/DatasetInformation.xlsx`

## Citing ZOAdaQN

Please cite the following [paper](https://arxiv.org/abs/2109.12213):

```
@article{ZOAdaQN2022,
	Author      = {Raghu Bollapragada and Stefan M. Wild},
	Title       = {Adaptive Sampling Quasi-{N}ewton Methods for Zeroth-Order Stochastic Optimization},
	Journal = {ArXiv},
	Year        = {2022},
	ArxivUrl    = {https://arxiv.org/abs/2109.12213},
	Url         = {https://arxiv.org/abs/2109.12213},
}
```

## Contributing to ZOAdaQN 

Contributions are welcome in a variety of forms; please see [CONTRIBUTING](CONTRIBUTING.rst).


## License 

All code included in ZOAdaQN is open source, with the particular form of license contained in the top-level 
subdirectories.  If such a subdirectory does not contain a LICENSE file, then it is automatically licensed 
as described in the otherwise encompassing ZOAdaQN [LICENSE](/LICENSE).  

## Resources

To seek support or report issues, e-mail:

 * ``poptus@mcs.anl.gov``

