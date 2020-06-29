# ![](https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/cheese.png) CheeseSquare | PERSONALITY

The file `IdentityDemos` is a good starting point as it includes several examples for using the code. The script shows how to:

1. Compute the identity-domains and showing their behavioral correlates
2. Use pre-computed Identity-Domains to produce a map of the personality space

To use your own data that was processed using CheeseSquare, concatenate tables provided by the behavioral analysis ([Behavior](https://github.com/AnonyMouseNeuro/CheeseSquare/tree/master/Behavior); use the 'obj.Profile.Tables.Individual' table). In the table, make sure each group gets a unique 'GroupNumer' value and each mouse assigned a unique 'MouseNumber'.

### Files

| File                | Description                                         |
| ------------------- | --------------------------------------------------- |
| `IdentityDemos.m`   | Includes a set of examples for using the code       |
| `IdentityDomains.m` | Main class for computing the IDs                    |
| `DimReduction.m`    | A collection of Dimensionality Reduction algorithms |
| `Auxiliary.m`       | Auxiliary functions                                 |
| `Color.m`           | A collection of colors and tools for colors         |
| `data_table.mat`    | Sample table                                        |
| `ID.mat`            | Precomputed IDs and Archetypes                      |

### Prerequisites

The code was tested on Matlab 2017a on both a Windows and MacOS machines, but should in principal work on earlier version of Matlab as well as other operating systems

### 