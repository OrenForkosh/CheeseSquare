# ![](https://github.com/OrenForkosh/CheeseSquare/raw/master/images/cheese.png) CheeseSquare


We present here all the source-code and hardware design needed to track position, behavior, and perosnalities of mice (as well as other animals). 

This code provides tools for infering personality traits from behavior. While personality traits offer considerable insight into the biological basis of individual differences, existing approaches toward understanding personality across species rely on subjective criteria and limited sets of behavioral readouts, resulting in noisy and often inconsistent outcomes. In a recent paper, we introduced a mathematical framework for studying individual differences along dimensions with maximum consistency and discriminative power. This framework was validated in mice, using data from a system for high-throughput longitudinal monitoring of group-housed male mice that yields a variety of readouts from all across an individual’s behavioral repertoire. The code is implemented in Matlab.

In our experiments, we used groups of four mice, that were marked with dyes of four different colors for identification purposes. The mice were housed in an enriched semi-naturalistic environment (see design [here](https://github.com/OrenForkosh/CheeseSquare/tree/master/Arena)) where they could move and interact freely over multiple days. Each arena contained a closed nest, two feeders, two water bottles, two ramps, an open shelter, and an S-shaped separation wall in the center. 

<p align="center">
<img src="https://github.com/OrenForkosh/CheeseSquare/raw/master/images/arena%20schematics.png" width="33%" />
</p>


All of the mice movements were automatically tracked and their behavior analyzed; with beahvioral readouts consisting of both individual (for example, locomotion, exploration and foraging patterns) and social (for example, approaches, contacts and chases). A total of 60 features per mouse per 12-h active phase were collected. For mice, the active phase is the dark time, so we used color light-sensitive cameras. In our case we used Sony's IMX174 sensor mounted on a [MANTA G-235](https://www.alliedvision.com/en/products/cameras/detail/Manta/G-235.html) by Allied Vision (although similar cameras are also available by Basler, Flir, and others).



The code is currently written in Matlab, although a Python version is also in development. 

The package is divided into these components:

- [Preprocessing](https://github.com/OrenForkosh/CheeseSquare/tree/master/Preprocessing). The code needed to prepare the video files for tracking. This process is GUI based and involves marking the bounds of the arena, automatically detecting dark/light transition, and marking the animals (open the 'Preprocessing' folder for more information)
- [Tracking](https://github.com/OrenForkosh/CheeseSquare/tree/master/Tracking). Used to track the position of the mice based on their colors.
- [Behavior](https://github.com/OrenForkosh/CheeseSquare/tree/master/Behavior). Using the tracked trajectories to infer behaviors such as chases, exploration, and more.
- [Personality](https://github.com/OrenForkosh/CheeseSquare/tree/master/Personality). Extracting traits that underly behavior which: (1) capture and represent a continuous gradient of differences between individuals of the same species and (2) tend to be stable for individuals over time.

In addition, we provide the design of the social arenas in the [Arena](https://github.com/OrenForkosh/CheeseSquare/tree/master/Arena) section.

The two additional folders [CheeseSquare](https://github.com/OrenForkosh/CheeseSquare/tree/master/CheeseSquare) and [Basics](https://github.com/OrenForkosh/CheeseSquare/tree/master/Basics) provide the base class that contains the tracking, behavior, personality, and meta data for each mouse, as well as basic utilities (accordingly).
