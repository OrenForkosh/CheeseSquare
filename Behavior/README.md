# ![](https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/cheese.png) CheeseSquare | BEHAVIOR

Note! make sure to run the [Tracking](https://github.com/AnonyMouseNeuro/CheeseSquare/tree/master/Tracking) algorithm first.

The behavioral analysis computes various behavioral and social readouts for each mouse. To run the analysis on a tracked video file, run the following command in Matlab

```matlab
>> obj = CheeseNewProfiler('Videos/Personality.exp0001.day01.cam01.avi')
```

the results of the behavioral analysis are then stored in 'obj.Profile'.

The algorithm measures behaviors such as chases, approaches, exploration, time in the feeder, and so on. It also evaluates the social structure of the group such as the dominance rank of each mouse. For example, 

```matlab
obj.Profile.Tables.Individual
```

 contains the following table (we show here only a subset of the behaviors for illustration):

| **Day** | **GroupID** | **GroupNumber** | **GroupType**   | **MouseID** | **MouseNumber** | **SocialRank** | **NormalizedDavidScore** | **FractionOfTimeOutside** | **VisitsOutsideRate** | **ForagingCorrelation** | **ContactRate**  | **MedianContactDuration** | **FractionOfChasesPerContact** | **FractionOfEscapesPerContact** | **FractionOfFollowsPerContact** | **FractionOfBeingFollowedPerContact** |
| ------- | ----------- | --------------- | --------------- | ----------- | --------------- | -------------- | ------------------------ | ------------------------- | --------------------- | ----------------------- | ---------------- | ------------------------- | ------------------------------ | ------------------------------- | ------------------------------- | ------------------------------------- |
| **1**   | 1           | 1               | personality_mix | 1           | 1               | 2              | 1.50208333333333         | 0.694917284686662         | 37.1926383649726      | 0.285678433252205       | 89.7638058516418 | 4.16                      | 0.015828677839851              | 0.015828677839851               | 0.0800744878957169              | 0.0577281191806331                    |
| **1**   | 1           | 1               | personality_mix | 2           | 2               | 2              | 0.999305555555556        | 0.607259668692992         | 43.9625343370238      | 0.364166440632409       | 74.3852767299452 | 4                         | 0.00224719101123596            | 0.0168539325842697              | 0.0348314606741573              | 0.0662921348314607                    |
| **1**   | 1           | 1               | personality_mix | 3           | 3               | 3              | 0.898611111111111        | 0.666173550787221         | 35.1031642995247      | 0.286989493177304       | 84.4147522440951 | 3.92                      | 0.00099009900990099            | 0.0138613861386139              | 0.0198019801980198              | 0.0653465346534653                    |
| **1**   | 1           | 1               | personality_mix | 4           | 4               | 1              | 2.6                      | 0.768517847823139         | 44.6311660379671      | 0.357746340707672       | 97.45307041249   | 4                         | 0.0222984562607204             | 0                               | 0.0668953687821612              | 0.0240137221269297                    |



##### Behavioral readouts

For the personality analysis, we collected a total of 60 different readouts for each mouse on each day. Due to the linearity of LDA used to compute the personality traits, some of the behaviors we measured were computed with several different normalizations. The most common normalizations we used were: the total time in the arena (12 hours; abbreviated 'total'), the time outside the nest ('outside'), and for interactions the total number of contacts ('contacts').

Pairwise:

·    Time outside [1]: Fraction of time that the mouse spends outside of the nest. Normalizations: Total time (%).

·    Frequency of visits outside [2]: The rate at which the mouse exits the nest. Normalizations: Total time (1/hour).

·    Foraging correlation [3]: The correlation between the times that the mouse is outside the nest and the times that another mouse is outside the nest, averaged over all mice. For example, the foraging correlation between two mice would equal one if the mouse is always outside the nest when the other mouse is outside, and also enters the nest whenever the other mouse enters the nest. The correlation would be -1 whenever the mouse is outside, the other mouse is inside the nest. Normalizations: [3] none (au).

·    Contact rate [4, 5]: Number of contacts the mouse had. A contact is defined as two mice being less than 10 cm apart while both are outside the nest. Normalizations: [4] Total time (1/hour), [5] Time outside (1/hour).

·    Time in contact [6]: Fraction of time that a mouse is in contact with other mice while outside the nest. Normalizations: [6] Time outside (1/hour).

·    Median\Mean contact duration [7, 8]: Median or mean duration of contacts. The contact duration does not include the times when the mouse approached, moved away from, or chased the other mouse. Normalizations: [7, 8] none (sec)

·    Follow [12, 18, 24]: A follow is a contact that ended with one mouse going after another mouse until disengagement. Follows can be either aggressive (chases) or non-aggressive. Normalizations: [12] Number of contacts (au), [18] Time outside (1/hour), [24] Total time (1/hour).

·    Being-followed [13, 19, 25]: Number of times a mouse is followed at the end of a contact. It can be either in an aggressive (chases) or non-aggressive manner. Normalizations: [13] Number of contacts (au), [19] Time outside (1/hour), [25] Total time (1/hour).

·    Chase [10, 16, 22]: Chases are interactions that ended with the mouse pursuing another mouse in an aggressive manner. Aggressiveness was determined using a classifier that was trained on labeled samples (see methods). Normalizations: [10] Number of contacts (au), [16] Time outside (1/hour), [22] Total time (1/hour).

·    Escape [11, 17, 23]: Number of time that the mouse was aggressively chased by another mouse. Normalizations: [11] Number of contacts (au), [17] Time outside (1/hour), [23] Total time (1/hour).

·    Non-aggressive follow [14, 20, 26]: Number of times the mouse has followed another mouse at the end of a contact in a non-aggressive way. Normalizations: [14] Number of contacts (au), [20] Time outside (1/hour), [26] Total time (1/hour).

·    Non-aggressively being-followed [15, 21, 27]: Number of times the mouse was followed by another mouse at the end of a contact in a non-aggressive way. Normalizations: [15] Number of contacts (au), [21] Time outside (1/hour), [27] Total time (1/hour)

·    Approach [28, 29, 30, 31, 32]: An approach is a directed movement of the mouse towards another mouse that ends in contact. Not all interactions necessarily start with an approach, while others might start mutually with both mice approaching each other. Normalizations: [28] none (au), [29] Time outside (1/hour), [31] Number of contacts (au), [32] total (1/hour), [30] Time outside with one or more mice (1/hour).

·    Being-approached [33, 34, 35]: Number of times the mouse was approached by another mouse. Normalizations: [33] Number of contacts (au), [34] Time outside (1/hour), [35] none (au).

·    Approach-escape [36]: Fraction of contacts in which the mouse initiated the contact and ended up being chased. Normalizations: [36] Number of aggressive contacts (au).

·    Difference between approaches and chases [9]: The total number of chases is subtracted from the total number of approaches. Normalizations: [9] none (au).

Individual:

·    ROI exploration [37, 38]: Quantifies the amount of exploration the mouse is doing. Measured as the entropy of the probability of being in each of the 10 regions-of-interest (ROIs). Mice that spend the same amount of time in all regions will get the highest score, while mice that spend all their time in a single ROI will be scored zero. When normalized to the time outside, the computation of the entropy differed also by ignoring the probability of being inside the nest. Normalizations: [37] none (bits), [38] Time outside (bits/hour).

·    Grid exploration [59]: Quantifies the amount of exploration the mouse is doing. Analogously to ‘ROI Exploration’, grid exploration was determined using entropy, however, instead of looking at the ROIs, we divided the arena into a 6x6 grid (10 cm by 10 cm; a total of 36 possible locations). Normalizations: [59] none (bits).

·    Predictability [60]: Measures how predictable the paths that the mouse takes as the mutual information between its current and previous location in the arena. For that, the arena was divided into a 6x6 grid (10 cm by 10 cm; a total of 36 possible locations), and for each cell we computed the probabilities of it moving to any of the adjacent cells. Normalizations: [60] none (bits).

·    Distance [58]: The total distance traveled by the mouse while outside the nest. To smooth the tracking, the mice locations were sampled once every second. Normalizations: [58] none (m).

·    Median\Mean speed [54, 55]: Median or mean speed while outside the nest. To smooth the computation of the speed, the locations of mice were sampled once every second. Normalizations: [54, 55] none (m/sec).

·    Tangential velocity [56]: The tangential component of the speed, or the part of speed perpendicular to the previous direction of movement. Normalizations: [56] none (m/sec).

·    Angular velocity [57]: The rate of change in the direction of the mouse. Normalizations: [57] none (rad/sec).

·    Food or water [39, 40]: Time spent in the feeder or water zones. Normalizations: [39] Total time (au), [40] Time outside (au)

·    Food [41]: Time spent next to the feeders. Normalizations: [41] Time outside (au).

·    Water [42]: Time spent next to the water bottles. Normalizations: [42] Time outside (au).

·    Feeder preference [43]: Time spent at (or near) the feeder adjacent to the nest (feeder 1) relative to the further-away feeder (feeder 2). Normalizations: [43] none (au).

·    Water preference [44]: Time spent near the water bottle adjacent to the nest (water 1) relative to the further-away water bottle (water 2). Normalizations: [44] none (au).

·    Elevated area [45, 46]: Time spent on an elevated object in the arena: ramps or block. Normalizations: [45] Total time (au), [46] Time outside (au).

·    Open area [47]: Time spent in the open area (outside of the nest and any of the ROIs). Normalizations: [47] Time outside (au).

·    Shelter [48]: Time spent in the shelter, which is a box without a roof, i.e. closed only on its sides. Normalizations: [48] Time outside (au).

·    Ramps [49]: Time spent on the elevated ramps. Normalizations: [49] outside (au).

·    S-wall [50]: Time spent in the S-wall. Normalizations: [50] Time outside (au).

·    Distance from walls [51]: Average distance from the walls while in the open area. Normalizations: [51] none (cm).

·    Distance from nest [52]: Average distance from the nest (while outside of the nest). Normalizations: [52] none (cm).

·    Alone outside [53]: Fraction of time the mouse is outside while all other mice are in the nest. Normalizations: [53] Total time (au).