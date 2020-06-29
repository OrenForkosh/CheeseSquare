# ![](https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/cheese.png) CheeseSquare | TRACKING

Note! To run the tracking algorithm, make sure to run the [Preprocessing](https://github.com/AnonyMouseNeuro/CheeseSquare/tree/master/Preprocessing) step first.

Tracking is initiated by running the 'CheeseTrackRun' command (in Matlab) on the preprocessed video file. For example,

```matlab
>> CheeseTrackRun('Videos/Personality.exp0001.day02.cam01.avi')
```

The tracking algorithm is divided into two main steps:

1. Frame segmentation - The algorithm uses image segmentation to find the possible locations of each mouse in every frame. Note that this step is divided across cpus so having a machine with many cores can speed this step of the algorithm considerably.
2. Path tracking - finding the most likely path trajectory for each mouse based on the blobs that were detected by the segmentation step.

