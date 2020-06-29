# ![](https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/cheese.png) CheeseSquare | PREPROCESSING

The preprocessing is done using a graphical wizard. To start the wizrd run the 'CheeseInit' command from the Matlab command line.

Note that for optimal use the wizard assumes that all files have the following name scheme:

 `(experiment name).exp(group number).day(day number).cam(camera number).avi`

For example, for the personality baseline test and group number 1 which we recorded four days the file names are: (a) Personality.exp0001.day01.cam01.avi, (b) Personality.exp0001.day02.cam01.avi, (c) Personality.exp0001.day03.cam01.avi, and (d) Personality.exp0001.day04.cam01.avi.



The wizards is made of the following stages:

##### 1. Init

Allows you to choose the video file to process. Fill the filename in the 'source file' field (or choose it by pressing the '...' button), and then press 'Load'. After the preprocessing is done you can use the 'Segment' button to run segmentation on the current frame (for testing purposes).

<center>
  <img src="https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/preprocessing/init.png" alt="drawing" width="70%"/>
</center>



##### 2. Scale

Used to defined the boundries of the arena (the base of the walls). Note to fill in the correct width and height in centimeters (if you use our design it should be 60x60 cms). This step is important so that the algorithm would be able to translate pixels to actuall distance (in centimeters) when measuring behaviors such as travel distance and speed.

<center>
  <img src="https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/preprocessing/scale.png" alt="drawing" width="70%"/>
</center>



##### 3. Background

Automatically extracts a background image by computing the median of a set of random frames (the number is set by the 'Number of frames' parameter). Ensure that none of the mice is visible following this stage (this usually happens if a mouse remained static for a long period).

This step of the wizrd also determines the onset and end of the dark phase by comparing the average brightness to a 'Darkness' value.

<center>
  <img src="https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/preprocessing/background.png" alt="drawing" width="70%"/>
</center>



##### 4. Label Mice

Gives the color segmentation algorithm several examples of each mouse. Mice that were detected by the algorithm would be enclosed by a white curve. To mark a mouse enter its number (by pressing the number on the keyboard) and click on its image. Avoid marking a mouse when it is in contact with another mouse (when the white curve surrounds both mice). We usually mark 20 utterances of each mouse. To get a new image press the 'New frame' button.

<center>
  <img src="https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/preprocessing/colors.png" alt="drawing" width="70%"/>
</center>

##### 5. Mark Regions

Mark various regions-of-interest (ROIs) within the arena. A ROI can be defined a nest when mice become hidden once within this region. In addition, when the ROI does not appear in the experiment it can be marked as non available.

<center>
  <img src="https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/preprocessing/regions.png" alt="drawing" width="70%"/>
</center>

##### 6. Save

The final step of the wizard creates a '.obj.mat' file with all the preprocessing data for the video file and all video recordings of the same group on other days.

<center>
  <img src="https://github.com/AnonyMouseNeuro/CheeseSquare/raw/master/images/preprocessing/save.png" alt="drawing" width="70%"/>
</center>