# Overall Information:

- This repository contains vibe-coded Python and C files as tools to explore the mathematical object defined in the PDF
- Claude Code, Grok, and some ChatGPT Plus were used to generate most code


# File Descriptions:

## "The Vector Sequence Framework.pdf"

> PDF paper that this repository accompanies. Beginning in number theory, this paper gives context for the following files in this repository. All graphs and data referenced and used in this paper came from output of one or more of the following Python and C files in this repository.
 
## "vector_sequence_GUI.py"
 
> GUI that takes input k, sequence type, starting x_0, ending x_0, how many v_n to check, and how many v_n to print, and calculates the associated vector sequences outlined in the PDF. It gives information for individual vector sequences including long-term behavior (zeroed, grown, cycling), which cycle it is in, and more.

## "sequence_simulator.py"
  
> GUI that allows you to visualize 1-4 vector sequence(s) evolve overtime. Given x_0 and sequence type it plays an animation of what the sequence looks like with faded data points representing the history of the vector sequence. For the general Hailstone function a $q$ must be specified

## "ZGC_grapher.c" and "ZGC_grapher.py"
  
> Two files that work in tandem. Open VS Code and in the terminal first run the C file to calculate the data points (based on input parameters including x_0_start, x_0_end, # data points, k_start, k_end, and sequence (i.e. Hailstone, Aliquot, Euler-Totient,...) and export them to a CSV file. Then run the Python file to plot them with mathplotlib. The C file (data generation) is optimized for speed and the Python file (visualization) is optimized for visual clarity. The plot then shows the k vs. Z#, G#, and C# data ploints as a line graph and in the title lists the parameters used in the C file to generate the data.

## "ZGC_animate.py"
  
> Python file that runs in the terminal to create an animation for multiple $q$ values for when $S$ is the general Hailstone function. Required parameters for the animation include $q_{start}$, $q_{end}$, $k_{start}$, $k_{end}$, $x_{start}$, $x_{end}$, and max $v_n$. Optional parameters include a $q$-step value, time per $q$ frame in the video, selection between GIF or MP4 output, number of workers computing the data, DPI, number of colors, video width, and video height. 

## "ZGC_cache.zip"

> This ZIP cache file holds saved CSV files with data for $q_{start}=1$ and $q_{end}=2500$ for $k_{start}=1$, $k_{end}=q+1.001$, $x_{start}=1$, $x_{end}=1000$, and max $v_n=1000$ for ease of reproducing 'animation.mp4'. For directions on how to create your own animations see the 'The Vector Sequence Framework.pdf' file for context and/or the 'ZGC_animate.py' file description and internal Python file comments.
