# Overall Information:

- This repository contains vibe-coded Python and C files as tools to explore the mathematical object defined in the PDF
- Claude Code, Grok, and some ChatGPT Plus were used to generate most code
- The ZIP cache file holds saved CSV files with data for $q_{start}=1$ and $q_{end}=2500$ for $k_{start}=1$, $k_{end}=q+1.001$, $x_0$, ${x_0}_{end}=1000$, max $v_n=1000$ for ease of reproducing the MP4 animation included. For directions on how to create your own animations see the PDF file for context, file description below, and comments in the Python file



# File Descriptions:
 
## "vector_sequence_GUI.py"
 
> GUI that takes input k, sequence type, starting x_0, ending x_0, how many v_n to check, and how many v_n to print, and calculates the associated vector sequences outlined in the PDF. It gives information for individual vector sequences including long-term behavior (zeroed, grown, cycling), which cycle it is in, and more.

## "sequence_simulator.py"
  
> GUI that allows you to visualize 1-4 vector sequence(s) evolve overtime. Given x_0 and sequence type it plays an animation of what the sequence looks like with faded data points representing the history of the vector sequence

## "ZGC_grapher.c" and "ZGC_grapher.py"
  
> Two files that work in tandem. Open VS Code and in the terminal first run the C file to calculate the data points (based on input parameters including x_0_start, x_0_end, # data points, k_start, k_end, and sequence (i.e. Hailstone, Aliquot, Euler-Totient,...) and export them to a CSV file. Then run the Python file to plot them with mathplotlib. The C file (data generation) is optimized for speed and the Python file (visualization) is optimized for visual clarity. The plot then shows the k vs. Z#, G#, and C# data ploints as a line graph and in the title lists the parameters used in the C file to generate the data.
