- This repository contains one LaTeX file giving context for the 3 Python files.
- This repository contains vibe-coded Python and C files as tools to explore the mathematical object defined in the LaTeX doc
- Claude Code, Grok, and some ChatGPT Plus were used to generate most code

Here are descriptions for each file:
 
"balancing_k.py"
  GUI that takes input k, starting x_0, ending x_0, how many v_n to check, and how many v_n to print, and calculates the associated vector sequences outlined in the LaTeX document. It gives information for individual vector sequences including long-term behavior (zeroed, grown, cycling), which cycle it is in, and more.

"collatz_vector_multigraph_simulator_fixed.py"
  GUI that allows you to visualize any given vector sequence evolve overtime. Given x_0 it plays an animation of what the sequence looks like with a faded history feature to keep the full picture.

"ZGC_grapher.c" and "ZGC_grapher.py"
  Two files that work in tandem. Open VS Code and in the terminal first run the C file to calculate the data points (based on input parameters including x_0_start, x_0_end, # data points, k_start, k_end, and sequence (i.e. Hailstone, Aliquot, Euler-Totient,...) and export them to a CSV file. Then run the Python file to plot them with mathplotlib. The C file (data generation) is optimized for speed and the Python file (visualization) is optimized for visual clarity. The plot then shows the k vs. Z#, G#, and C# data ploints as a line graph and in the title lists the parameters used in the C file to generate the data.
