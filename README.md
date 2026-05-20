This repository contains one LaTeX file giving context for the 3 Python files. Here are descriptions for each Python file...

"balancing_k.py" - GUI that takes input k, starting x_0, ending x_0, how many v_n to check, and how many v_n to print, and 
                   calculates the associatedvector sequences outlined in the LaTeX document. It gives information for 
                   individual vector sequences including long-term behavior (zeroed, grown, cycling), which cycle it is in, 
                   and more.

"collatz_vector_multigraph_simulator_fixed.py" - GUI that allows you to visualize any given vector sequence evolve overtime. 
                                                 Given x_0 it plays an animation of what the sequence looks like with a 
                                                faded history feature to keep the full picture.

"k_vs_Z_and_G_graph_generator.py" - Allows you to see what k vs. (Z# and G#) looks like for any range of k with a specified
                                    number of data points, x_0_start, x_0_end, and max v_n check
