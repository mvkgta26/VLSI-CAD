**MAIN: Source.cpp** <br>
**SAMPLE INPUTS FOLDER: INPUT NETLISTS SAMPLES/** <br>
**FINAL SOLUTION: QP3Solution.txt**
# VLSI CAD PLACER: QUADRATIC WIRELENGTH MODEL
A C++ program for an analytical-­style placer that takes in a NETLIST as INPUT (in the form of .txt file). It uses Quadratic Wirelength model and large matrix solves to optimize the placement of very large number of gates on a chip surface, so that the wires/nets that would later be routing these gates will be short. This placer will also use a single round of recursive partitioning strategy to deal with the fact that the gates will be clustered in overlapping, nonphysical ways, after the first quadratic placement solution. The output will be given as a list of all the gates in the netlist along with their coordinates on abstract chip surface (on a .txt file). Netlist for some basic industrial benchmarks were run on this program.
