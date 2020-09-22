# Split_Operator_Method

# Summary

This code solve the time dependency of a Gaussian package with ´Schrödinger  Equation". For numerical results i've use the Split-Operator Method in 2 dimensions.

The README summary will be shown as seen below:

<p align="center"><img src="/tex/ca7bbcfc06e90bbbb2599b398b8b2075.svg?invert_in_darkmode&sanitize=true" align=middle width=359.22475545pt height=115.06849364999998pt/></p>

# Packages-C++
Well, for this code we use a lot of packages:
<p align="center"><img src="/tex/1b8b2b657d49ced08c0a7bce7c52ba24.svg?invert_in_darkmode&sanitize=true" align=middle width=676.89547695pt height=185.57078264999998pt/></p>

\section{Packages-Python3}
Well, we plot al data in Python3 with the fallowing packages:

<p align="center"><img src="/tex/d76a6f8a0d485bf7a425c157414b20ff.svg?invert_in_darkmode&sanitize=true" align=middle width=675.616821pt height=165.84475544999998pt/></p>

# Basic Constants and Considerations 
The code have three important pre-processing constants: BASE is about the points for discretization as we want go to reciprocal space the FFT is really much faster when the number of points is a power of two; is the usual hbar in quantum mechanics; mass is the mass of the gaussian package.

All functions are using reference passing, not a object passing (for speed reasons) also in this field of thinking, in the code I've avoid to use if's statements like:

FIGURA

\section{Algorithm}

Let's think in a PDE like:

<p align="center"><img src="/tex/87a5a7ba313bd16eac2ab6afd6174237.svg?invert_in_darkmode&sanitize=true" align=middle width=490.9215762pt height=40.11819404999999pt/></p>

We can write this Hamiltonian in two parts, a real-space part and a reciprocal parte (´r" and ´k" respectively) like ´<img src="/tex/0a2aa24f45ff52f5fef3e9a7d9bc1b20.svg?invert_in_darkmode&sanitize=true" align=middle width=98.88298694999997pt height=22.465723500000017pt/>" where ´<img src="/tex/83a21f7a8c4d5ab4a44f6653c07462b5.svg?invert_in_darkmode&sanitize=true" align=middle width=98.86921274999999pt height=33.45973289999998pt/>" and ´<img src="/tex/1a4b9cc06116da2789e3d8044d57acde.svg?invert_in_darkmode&sanitize=true" align=middle width=76.6741866pt height=24.65753399999998pt/>". So, taking a initial condition for ´<img src="/tex/1c899e1c767eb4eac89facb5d1f2cb0d.svg?invert_in_darkmode&sanitize=true" align=middle width=36.07293689999999pt height=21.18721440000001pt/>" (In my code a gaussian package):

<p align="center"><img src="/tex/584b83fbd098870e9536942b39dfe321.svg?invert_in_darkmode&sanitize=true" align=middle width=524.38206975pt height=23.6529876pt/></p>

Taking ´<img src="/tex/ef16ee29cc72048df1bd5d3e88675505.svg?invert_in_darkmode&sanitize=true" align=middle width=44.62890134999999pt height=22.831056599999986pt/>' and using the Baker-Campbell-Housdorff formula:
<p align="center"><img src="/tex/4209ba3c7584d62744d15ddc5370f3b6.svg?invert_in_darkmode&sanitize=true" align=middle width=533.84799435pt height=39.452455349999994pt/></p>

The above formula have an error of `<img src="/tex/d768ed002d227252dd0a0adc27b91109.svg?invert_in_darkmode&sanitize=true" align=middle width=21.04460819999999pt height=26.76175259999998pt/>' order, in order to make that smaller we cand split the space step into two half-steps, what leave us to a ´<img src="/tex/d0e485049ae83c290f34b353d3ec90d7.svg?invert_in_darkmode&sanitize=true" align=middle width=21.04460819999999pt height=26.76175259999998pt/>" error:

<p align="center"><img src="/tex/a1ef69b6a25fc82a0e0e2f12e4c5b434.svg?invert_in_darkmode&sanitize=true" align=middle width=547.3717128pt height=29.77187565pt/></p>

Is \textbf{important} make that obvious, here `<img src="/tex/b82909b02a40877fd70c2350d229b92a.svg?invert_in_darkmode&sanitize=true" align=middle width=19.30029914999999pt height=27.6567522pt/>' where `H' is a matrix we are not expanding in terms of taylor series, but applying exponential of each element of the matrix.

We can easily solve that with Fourier and Inverse Fourier trasform (below <img src="/tex/5815709515a245ef0734c68053c5128b.svg?invert_in_darkmode&sanitize=true" align=middle width=91.17132089999998pt height=34.063933200000015pt/> and <img src="/tex/ac2cd6f31708991900e2154f7b667b66.svg?invert_in_darkmode&sanitize=true" align=middle width=92.54120204999998pt height=36.52973610000002pt/> ):
<p align="center"><img src="/tex/04ca66e1234a7b36fcbfbd723d837c3d.svg?invert_in_darkmode&sanitize=true" align=middle width=555.07112265pt height=39.452455349999994pt/></p>

# How to Use?
Well, first run the command `./pre\_requisites' to install all pre-requisites and create all necessary folders.\newline

To run the ´.cpp' file you can execute ´make && ./BINARY'  or just run `g++ schr\_main.cpp -larmadillo -fopenmp -lm -lfftw3 -llapacke -Ofast -lblas -o BINARY && ./BINARY'.\newline

For the python3 program you just need to write `python3 plotar.py'.\newline

