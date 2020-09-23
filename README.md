# Split_Operator_Method

# Summary

This code solves the time dependency of a Gaussian package by using "Schrödinger  Equation". For numerical results, I've used the Split-Operator Method in 2 dimensions.

The README summary will be shown as seen below:

<p align="center"><img src="/tex/ca7bbcfc06e90bbbb2599b398b8b2075.svg?invert_in_darkmode&sanitize=true" align=middle width=359.22475545pt height=115.06849364999998pt/></p>

# Packages-C++
Well, for this code we've used a lot of packages:
<p align="center"><img src="/tex/26ba1b35fdd4f0f6fec154645331119c.svg?invert_in_darkmode&sanitize=true" align=middle width=675.61699095pt height=205.29680985pt/></p>

# Packages-Python3
Well, we plot all data in Python3 by the use of the following packages:

<p align="center"><img src="/tex/5285939929b9f406f3e79a4fed1ea353.svg?invert_in_darkmode&sanitize=true" align=middle width=675.61703715pt height=185.57078264999998pt/></p>

# Basic Constants and Considerations 
The code has three important pre-processing constants: 
BASE is used for the space discretization. Therefore, the number of points must be the power of two(<img src="/tex/ca5b0a2028e3571d911e4fd5ea301a66.svg?invert_in_darkmode&sanitize=true" align=middle width=53.262831599999984pt height=22.465723500000017pt/>), once we want to use the "divide and conquest" algorithm; hbar is the usual hbar <img src="/tex/7fb5a0da6f3b27af9ce94b52af075052.svg?invert_in_darkmode&sanitize=true" align=middle width=33.66627329999999pt height=28.92634470000001pt/> in quantum mechanics; mass is the mass of the gaussian package.

All functions are using a reference passing in order to improve code performance. Also, with the "as fast as possible" philosophy, I have to avoid the use of "if's", like:

![Comaparision](comp.png)

# Algorithm

Let's think in a PDE like:

<p align="center"><img src="/tex/3eef763b3891cc1b33c590cb78f05974.svg?invert_in_darkmode&sanitize=true" align=middle width=281.56922145pt height=40.11819404999999pt/></p>

We can write this Hamiltonian in two parts, a real-space part and a reciprocal part ('r" and 'k" respectively) like "<img src="/tex/0a2aa24f45ff52f5fef3e9a7d9bc1b20.svg?invert_in_darkmode&sanitize=true" align=middle width=98.88298694999997pt height=22.465723500000017pt/>" where "<img src="/tex/83a21f7a8c4d5ab4a44f6653c07462b5.svg?invert_in_darkmode&sanitize=true" align=middle width=98.86921274999999pt height=33.45973289999998pt/>" and '<img src="/tex/1a4b9cc06116da2789e3d8044d57acde.svg?invert_in_darkmode&sanitize=true" align=middle width=76.6741866pt height=24.65753399999998pt/>". So, taking a initial condition for "<img src="/tex/1c899e1c767eb4eac89facb5d1f2cb0d.svg?invert_in_darkmode&sanitize=true" align=middle width=36.07293689999999pt height=21.18721440000001pt/>" (In my code a gaussian package):

<p align="center"><img src="/tex/92b67306e0e9e6e2aa6c52fa9066970b.svg?invert_in_darkmode&sanitize=true" align=middle width=348.49030095pt height=23.6529876pt/></p>

Taking "<img src="/tex/ef16ee29cc72048df1bd5d3e88675505.svg?invert_in_darkmode&sanitize=true" align=middle width=44.62890134999999pt height=22.831056599999986pt/>" and using the Baker-Campbell-Housdorff formula:
<p align="center"><img src="/tex/82209986de72a2fa350433f89720bf10.svg?invert_in_darkmode&sanitize=true" align=middle width=367.42210725pt height=39.452455349999994pt/></p>

The above formula has an error of "<img src="/tex/2ad4fa62e890e7471af78fd2eebd4a2d.svg?invert_in_darkmode&sanitize=true" align=middle width=48.195394499999985pt height=26.76175259999998pt/>". Aiming at making it smaller, we can split the space step into two half-steps, which lead us to a "<img src="/tex/b876e42142198eef54b7bbc7ae63a804.svg?invert_in_darkmode&sanitize=true" align=middle width=48.195394499999985pt height=26.76175259999998pt/>" situation:

<p align="center"><img src="/tex/584808948b48b262c9bb2755389de9fc.svg?invert_in_darkmode&sanitize=true" align=middle width=394.4696052pt height=29.77187565pt/></p>

In this "<img src="/tex/b82909b02a40877fd70c2350d229b92a.svg?invert_in_darkmode&sanitize=true" align=middle width=19.30029914999999pt height=27.6567522pt/>" ,where 'H' is a matrix, we are not expanding by using Taylor's series, but we are applying the exponential to each element of the matrix.

We can easily solve the above equation with Fourier and Inverse Fourier trasform (below <img src="/tex/5815709515a245ef0734c68053c5128b.svg?invert_in_darkmode&sanitize=true" align=middle width=91.17132089999998pt height=34.063933200000015pt/> and <img src="/tex/ac2cd6f31708991900e2154f7b667b66.svg?invert_in_darkmode&sanitize=true" align=middle width=92.54120204999998pt height=36.52973610000002pt/> ):
<p align="center"><img src="/tex/a0e4525bc93f3d59716e6ef470d05104.svg?invert_in_darkmode&sanitize=true" align=middle width=409.86834239999996pt height=39.452455349999994pt/></p>

# How to Use?
Well, first we must run the command "./pre\_requisites" to install all pre-requisites and create all necessary folders.

To run the '.cpp' file, you can execute "make && ./BINARY" (to compile and execute)  or just run the fallowing comand 'g++ schr\_main.cpp -larmadillo -fopenmp -lm -lfftw3 -llapacke -Ofast -lblas -o BINARY && ./BINARY' (again to compile and execute).

For the python3 program, you just need to write and run 'python3 plotar.py'.

# Results:

Some results for you have an idea about how plots are make:

A wall potential:
![Comaparision](figures/Pacote_500s_0.4.gif)
![Comaparision](figures/3dbarrier.gif)

Almost infinite one:


