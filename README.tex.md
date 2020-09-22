# Split_Operator_Method

# Summary

This code solve the time dependency of a Gaussian package with ´Schrödinger  Equation". For numerical results i've use the Split-Operator Method in 2 dimensions.

The README summary will be shown as seen below:

\begin{itemize}
    \item Packages (Python3 and C++)
    \item Basic Constants and Important Considerations;
    \item Analysis of Algorithm and Method;
    \item How to use all possible files (.sh, .cpp and .py)
\end{itemize}

# Packages-C++
Well, for this code we use a lot of packages:
\begin{itemize}
    \item Armadillo: is the mainly library in my code, the front-end package for the hole algebra of the system (front-end for lapacke,cblas, fftw3 and openMP);
    \item OpenMP: parallel library for all for's statements in the code, it improve the time execution a lot;
    \item Chrono: time meansurement of execution;
    \item FFTW3: fast fourier transform, it's not used explicitly, but armadillo use in the back-end;
    \item Lapacke and Cblas: algebra libraries for armadillo back-end;
\end{itemize}

\section{Packages-Python3}
Well, we plot al data in Python3 with the fallowing packages:

\begin{itemize}
    \item Seaborn: used for the heatmap plot;
    \item Matplotlib: used in back-end in Seaborn and necessary to save images;
    \item Numpy:  to load the data in Python3;
    \item Joblib: to run the image generation in parallel, important because it uses a lot of time;
    \item FFMPEG: is not a Python3 library, but inside the Python3 code we call a command to generate the video from plots;
\end{itemize}

# Basic Constants and Considerations 
The code have three important pre-processing constants: BASE is about the points for discretization as we want go to reciprocal space the FFT is really much faster when the number of points is a power of two; is the usual hbar in quantum mechanics; mass is the mass of the gaussian package.

All functions are using reference passing, not a object passing (for speed reasons) also in this field of thinking, in the code I've avoid to use if's statements like:

![Comaparision](pictures_to_readme/comp.png)

# Algorithm

Let's think in a PDE like:

\begin{equation}
    i\hbar \frac{\partial \psi(\textbf{r},t)}{dt} = \left(-\frac{\hbar^2}{2m} \nabla^2 + V(\textbf{r}) \right)\psi(\textbf{r},t)
\end{equation}

We can write this Hamiltonian in two parts, a real-space part and a reciprocal parte ('r" and 'k" respectively) like '$H = H_k + H_r$" where '$H_k = -\frac{\hbar^2}{2m} \nabla^2$" and '$H_r = V(\textbf{r}) $". So, taking a initial condition for '$t=0$" (In my code a gaussian package):

\begin{equation}
    \psi(\textbf{r},t+dt) = e^{-\frac{iHdt}{\hbar}}\psi(\textbf{r},t) = e^{-\frac{i(H_k+H_r)dt}{\hbar}}\psi(\textbf{r},t) 
\end{equation}

Taking '$dt \approx 0$' and using the Baker-Campbell-Housdorff formula:
\begin{equation}
    \psi(\textbf{r},t+dt) = \left( e^{\frac{-iH_rdt}{\hbar}}e^{\frac{-iH_kdt}{\hbar}}e^{\frac{-[iH_r,H_K]dt^2}{2}} \right) \psi(\textbf{r},t)
\end{equation}

The above formula have an error of '$dt^2$' order, in order to make that smaller we cand split the space step into two half-steps, what leave us to a '$dt^3$" error:

\begin{equation}
    \psi(\textbf{r},t+dt) = \left( e^{\frac{-iH_rdt}{\hbar}}e^{\frac{-iH_kdt}{2\hbar}}e^{\frac{-iH_rdt}{2 \hbar}} \right) \psi(\textbf{r},t) +\mathcal{O}(dt^3)
\end{equation}

Is \textbf{important} make that obvious, here '$e^{H}$' where 'H' is a matrix we are not expanding in terms of taylor series, but applying exponential of each element of the matrix.

We can easily solve that with Fourier and Inverse Fourier trasform (below $U_r = e^{-\frac{iH_rdt}{\hbar}}$ and $U_k = e^{-\frac{iH_kdt}{\hbar}}$ ):
\begin{equation}
    \psi(\textbf{r},t) = 
    \left\{
    U_r\left(\frac{dt}{2}\right) \mathcal{F}^-1
    \left[
    U_k(dt) \mathcal{F} 
    \left( U_r\left(\frac{dt}{2} \right)\psi(\textbf{r},t)
    \right)
    \right]
    \right\}
\end{equation}

# How to Use?
Well, first run the command './pre\_requisites' to install all pre-requisites and create all necessary folders.\newline

To run the '.cpp' file you can execute 'make && ./BINARY'  or just run 'g++ schr\_main.cpp -larmadillo -fopenmp -lm -lfftw3 -llapacke -Ofast -lblas -o BINARY && ./BINARY'.\newline

For the python3 program you just need to write 'python3 plotar.py'.\newline

