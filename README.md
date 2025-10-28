# fits1

This exercise will begin to introduce non linear data fitting.  The examples will be based on the ROOT libraries.  Similar functionality for everything we see can be found in the numpy + scipy + lmfit + matplotlib modules. But the examples we will follow over the next few projects offer several advantages:
* far more natural histogramming tools
* completely automated fitting options ("one liners")
* convenient methods for defining custom functions and controlling the convergence of the fit
* detailed and consistent documentation
* and a low level interface to modify the objective function, running in fully optimized compiled code

You are welcome to modify the provided code for your projects and to use other packages.  Where applicable alternate examples will be included. 

* **fit1.C**: C++ function to generate random data according to a normal distribution with mean=20, sigma=10. <br> A fit is performed using the built in Gaussian model in ROOT.  Then the parameter values, their uncertainteis, and the p-value for the fit are extracted.  To run this code type ```root fit1.C``` or if you are already running ROOT, type ```.X fit1.C```  
* **fit1.py**: The same code using the python interface, run this example using ```python fit1.C```.
* For a contrast see **fit1mpl.py** for a version using matplotlib+scipy.  
* readhist.C(py):  Examples for reading the histogram files given in this example 
* ParamUnceratinties.ipynb : a guided tutorial towards most of what you will be coding in this week's exercise.
* LLexample.ipynb : a notebook giving an example for calculating (N)LLs
* TH1hist2Numpy.ipynb : an example for converting a ROOT histogram to numpy arrays

Note that from ROOT you can type ```new TBrowser()``` or in Python r.TBrowser() to get a graphical browser that allows you to look at what's contained in the TFiles.

p-value estimated = 0.466400