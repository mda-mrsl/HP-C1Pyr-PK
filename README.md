# HP-C1Pyr-PK
Matlab code for pharmacokinetic modeling of hyperpolarized [1-13C]-pyruvate.

These scripts calculate the relationship between [1-13C]-pyruvate and [1-13C]-lactate, but can easily be extended to include other metabolites.

There are three base models, of differing complexity, which include a different number of physical compartments:

The scripts titled P2L3*.m describe signal evolution between two chemical pools and three physical compartments.  The three-compartment model is the most physiologically accurate among these models, but it is also the most computationally intensive, requiring the largest number of descriptive parameters.

The scripts titled P2L2*.m are for the PK model with two chemical pools and two physical compartments.  Here, HP pyruvate and lactate in the extravascular/extracellular space is assumed to be very well mixed with HP pyruvate in the intracellular space - and separate from pyr/lac in the vascular supply.

The scripts titled P2L1*.m represent a simple precursor-product model, which does not consider physical compartmentalization of imaging agents.  Equivalently, this assumes that all compartments (vascular; extravascular/extracellular; and cellular) are well mixed.

Each of these groups of scripts are also accompanied by a test script (TestP2L#.m) that shows how to work with these models - from setup and generation of synthetic data to PK analysis by fitting.

Please do not hesitate to contact with any questions, comments, concerns, or suggestions for improvement:
Jim Bankson
jbankson@mdanderson.org


Python scripts are also available as a package on pypi repository. Install via 
pip install pkhpc1pyr

Python package code available at github.com/RytyB/HP-C1Pyr-PK-python

For questions concerning the python implementation please contact
Ryan Boyce 
rboyce@mdanderson.org
