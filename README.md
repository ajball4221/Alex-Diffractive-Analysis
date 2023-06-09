# Alex-Diffractive-Analysis
My cloned CC-NuE-XSec folder that contains all of the code I've modified and built for my diffractive pion production analysis

## Steps to Getting a Cross Section Using my Code:

1. Preweighted Event Selection: CC-NuE-XSec/selection/eventSelection.py --unweighted (This gets you preweighted ROOT histograms that just have a fudge factor to make IUE agree better and you could add your own fudge factors in fitBackground.py if you restructure it a bit. I put some scripts in here to run these.)
2. Combining Event Selection Output: CC-NuE-XSec/combine_file.py (Can only run this off of things made by eventSelection or gridSelection, it just puts together all the ROOT files into one. Run this separately for data and mc/dfr to get a data file and an mc file)
3. Getting Weights from Fitter based on Eel distribution: CC-NuE-XSec/tools/fitBackground.py (not needed but good to look at weights sometimes)
4. Get Preweighted Plots: CC-NuE-XSec/selection/plotSelection.py (Also not needed but good to make sure nothing has gone horribly wrong at this point before you waste hours running event selection again)
5. Move the preweighted ROOT files into the tools directory and have the fitBackground.py point to them (this is clunky but the grid can only access the CC-NuE-XSec directory
6. Final Event Selection: CC-NuE-XSec/selection/eventSelection.py (now don't add the --unweighted flag, there's other bash scripts for these)
7. Do combine_file and plotSelection again to see results
8. Efficiency Correction: CC-NuE-XSec/efficiency/efficiencyCalculator.py (this corrects for the fact that many of the events that happen are cut out in order to get the true cross section eventually)
9. Unfolding: CC-NuE-XSec/unfolding/matrixFitting.py (this converts detector distributions into truth distributions either using an inverse matrix method or an exponential fit, chef's choice)
10. Getting Cross Section: CC-NuE-XSec/xsec/AlexXSection.py (this normalizes the unfolded distributions to get the true cross sections)

## Other Helpful Files You Will Need:
