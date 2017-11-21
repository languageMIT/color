### Tsimane' Color Project

Data and standalone analysis pipeline for generating figures and data tables for the color paper:

Edward Gibson, Richard Futrell, Julian Jara-Ettinger, Kyle Mahowald, Leon Bergen, Sivalogeswaran Ratnasingam, Mitchell Gibson, Steven T. Piantadosi, and Bevil Conway. 2017. [Color naming across languages reflects color use](http://www.pnas.org/content/early/2017/09/12/1619666114.full). *Proceedings of the National Academy of Sciences* 114(40): 10785-10790.

The data is anonymized: participant names and locations are replaced with numeric codes.

The analysis is done with a large R script and then a small python script. To install the R dependencies, open the R interpreter and run `install.packages(c("MASS", "lme4", "reshape2", "plyr", "stringr", "hexbin", "Hmisc", "tidyverse"))`. To install the python dependencies, on the command line do `pip install pandas`, `pip install scipy`, `pip install matplotlib`. 

Now to run the analysis, on the command line, do: `Rscript allfigs.R`, then `python contours.py`. Results will be stored in the directory `output`.


