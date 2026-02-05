## `postprocessing/` Folder Overview

This folder contains the python files the are used for the postprocessing of the simulation data. In particular:

- `BGK_read.py`: Contains the functions to read the simulation outputs `.txt` files and store the data in `numpy` tensors

- `BGK_plot.py`: Contains the functions to plot the contents of the simulation outputs `.txt` files

- `ytrehus_sol.py`: Contains a set of functions that implement the thoretical result on evaporation presented by Ytrehus in his article (check report)

In the [Google Drive Folder](https://drive.google.com/drive/folders/14fjwMq2cQmRHNADEw7G7EoBn30CHxR80?usp=sharing) a virtual enviroment (`Bgk-venv`) can be found that contains all the necessary dependencies. After a personalized plotting file is created with the desired functions run:

```bash
source path/to/venv/Bgk-venv/bin/activate
python -m postprocessing.my_plot
```
