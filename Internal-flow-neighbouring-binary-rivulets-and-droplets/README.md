# Internal-flow-neighbouring-binary-rivulets-and-droplets
This repository contains the code used to reproduce the simulations reported in the manuscript “Internal flow and concentration in neighbouring evaporating binary droplets and rivulets” ([link](arXiv:2508.08820)). The paper and dataset explore how proximity between neighbouring binary droplets or rivulets alters evaporation-driven flows and solute distributions.

Citation (currently for preprint):
```
@article{dekker2025internal,
  title={Internal flow and concentration in neighbouring evaporating binary droplets and rivulets},
  author={Dekker, Pim J and Rocha, Duarte and Diddens, Christian and Lohse, Detlef},
  journal={arXiv preprint arXiv:2508.08820},
  year={2025}
}
```

## Prerequisites
Assuming you have installed Python 3.9 to 3.13, you can install pyoomph via:
  ```
  python -m pip install pyoomph
  ```
On Apple Silicon (M1–M4) run the above in a Rosetta 2 terminal. If installation fails, consult the pyoomph installation guide: https://pyoomph.readthedocs.io/en/latest/tutorial/installation.html. For persistent issues, contact the authors.

## Introduction
We study evaporation-driven internal flow and solute concentration in neighbouring sessile binary droplets and rivulets. When two liquid bodies evaporate near each other, the vapour concentration between them rises relative to their exposed outer sides, producing a shielding effect that modifies local evaporation rates. This work quantifies how flow and concentration asymmetries depend on key control parameters: the separation distance, contact angle, and solutal Marangoni forces.

## Guide through scripts 
The repository contains the following scripts:
- `transient_rivulets.py`: Time-dependent simulation of two neighbouring rivulets based on the model in Sec. 2.1; reproduces the results shown in Fig. 1 of the manuscript.
- `quasi_stationary_rivulets.py`: Quasi-stationary simulation of two neighbouring rivulets base on the model in Sec 3.1, reproduces the results shown in Fig. 3 of the manuscript.
- `quasi_stationary_lubrication_droplets.py`: Quasi-stationary simulation of two neighbouring droplets based on the lubrication model in Sec 3.2, reproduces the results shown in Fig. 8 of the manuscript.

## Usage
To run the simulations, navigate to the directory containing the desired script and execute it using Python. For example:
```
python3 transient_rivulets.py
```
This will run the transient simulation for neighbouring rivulets as described in the manuscript.
Results will be saved in the output folder with the same name as the script.