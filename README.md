<img src="https://drive.google.com/uc?export=view&id=1VjogeXABRBYJygsGfP8bTxbuLgi1Yevl" width = "890" height = "375"> <img src="https://drive.google.com/uc?export=view&id=1QspeLL4wVjzzIGHuyYVMfy2-77yHECNT" width = "350" height="150">

# Unbinned-Likelihood-Analysis

This project aims to conduct an unbinned log-likelihood analysis to calculate the 30-day signal significance of the Hartlepool nuclear reactor complex at the proposed Boulby mine antineutrino detector. Moreover, this powerful technique provides a base for a new analysis chain that expresses powerful background discrimination, which can be implemented to filter events from an active detector. For more information on the theory behind the analysis, please see my master's thesis, which can be found [here](https://drive.google.com/file/d/1bwOUjMAag0bPoYmGG7djkz9zt5Y5vxPh/view?usp=sharing).

## Prerequisities

Simulations of events at the Boulby detector are output as root files for analysis in the object-based data analysis framework, ROOT, developed by particle physicists [1](https://doi.org/10.1016/S0168-9002(97)00048-X). To install ROOT, please follow the instructions on (<https://root.cern/install/>).

## Download

To download this analysis chain just clone this repository and then move to that directory to execute.

```bash
git clone https://github.com/JamArm99/Unbinned-Likelihood-Analysis.git
cd Unbinned-Likelihood-Analysis
```

## Excecute

To execute **_unbin_like.C_**, your simulation root files must be added to the same directory for it to find them. Once complete, open the header file **_unbin_like.h_** to change the file names inside the data_files vector. The analysis is sensitive to the order of the components; the correct order should be: Hartlepool 1 (Big), Hartlepool 2 (Small), combined singles, world reactors, geoneutrinos, nitrogen - 17, and lithium - 9.

Ensure the ROOT shell is sourced.

```bash
source /path-to-root/thisroot.sh
```

Then the analysis can be excecuted either from the terminal

```bash
root unbin_like.C
```

or within ROOT itself

```bash
root
root[0] .X unbin_like.C()
```

## Outputs