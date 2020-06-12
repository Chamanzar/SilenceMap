# SilenceMap
SilenceMap is a novel algorithm for detection of neural “silences” using noninvasive scalp electroencephalography (EEG) signals. Regions of “silences” are defined as parts of the brain tissue with little or no neural activity, e.g., ischemic, necrotic, or lesional tissue in stroke, traumatic brain injuries (TBIs), intracranial hematoma, or even tumors in the brain. SilenceMap uses a novel hemispheric baseline approach, and with the aid of a convex spectral clustering (CSpeC) framework, provides fast detection and localization of the regions of silence in the brain based on a relatively small amount of scalp EEG data. SilenceMap was introduced in our recent paper[1].

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

We used the FreeSurfer (v6.0.0) and MNE (v0.14.0) open-source software for processing the MRI scans in this study. We preprocessed the recorded EEG signals using EEGLAB toolbox (v2019.0) in MATLAB. We calculated the lead-field/forward matrix using the FieldTrip MATLAB toolbox (fieldtrip-20170828). The SilenceMap algorithm was developed in MATLAB (R2018b), using standard toolboxes, and the CVX MATLAB package (v2.1). 

### Installing

You need to install the required toolboxes/software, i.e., MATLAB, CVX toolbox, EEGLAB toolbox, and FieldTrip MATLAB toolbox. The anonymized raw EEG dataset and MRI scans of the participants in this research are made available online on KiltHub, Carnegie Mellon University’s online data repository (DOI: 10.1184/R1/12402416) [2].  

## Running the tests

To test the performance of the SilenceMap algorithm, we have included example codes in this project which simulate different regions of silence at random locations in the brain model, and apply the SilenceMap algorithm to localize the simulated region. To make it easy for you to run this example, we have preprocessed the MRI scans and included the extracted headmodels and the corresponding leadfield matrices.  

### EEG preprocessing
[EEG_prep/EEG_epoched_pruning.m]



### MRI preprocessing and leadfield extraction 

### SilenceMap

Two example files are included in this project, which show how to use and test the SilenceMap algorithm: 1) "SilenceMap_with_baseline.m" which tests the SilenceMap algorithm with hemispheric baseline, and 2) "SilenceMap_without_baseline.m" which tests the SilenceMap algorithm without hemispheric baseline. Each of these two codes generates, tests, and saves different simulated regions of silence, at random locations on in the brain model. 



All you need to do is to put the required files in a folder, i.e., "SilenceMap.m", "", and "", and simply run the following lines in MATLAB: 

```
cd /path-to-the-folder/
run SilenceMap_with_baseline.m
run SilenceMap_wo_baseline.m
```


## Built With

* [MATLAB](https://www.mathworks.com/products/matlab.html) - The standard toolboxes are required.
* [CVX](http://cvxr.com/cvx/) - Matlab Software for Disciplined Convex Programming.
* [FieldTrip](http://www.fieldtriptoolbox.org) - FieldTrip is the MATLAB software toolbox for MEG, EEG, iEEG and NIRS analysis.
* [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) - EEGLAB is an interactive Matlab toolbox used here for preprocessing EEG signals.
* [FreeSurfer](https://surfer.nmr.mgh.harvard.edu) - An open source software suite for processing and analyzing (human) brain MRI images.
* [MNE](https://mne.tools/stable/index.html) - Open-source Python software used here for processing the MRI scans and extraction of different layers of the brain.



## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/Chamanzar/SilenceMap/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

## Refrences

[1] A. Chamanzar, M. Behrmann, and P. Grover, "Neural silences can be localized rapidly using noninvasive scalp EEG", to be published, 2020. 

[2] A. Chamanzar, M. Behrmann, and P. Grover, "Pediatric patients with lobectomy (MRI and EEG)", Carnegie Mellon University. Dataset. https://doi.org/10.1184/R1/12402416

