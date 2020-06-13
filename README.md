# SilenceMap
SilenceMap is a novel algorithm for detection of neural “silences” using noninvasive scalp electroencephalography (EEG) signals. Regions of “silences” are defined as parts of the brain tissue with little or no neural activity, e.g., ischemic, necrotic, or lesional tissue in stroke, traumatic brain injuries (TBIs), intracranial hematoma, or even tumors in the brain. SilenceMap uses a novel hemispheric baseline approach, and with the aid of a convex spectral clustering (CSpeC) framework, provides fast detection and localization of the regions of silence in the brain based on a relatively small amount of scalp EEG data. SilenceMap was introduced in our recent paper[1].

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

We used the FreeSurfer (v6.0.0) and MNE (v0.14.0) open-source software for processing the MRI scans in this study. We preprocessed the recorded EEG signals using EEGLAB toolbox (v2019.0) in MATLAB. We calculated the lead-field/forward matrix using the FieldTrip MATLAB toolbox (fieldtrip-20170828). The SilenceMap algorithm was developed in MATLAB (R2018b), using standard toolboxes, and the CVX MATLAB package (v2.1). 

### Installing

You need to install the required toolboxes/software, i.e., MATLAB, CVX toolbox, EEGLAB toolbox, FreeSurfer, MNE, and FieldTrip MATLAB toolbox. The anonymized raw EEG dataset and MRI scans of the participants in this research are made available online on KiltHub, Carnegie Mellon University’s online data repository (DOI: 10.1184/R1/12402416) [2].  

## Running the tests

To test the performance of the SilenceMap algorithm, we have included example codes in this project which simulate different regions of silence at random locations in the brain model, and apply the SilenceMap algorithm to localize the simulated region. To make it easy for you to run this example, we have preprocessed the MRI scans and included the extracted headmodels and the corresponding leadfield matrices.  

### EEG preprocessing
In [EEG_prep](EEG_prep) folder, there is a matlab code named [EEG_epoched_pruning.m](EEG_prep/EEG_epoched_pruning.m), which includes all of the preprocessing steps we have used for the recorded EEG dataset. This code requires [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) toolbox. You can download any of the EEG files from our opensource dataset [here](https://doi.org/10.1184/R1/12402416), and put it inside the [EEG_prep](EEG_prep) folder, and run [EEG_epoched_pruning.m](EEG_prep/EEG_epoched_pruning.m), with proper modifications, to preprocess the EEG data.    

### MRI preprocessing and leadfield extraction 

To Prerprocess the MRI scans using FreeSurfer and MNE, please go to [MRI_prep_leadfield_ext](MRI_prep_leadfield_ext) folder and follow the instructions in [MRI_workflow.txt](MRI_workflow). To extract the leadfiled matrix (forward model), please follow the instructions in [Leadfield_workflow.txt](Leadfield_workflow.txt). As an example, we have included the extracted headmodels and their corresponding leadfield matrices for the participant OT (see [here](https://doi.org/10.1184/R1/12402416) for more information), for a low and a high resolution source grid. MRI scans of participants in our study can be downloaded from [here](https://doi.org/10.1184/R1/12402416).      

### SilenceMap

Two example files are included in [SilenceMap](SilenceMap), which show how to use and test the SilenceMap algorithm: 1) "[SilenceMap_with_baseline.m](SilenceMap_with_baseline.m)" which tests the SilenceMap algorithm with hemispheric baseline, and 2) "[SilenceMap_without_baseline.m](SilenceMap_without_baseline.m)" which tests the SilenceMap algorithm without hemispheric baseline. Each of these two codes generates, tests, and saves different simulated regions of silence, at random locations in the brain model. In addtion, we have implemented the state-of-the-art source localization algorithms, i.e., sLORETA, MNE, and MUSIC, with proper modifications for the silence localization task (see [1] for more information). These implementations are available as a matlab function in [modified_src_loc.m](modified_src_loc.m), and an example of how to use this function is available in [sLoreta_MNE_MUSIC_simulation_vF.m](sLoreta_MNE_MUSIC_simulation_vF.m). [hemispheric_base.m](hemispheric_base.m) in [SilenceMap](SilenceMap) calculates the source contribution measure with hemispheric baseline, [CSpeC.m](CSpeC.m) includes the implementation of the convex spectral clustering (CSpeC) framework used in the SilenceMap algorithm, and [plot_source_space_signal_vF.m](plot_source_space_signal_vF.m) is a function to plot the localized regions of silence in the brain model.   


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

First/Corresponding Author Contact Information  

    Name: Alireza Chamanzar  
    Institution: Carnegie Mellon University  
    Address: Hamerschlag Hall B200, 5000 Forbes Ave., Pittsburgh PA 15213 United States  
    Email: achamanz@andrew.cmu.edu  
    
Second Author Contact Information  

    Name: Marlene Behrmann 
    Institution: Carnegie Mellon University
    Address: 331 Baker Hall, 5000 Forbes Ave., Pittsburgh PA 15213 United States
    Email: behrmann@cmu.edu

Third/Corresponding Author Contact Information  

    Name: Pulkit Grover
    Institution: Carnegie Mellon University
    Address: Hamerschlag Hall B202, 5000 Forbes Ave., Pittsburgh PA 15213 United States
    Email: pgrover@andrew.cmu.edu

## License

This project is licensed - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

This work was supported in part by the Chuck Noll Foundation for Brain Injury Research, and a CMU BrainHub award. The authors of this manuscript thank Chaitanya Goswami, Praveen Venkatesh, Ashwati Krishnan, Anne Margarette S. Maallo, Sarah M. Haigh, and Maysamreza Chamanzar for helpful discussions, and Ashwati Krishnan and Patricia Brosseau for help in data collection.

## Refrences

[1] A. Chamanzar, M. Behrmann, and P. Grover, "Neural silences can be localized rapidly using noninvasive scalp EEG", to be published, 2020. 

[2] A. Chamanzar, M. Behrmann, and P. Grover, "Pediatric patients with lobectomy (MRI and EEG)", Carnegie Mellon University. Dataset. https://doi.org/10.1184/R1/12402416

