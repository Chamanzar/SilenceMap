# SilenceMap
SilenceMap is a novel method for detection of neural “silences” using noninvasive scalp electroencephalography (EEG) signals. Regions of “silences” are defined as parts of the brain tissue with little or no neural activity, e.g., ischemic, necrotic, or lesional tissue in stroke, traumatic brain injuries (TBIs), intracranial hematoma, or even tumors in the brain. SilenceMap uses a novel hemispheric baseline approach, and with the aid of a convex spectral clustering (CSpeC) framework, provides fast detection and localization of the regions of silence in the brain based on a relatively small amount of scalp EEG data. SilenceMap was introduced in our recent paper [1].

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

We used the FreeSurfer (v6.0.0) and MNE (v0.14.0) open-source software for processing the MRI scans in this study. We preprocessed the recorded EEG signals using EEGLAB toolbox (v2019.0) in MATLAB. We calculated the lead-field/forward matrix using the FieldTrip MATLAB toolbox (fieldtrip-20170828). The SilenceMap algorithm was developed in MATLAB (R2018b), using standard toolboxes, and the CVX MATLAB package (v2.1). 

### Installing

You need to install the required toolboxes/software, i.e., MATLAB, CVX toolbox, EEGLAB toolbox, and FieldTrip MATLAB toolbox.   

## Running the tests

Explain how to run the automated tests for this system


### And coding style tests

Explain what these tests test and why

```
Give an example
```


## Built With

* [MATLAB](https://www.mathworks.com/products/matlab.html) - The standard toolboxes are required
* [CVX](http://cvxr.com/cvx/) - Matlab Software for Disciplined Convex Programming
* [FieldTrip](http://www.fieldtriptoolbox.org) - FieldTrip is the MATLAB software toolbox for MEG, EEG, iEEG and NIRS analysis.
* [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) - Used to generate RSS Feeds
* [FreeSurfer](https://rometools.github.io/rome/) - Used to generate RSS Feeds
* [MNE](https://rometools.github.io/rome/) - Used to generate RSS Feeds





## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

