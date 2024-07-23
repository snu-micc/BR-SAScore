# BR-SAScore
Implementation of BR-SAScore developed by Prof. Yousung Jung group at Seoul National University (contact: yousung@gmail.com).

## Contents
- [Developer](#developer)
- [OS Requirements](#os-requirements)
- [Python Dependencies](#python-dependencies)
- [Installation Guide](#installation-guide)
- [Usage](#usage)
- [Data](#data)
- [Reproduce the results](#reproduce-the-results)
- [Publication](#publication)
- [License](#license)

## Developer
Shuan Chen (shuan.micc@gmail.com)<br>

## OS Requirements
This repository has been tested on both **Linux** and **Windows** operating systems.

## Python Dependencies
* Python (version >= 3.6)
* Numpy (version >= 1.16.4)
* rdkit (version >= 2019)

## Installation Guide
### From Github
```
pip install BRSAScore
```

### From Github
```
git clone https://github.com/snu-micc/BR-SAScore.git
cd BR-SAScore
pip install -e .
```

## Usage
```
from BRSAScore import SAScorer
scorer = SAScorer()
smi = 'CC(OC1=CC=CC=C1C(O)=O)=O' # Aspirin
score, contribution = scorer.calculateScore(smi)
```
The expected output of `score` should be
```
0.8789047205405656
```
See `Demo.ipynb` for the examples of estimating the synthetic accessibility of 18 structurally complex molecules shown in the paper.


## Data
#### USPTO reaction dataset
The reactions of of the training reactions `train_data.csv` are downloaded from the Dropbox link available at the GitHub repo of [GLN](https://github.com/Hanjun-Dai/GLN), originally from [Figshare](https://doi.org/10.6084/m9.figshare.25046471.v1).

#### Building blocks
The zip file of eMolecules building blocks `origin_dict.csv` can be downloaded from the Dropbox link available the GitHub repo of [Retro*](https://github.com/binghong-ml/retro_star).

## Reproduce the results
### [1] Download reaction and building-block data
Downlaod the data can put them in `./data/` and rename them to `uspto.csv` and `emolecules.csv` with column titled `reactants>reagents>production` and `SMILES`.

Feel free to change the (atom-mapped) reaction data and building blocks data as long as the data is correctly formatted.
If the reaction data is not atom-mapped, we recommend to use [LocalMapper](https://github.com/snu-micc/LocalMapper/tree/main) to prepare high-quality the atom-mappings for your reactions.

### [2] Prepare the BRScores from reaction and building-block data
Prepare the BRSAScores by
```
python scripts/prepare_Scores.py -r uspto -b emolecules
```

### [3] Use BRSAScocre
```
from BRSAScore import SAScorer
reaction_from = 'uspto'
buildingblock_from = 'emolecules'
scorer = SAScorer(reaction_from=reaction_from, buildingblock_from=buildingblock_from)
```
To use the BRSAScore on you own data, replace the `reaction_from` and `buildingblock_from`

## Publication
@article{chen2021deep,
  title={Estimating the synthetic accessibility of molecules with building block and reaction-aware SAScore},
  author={Chen, Shuan and Jung, Yousung},
  journal={Journal of Cheminformatics},
  volume={},
  number={},
  pages={},
  year={2024},
  publisher={Springer}
}

## License
This project is covered under the **MIT Liscence**.
