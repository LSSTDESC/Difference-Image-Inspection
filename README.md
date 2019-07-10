# Difference-Image-Inspection
An inspection and comparison of DC2 image differencing algorithms

This repository is currently a dumping ground as we collect ongoing work from various places. Better organization to come.

## Project Goals

1. Compare image subtraction techniques and understand in what ways each algorithm fails or misperforms on DC2 data.
2. Understand these failures on both a qualitative and mathematical level.
3. Draw informed conclusions on challenges we will face on DC2 and what can be done about them.
4. Possibly suggest improvements to the image subtraction pipeline that will be implimented for LSST

## References / Resources:

#### Papers
- [Alard & Lupton 1998](https://ui.adsabs.harvard.edu/abs/1998ApJ...503..325A/abstract): A Method for Optimal Image Subtraction
- [Zackay, et al. (2016) ](https://ui.adsabs.harvard.edu/abs/2016ApJ...830...27Z/abstract): Proper Image Subtraction—Optimal Transient Detection, Photometry, and Hypothesis Testing

#### Internal Documentation
- [DC2 Data Products Overview](https://confluence.slac.stanford.edu/display/LSSTDESC/DC2+Data+Product+Overview)
- [DMTN-021](https://dmtn-021.lsst.io): Implementation of Image Difference Decorrelation
- [DMTN-061](https://dmtn-061.lsst.io): State of image subtraction in the LSST stack

#### Repos
- [Difference Image Analysis Pipeline](https://github.com/LSSTDESC/dia_pipe)
- [DC2 Analysis Tutorials](https://github.com/LSSTDESC/DC2-analysis)

## Running Image Subtraction at NERSC

DMTN 61 outlines how to run image subtraction at NERSC, but is a little out of date. As an alternative example, the following config file can be used:

```python
import os
from lsst.utils import getPackageDir
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask
from lsst.ip.diffim.getTemplate import GetCalexpAsTemplateTask

config.imageDifference.refObjLoader.retarget(LoadIndexedReferenceObjectsTask)
config.imageDifference.refObjLoader.load(os.path.join(getPackageDir('obs_lsstCam'), 'config', 'filterMap.py'))
config.imageDifference.kernelSourcesFromRef = True
config.ccdKey = 'detector'

# config.makeDiffim.doWriteSubtractedExp=True
# config.makeDiffim.doWriteMatchedExp=True
# config.makeDiffim.doDecorrelation=True
config.imageDifference.subtract='al'
config.imageDifference.subtract['zogy'].zogyConfig.inImageSpace=False
# config.getTemplate.retarget(GetCalexpAsTemplateTask)
```

Image subtraction can then be run for a specific visit as follows:

```bash
source /global/cscratch1/sd/rearmstr/example_diffim/setup.sh

imageDifferenceDriver.py /global/cscratch1/sd/rearmstr/example_diffim/Run1.2_data/rerun/coadd-v4 \
    --output test_imdiff --id visit=431306 detector=28 \
    -C diffimConfig.py --config imageDifference.subtract='zogy' \
    --cores 4
```

where the `--config imageDifference.subtract='zogy'` argument can be dropped to do an Alard & Lupton subtraction.
