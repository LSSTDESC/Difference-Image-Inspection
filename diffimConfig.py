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
