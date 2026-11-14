import importlib.util,sys
from pathlib import Path
sys.path.insert(0,str(Path.home()/'tools/aipharma-prep-python'))
for name in ['openmm','pdbfixer','openff.toolkit','openmmforcefields','parmed','mdtraj']:
 try:print(name,importlib.util.find_spec(name))
 except Exception as e:print(name,str(e))
import openmm
print('platforms',[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])
