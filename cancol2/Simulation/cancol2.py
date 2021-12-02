
from cc3d import CompuCellSetup
        

from cancol2Steppables import cancol2Steppable
CompuCellSetup.register_steppable(steppable=cancol2Steppable(frequency=1))



# Not sure about next 2 lines
pyAttributeAdder,dictAdder=CompuCellSetup.attachDictionaryToCells(sim)

OrientedGrowthPlugin = CompuCell.getOrientedGrowthPlugin()
# ----------



from cancol2Steppables import CellLayoutSteppable
CompuCellSetup.register_steppable(steppable=CellLayoutSteppable(frequency=1))

from cancol2Steppables import VolumeParamSteppable
CompuCellSetup.register_steppable(steppable=VolumeParamSteppable(frequency=1))

from cancol2Steppables import MatrixDegradation
CompuCellSetup.register_steppable(steppable=MatrixDegradation(frequency=1))

from cancol2Steppables import MitosisSteppable
CompuCellSetup.register_steppable(steppable=MitosisSteppable(frequency=1))

from cancol2Steppables import CellMotilitySteppable
CompuCellSetup.register_steppable(steppable=CellMotilitySteppable(frequency=1))

from cancol2Steppables import SecretionSteppable
CompuCellSetup.register_steppable(steppable=SecretionSteppable(frequency=1))

from cancol2Steppables import OrientedConstraintSteppable
CompuCellSetup.register_steppable(steppable=OrientedConstraintSteppable(frequency=4))






CompuCellSetup.run()
