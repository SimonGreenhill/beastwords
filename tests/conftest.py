from pathlib import Path

import pytest

from beastwords.main import Converter, CovarionConverter, CTMCConverter

### --------------------------------------------------------------------------------------------------###
### Fixtures
### --------------------------------------------------------------------------------------------------###
@pytest.fixture
def covarion():
    return Converter.from_file(Path(__file__).parent / 'overall-covarion.xml')


@pytest.fixture
def covarionNMR():  # no mutation rate
    return Converter.from_file(Path(__file__).parent / 'overall-covarion-no_mutationrate.xml')


@pytest.fixture
def ctmc():
    return Converter.from_file(Path(__file__).parent / 'overall-ctmc.xml')


@pytest.fixture
def covarionPartSize2():  # repartitioned into 2 
    o = Converter.from_file(Path(__file__).parent / 'overall-covarion.xml')
    o.set_partitions(2)
    return o


@pytest.fixture
def ctmcPartSize3():  # repartitioned into 3 by partition size method
    o = Converter.from_file(Path(__file__).parent / 'overall-ctmc.xml')
    o.set_partitions(3)
    return o


@pytest.fixture
def covarionGroupSize2():  # repartitioned into 3 by group size method
    o = Converter.from_file(Path(__file__).parent / 'overall-covarion.xml')
    o.set_partitions("1-6,7-9")
    return o

