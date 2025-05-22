from pathlib import Path
import pytest

from lxml import etree

from beastwords.main import Converter, CovarionConverter, CTMCConverter

### --------------------------------------------------------------------------------------------------###
### Fixtures
### --------------------------------------------------------------------------------------------------###
@pytest.fixture
def partitions():
    return ('hand', 'foot', 'eye')


@pytest.fixture
def covarion():
    return Converter.from_file(Path(__file__).parent / 'overall-covarion.xml')


@pytest.fixture
def ctmc():
    return Converter.from_file(Path(__file__).parent / 'overall-ctmc.xml')



### --------------------------------------------------------------------------------------------------###
### Test Helpers
### --------------------------------------------------------------------------------------------------###
def has_id(tree, tag, id_value):
    """
    Returns True if the XML tree contains an element with the given tag and id.
    """
    return bool(tree.xpath(f".//{tag}[@id='{id_value}']"))


def get_all(tree, tag, id_value):
    return tree.xpath(f".//{tag}[starts-with(@id, '{id_value}')]")



### --------------------------------------------------------------------------------------------------###
### Methods etc
### --------------------------------------------------------------------------------------------------###

def test_model(covarion, ctmc):
    assert isinstance(covarion, CovarionConverter)
    assert isinstance(ctmc, CTMCConverter)

    assert covarion.model == 'BinaryCovarion', 'Got %s' % covarion.model
    assert ctmc.model == 'BinaryCTMC', 'Got %s' % ctmc.model


def test_patch(covarion):
    root = etree.Element("root")
    child = etree.SubElement(root, "child")
    child.set('id', 'a')
    child.set('x', 'b')
    tree = etree.ElementTree(root)
    
    new = covarion.patch(child, {'id': 'z', 'test': 'true'})
    assert new.get('id') == 'z'
    assert new.get('test') == 'true'
    assert new.get('x') == 'b'
    # old is not changed
    assert child.get('id') == 'a'


def test_patch_child_ids(covarion):
    prior = etree.Element("prior", id="GammaShapePrior.s:overall", name="distribution", x="@gammaShape.s:overall")
    exp = etree.SubElement(prior, "Exponential", id="Exponential.0", name="distr")
    mean = etree.SubElement(exp, "mean", id="Function$Constant.0", spec="Function$Constant", value="1.0")

    #  <prior id="GammaShapePrior.s:overall" name="distribution" x="@gammaShape.s:overall">
    #    <Exponential id="Exponential.0" name="distr">
    #      <mean id="Function$Constant.0" spec="Function$Constant" value="1.0"/>
    #    </Exponential>
    #  </prior>
    
    new = covarion.patch_child_ids(prior, 'test')
    assert new.get('id') == 'GammaShapePrior.s:test'
    
    exp = new.getchildren()[0]
    assert exp.get('id') == 'Exponential.0:test'

    mean = exp.getchildren()[0]
    assert mean.get('id') == 'Function$Constant.0:test'


def test_replace(covarion):
    xpath = ".//log[starts-with(@idref, 'mutationRate')]"
    assert len(covarion.root.xpath(xpath)) == 1
    covarion.replace(xpath, idref="mutationRate.s:{}")
    print(covarion)
    assert len(covarion.root.xpath(xpath)) == 3
    assert sorted([f"mutationRate.s:{p}" for p in covarion.partitions]) == \
        sorted([p.get('idref') for p in covarion.root.xpath(xpath)])
    
    # <sequence id="seq_Taxon1" spec="Sequence" taxon="Taxon1" totalcount="2" value="0111000?11"/>
    xpath = ".//sequence[starts-with(@id, 'seq_Taxon1')]"
    covarion.replace(xpath, partition="word.{}")
    assert len(covarion.root.xpath(xpath)) == 3
    assert sorted([f"word.{p}" for p in covarion.partitions]) == \
        sorted([p.get('partition') for p in covarion.root.xpath(xpath)])


@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_words(request, fixture):
    m = request.getfixturevalue(fixture)
    words = m.get_words()
    assert words[0] == ('_ascertainment_0', 'UserDataType.0')
    assert words[9] == ('eye_3', 'UserDataType.9')


def test_get_partition_range(covarion):
    # {'hand': [1, 2], 'foot': [3, 4, 5, 6], 'eye': [7, 8, 9]})
    covarion.partitions = {
        'nothing': [],
        'one_to_four': [1, 2, 3, 4],
        'five': [5],
        'split': [6, 7, 9, 10],  # no 8
        'complex': [11, 12, 15, 16, 18, 20, 21, 22, 23, 24, 25],  # 11-12,15-16,18,20-25
        'sorting': [30, 29, 27, 28, 26],
    }
    assert covarion.get_partition_range("unknown") == ""
    assert covarion.get_partition_range("nothing") == ""
    assert covarion.get_partition_range('one_to_four') == '1-4'
    assert covarion.get_partition_range('five') == '5'
    assert covarion.get_partition_range('split') == '6-7,9-10'
    assert covarion.get_partition_range('complex') == '11-12,15-16,18,20-25'
    assert covarion.get_partition_range('sorting') == '26-30'


@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_parse_word(request, fixture):
    m = request.getfixturevalue(fixture)
    assert m.parse_word("_ascertainment_0") == ["_ascertainment", "0"]
    assert m.parse_word("hand_1") == ['hand', '1']
    assert m.parse_word("hand_889") == ['hand', '889']
    assert m.parse_word("leg_foot_42") == ['leg_foot', '42']
    assert m.parse_word("hand_u_136013") == ['hand', '136013']


@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_partitions(request, fixture):
    partitions, ascertainment = request.getfixturevalue(fixture).get_partitions()
    assert partitions['hand'] == [1, 2]
    assert partitions['foot'] == [3, 4, 5, 6]
    assert partitions['eye'] == [7, 8, 9]
    assert len(partitions) == 3
    assert ascertainment == [0]

### --------------------------------------------------------------------------------------------------###
### State
### --------------------------------------------------------------------------------------------------###
@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_convert_state(request, fixture, partitions):
    m = request.getfixturevalue(fixture)
    # check we have the original ones
    assert has_id(m.tree, 'parameter', 'mutationRate.s:overall')
    assert len(get_all(m.tree, 'parameter', 'mutationRate.s')) == 1

    m._convert_state()

    assert len(get_all(m.tree, 'parameter', 'mutationRate.s')) == 3, "should have 3"
    assert not has_id(m.tree, 'parameter', 'mutationRate.s:overall'), "should have removed old one"
    for p in partitions:
        assert has_id(m.tree, 'parameter', f'mutationRate.s:{p}'), f"should have element for {p}"

    assert not has_id(m.tree, 'parameter', 'mutationRate.s:_ascertainment'), "should not have _ascertainment partition"
    
    
def test_convert_state_covarion(covarion):
    covarion._convert_state()
    assert has_id(covarion.tree, 'parameter', 'bcov_alpha.s:combined'), "should have element for combined"
    assert has_id(covarion.tree, 'parameter', 'bcov_s.s:combined'), "should have element for combine"
    assert has_id(covarion.tree, 'parameter', 'frequencies.s:combined'), "should have element for combined"


def test_convert_state_ctmc(ctmc, partitions):
    ctmc._convert_state()
    
    for p in partitions:
        assert has_id(
            ctmc.tree, 'parameter', f'gammaShape.s:{p}'), f"should have gammaShape element for {p}"
        assert has_id(
            ctmc.tree, 'parameter', f'freqParameter.s:{p}'), f"should have freqParameter element for {p}"
    
    assert not has_id(ctmc.tree, 'parameter', 'gammaShape.s:overall'), 'gammaShape.s:overall should be removed'
    assert not has_id(ctmc.tree, 'parameter', 'freqParameter.s:overall'), 'freqParameter.s:overall should be removed'


### --------------------------------------------------------------------------------------------------###
### Prior
### --------------------------------------------------------------------------------------------------###
@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_convert_prior(request, fixture, partitions):
    m = request.getfixturevalue(fixture)
    # check we have the original ones
    assert has_id(m.tree, 'prior', 'MutationRatePrior.s:overall')
    assert len(get_all(m.tree, 'prior', 'MutationRatePrior.s')) == 1
    
    m._convert_prior()
    assert len(get_all(m.tree, 'prior', 'MutationRatePrior.s')) == 3, "should have 3"
    assert not has_id(m.tree, 'prior', 'MutationRatePrior.s:overall'), "should have removed old one"
    
    for p in partitions:
        assert has_id(m.tree, 'prior', f'MutationRatePrior.s:{p}'), f"should have element for {p}"
    assert not has_id(m.tree, 'parameter', 'MutationRatePrior.s:_ascertainment'), "should not have _ascertainment partition"
    
    # have we updated the children ids?
    for p in partitions:
        # <OneOnX id="OneOnX.0:foot" name="distr"/>
        assert has_id(m.tree, 'OneOnX', f'OneOnX.0:{p}')


def test_convert_prior_covarion(covarion):
    # covarion._convert_prior()
    pass # nothing here yet


def test_convert_prior_ctmc(ctmc, partitions):
    ctmc._convert_prior()
    # <prior id="GammaShapePrior.s:foot" name="distribution" x="@gammaShape.s:foot">
    #     <Exponential id="Exponential.1" name="distr">
    #         <mean id="Function$Constant.1" spec="Function$Constant" value="1.0"/>
    #     </Exponential>
    # </prior>
    for p in partitions:
        assert has_id(ctmc.tree, 'prior', f"GammaShapePrior.s:{p}"), "should have element for eye"
        assert has_id(ctmc.tree, 'Exponential', f"Exponential.0:{p}")
        assert has_id(ctmc.tree, 'mean', f"Function$Constant.0:{p}")




### --------------------------------------------------------------------------------------------------###
### substModel
### --------------------------------------------------------------------------------------------------###
def test_convert_substmodel_covarion(covarion):
    sm = list(covarion._convert_substmodel())
    assert len(sm) == 1
    assert sm[0].get('id') == 'covarion:combined'
    assert sm[0].get('spec') == 'BinaryCovarion'
    assert sm[0].get('alpha') == '@bcov_alpha.s:combined'
    assert sm[0].get('switchRate') == '@bcov_s.s:combined'
    assert sm[0].get('vfrequencies') == '@frequencies.s:combined'

    assert has_id(sm[0], 'parameter', f"hiddenfrequencies.s:combined")
    assert has_id(sm[0], 'frequencies', f"dummyfrequencies.s:combined")


def test_convert_substmodel_ctmc(ctmc, partitions):
    sm = list(ctmc._convert_substmodel())
    assert len(sm) == 3
    
    # squash these into one element
    xml = etree.Element("root")
    for s in sm:
        xml.append(s)
    
    # we should have one sm per partition
    for p in partitions:
        t = xml.xpath(f".//substModel[@id='CTMC.s:{p}']")[0]
        
        assert t.get('id') == f'CTMC.s:{p}'
        assert t.get('spec') == 'GeneralSubstitutionModel'
    
        assert len(t.getchildren()) == 2
        # check subelements
        param = t.xpath('.//parameter')[0]
        assert param.get('id') == f"rates.s:{p}"
        assert param.get('spec') == "parameter.RealParameter"
    
        freqs = t.xpath('.//frequencies')[0]
        assert freqs.get('id') == f"estimatedFreqs.s:{p}"
        assert freqs.get('spec') == "Frequencies"
        assert freqs.get('frequencies') == f"@freqParameter.s:{p}"
    


### --------------------------------------------------------------------------------------------------###
### Operators
### --------------------------------------------------------------------------------------------------###
@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_convert_operators(request, fixture, partitions):
    m = request.getfixturevalue(fixture)
    m._convert_operators()
    # should have one for each partition
    for p in partitions:
        op = m.tree.xpath(f".//operator[@id='mutationRateScaler.s:{p}']")[0]
        assert op.get('spec') == 'ScaleOperator'
        assert op.get('parameter') == f"@mutationRate.s:{p}"
        assert op.get('scaleFactor') is not None
        assert op.get('weight') is not None


def test_convert_operators_covarion(covarion):
    covarion._convert_operators()

    opAlpha = covarion.tree.xpath(".//operator[@id='bcovAlphaScaler.s:combined']")
    opSwitch = covarion.tree.xpath(".//operator[@id='bcovSwitchParamScaler.s:combined']")
    opFD = covarion.tree.xpath(".//operator[@id='frequenciesDelta.s:combined']")
    
    assert opAlpha, 'bcovAlphaScaler not renamed'
    assert opSwitch, 'bcovSwitchParamScaler not renamed'
    assert opFD, 'frequenciesDelta not renamed'
    
    assert opAlpha[0].get('parameter') == "@bcov_alpha.s:combined", 'parameter -> @bcov_alpha.s:combined?'
    assert opSwitch[0].get('parameter') == "@bcov_s.s:combined", 'parameter -> @bcov_s.s:combined?'
    
    assert opFD[0].getchildren()[0].get('idref') == "@frequencies.s:combined", 'frequenciesDelta not updated'


def test_convert_operators_ctmc(ctmc, partitions):
    ctmc._convert_operators()
    for p in partitions:
        op = ctmc.tree.xpath(f".//operator[@id='FrequenciesExchanger.s:{p}']")
        assert op, f'FrequenciesExchanger.s:{p} not renamed'
        assert op[0].getchildren()[0].get('idref') == f"freqParameter.s:{p}", "%s != %s" % (
            op[0].getchildren()[0].get('idref'), f"freqParameter.s:{p}")



### --------------------------------------------------------------------------------------------------###
### Tree Likelihood
### --------------------------------------------------------------------------------------------------###
@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_convert_log(request, fixture, partitions):
    m = request.getfixturevalue(fixture)
    m._convert_log()
    for p in partitions:
        assert m.tree.xpath(f".//log[@idref='treeLikelihood.{p}']")
        assert m.tree.xpath(f".//log[@idref='mutationRate.s:{p}']")


def test_convert_log_covarion(covarion, partitions):
    covarion._convert_log()
    # just some renames
    assert covarion.tree.xpath(".//log[@idref='bcov_alpha.s:combined']")
    assert covarion.tree.xpath(".//log[@idref='bcov_s.s:combined']")
    assert covarion.tree.xpath(".//log[@idref='frequencies.s:combined']")
    

def test_convert_log_ctmc(ctmc, partitions):
    ctmc._convert_log()
    for p in partitions:
        assert ctmc.tree.xpath(f".//log[@idref='freqParameter.s:{p}']")
    
    

### --------------------------------------------------------------------------------------------------###
### Tree Likelihood
### --------------------------------------------------------------------------------------------------###
@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_convert_treelikelihood(request, fixture, partitions):
    m = request.getfixturevalue(fixture)
    m._convert_treelikelihood()
    for p in partitions:
        pass
    
    # check how many substmodels -> cov/ctmc diff
    
    # check the old substmodel was removed.
    #print(m)
    #assert False

#def test_convert_treelikelihood_covarion(covarion):
#    pass # nothing here yet
#
#
#def test_convert_treelikelihood_ctmc(ctmc):
#    pass