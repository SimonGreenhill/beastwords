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
### Convert Sequence / Ascertainment
### --------------------------------------------------------------------------------------------------###
@pytest.fixture
def apartitions():
    return {
        'partitions': {
            'eye': [1, 2, 3, 4],
            'foot': [5, 6, 7, 8, 9],
            'hand': [10, 11, 12],
        },
        'ascertainment': [1, 5, 10],
    }


def test_convert_sequences(covarion, apartitions):
    covarion._convert_sequences()
    assert len(covarion.partitions) == 3, 'number of partitions has changed!'
    # check each partition is correct (note alpha reordering)
    assert covarion.partitions['eye'] == apartitions['partitions']['eye']
    assert covarion.partitions['foot'] == apartitions['partitions']['foot']
    assert covarion.partitions['hand'] == apartitions['partitions']['hand']
    assert covarion.ascertainment == apartitions['ascertainment']


def test_convert_sequences_xmlsequences(covarion, apartitions):
    covarion._convert_sequences()
    # get seqs
    sequences = {s.get('taxon'): [_ for _ in s.get('value') if _.strip()] for s in covarion.root.xpath('.//sequence')}
    print(sequences)
    for values in sequences.values():
        assert len(values) == 12
    
    expected = {   #Xeye Xfoot Xhand
        'Taxon1':  '0?11 01000 011',
        'Taxon2':  '0??1 00100 011',
        'Taxon3':  '???? 00010 011',
    }
    for e in expected:
        assert e in sequences
        print(e, "".join(sequences[e]))
        print(e, "".join([_ for _ in expected[e] if _.strip()]))
        print()
        assert sequences[e] == [_ for _ in expected[e] if _.strip()]


def test_convert_sequences_xmluserdatatype(covarion, apartitions):
    covarion._convert_sequences()
    udt = [
        e.get('characterName') for e in covarion.root.xpath('./data/userDataType')[0].getchildren()
    ]
    # make sure we do not have the old _ascertainment character
    assert len(udt) == 12
    for label in udt:
        assert 'ascertain' not in label, 'looks like the ascertainment character is still there'
    
    for p in apartitions['partitions']:
        for site in apartitions['partitions'][p]:
            assert udt[site - 1].split("_")[0] == p


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
        assert has_id(m.tree, 'OneOnX', f'OneOnX:{p}'), f'Missing renamed OneOnX:{p}'


def test_convert_prior_covarion(covarion):
    covarion._convert_prior()
    ba = get_all(covarion.tree, 'prior', 'bcov_alpha_prior.s:combined')
    assert len(ba) == 1, 'bcov_alpha_prior.s: not renamed'

    bs = get_all(covarion.tree, 'prior', 'bcov_s_prior.s:combined')
    assert len(bs) == 1, 'bcov_s_prior.s: not renamed'

    assert ba[0].get('x') == '@bcov_alpha.s:combined'
    assert bs[0].get('x') == '@bcov_s.s:combined'


def test_convert_prior_ctmc(ctmc, partitions):
    ctmc._convert_prior()
    # <prior id="GammaShapePrior.s:foot" name="distribution" x="@gammaShape.s:foot">
    #     <Exponential id="Exponential.1" name="distr">
    #         <mean id="Function$Constant.1" spec="Function$Constant" value="1.0"/>
    #     </Exponential>
    # </prior>
    for p in partitions:
        assert has_id(ctmc.tree, 'prior', f"GammaShapePrior.s:{p}"), "should have element for eye"
        assert has_id(ctmc.tree, 'Exponential', f"Exponential:{p}")
        assert has_id(ctmc.tree, 'mean', f"Function$Constant:{p}")




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
    
    assert opFD[0].getchildren()[0].get('idref') == "frequencies.s:combined", 'frequenciesDelta not updated'


def test_convert_operators_ctmc(ctmc, partitions):
    ctmc._convert_operators()
    for p in partitions:
        op = ctmc.tree.xpath(f".//operator[@id='FrequenciesExchanger.s:{p}']")
        assert op, f'FrequenciesExchanger.s:{p} not renamed'
        assert op[0].getchildren()[0].get('idref') == f"freqParameter.s:{p}", "%s != %s" % (
            op[0].getchildren()[0].get('idref'), f"freqParameter.s:{p}")



### --------------------------------------------------------------------------------------------------###
### Log
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
    # <distribution id="treeLikelihood.foot" spec="TreeLikelihood" branchRateModel="@StrictClock.c:clock" tree="@Tree.t:tree" useAmbiguities="true">
    #     <data id="orgdata.foot" spec="FilteredAlignment" ascertained="true" excludeto="1" filter="-">
    #         <data id="foot" spec="FilteredAlignment" data="@words" filter="4-7"/>
    #         // COV
    #         <userDataType id="TwoStateCovarion.1" spec="beast.base.evolution.datatype.TwoStateCovarion"/>
    #         // CTMC
    #         <userDataType id="Binary.1" spec="beast.base.evolution.datatype.Binary"/>
    #     </data>
    #     // COV
    #     <siteModel id="SiteModel.s:foot" spec="SiteModel" gammaCategoryCount="1" mutationRate="@mutationRate.s:foot" substModel="@covarion">
    #     // CTMC
    #     <siteModel id="SiteModel.s:foot" spec="SiteModel" gammaCategoryCount="4" mutationRate="@mutationRate.s:foot" shape="@gammaShape.s:foot">
    #         <parameter id="gammaShape.s:foot" spec="parameter.RealParameter" estimate="false" name="shape">1.0</parameter>
    #         <parameter id="proportionInvariant.s:foot" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0">0.0</parameter>
    #     </siteModel>
    # </distribution>
    
    # check we don't have extra (i.e. the old one)
    assert len(m.tree.xpath(f".//distribution[starts-with(@id, 'treeLikelihood.')]")) == len(partitions)
    
    for i, p in enumerate(partitions):
        el = m.tree.xpath(f".//distribution[@id='treeLikelihood.{p}']")
        # check we have all the required attribs on <distribution
        assert el, f'missing treeLikelihood for {p}'
        assert el[0].get('spec') == 'TreeLikelihood'
        
        # partition 1 should have the full branchRateModel but the others should just have a branchRateModel
        # as an attribute
        if i == 0:
            assert el[0].get('branchRateModel') is None  # should not have the branchRateModel attribute
            assert len(el[0].xpath(".//branchRateModel")) == 1
        else:
            assert el[0].get('branchRateModel') == '@StrictClock.c:clock'
        
        assert el[0].get('tree') == '@Tree.t:tree'
        
        # data element
        data = el[0].xpath(".//data")
        assert len(data) == 2, f'missing data elements for treeLikelihood.{p}: {data}'
        
        assert p in data[0].get('id')
        assert data[0].get('spec') == 'FilteredAlignment'
        assert data[0].get('ascertained') == 'true'
        assert data[0].get('excludeto') == '1'
        assert data[0].get('filter') == '-'
        
        # nested data element
        assert data[1].getparent() == data[0]
        assert data[1].get('id') == p
        assert data[1].get('spec') == 'FilteredAlignment'
        assert data[1].get('filter') != '-'

        # site model
        siteModel = el[0].xpath('.//siteModel')
        assert siteModel, f'missing siteModel for {p}'
        assert siteModel[0].get('id') == f'SiteModel.s:{p}'
        assert siteModel[0].get('spec') == 'SiteModel'
        assert siteModel[0].get('gammaCategoryCount') == '1'
        assert siteModel[0].get('mutationRate') == f'@mutationRate.s:{p}'
        
        assert siteModel[0].get('substModel')   # Value is tested in covarion/ctmc models
        
        # check siteModel parameters
        gs, pi = siteModel[0].getchildren()
        assert gs is not None
        assert gs.get('id') == f"gammaShape.s:{p}"
        assert gs.get('spec') == "parameter.RealParameter"
        assert gs.get('name') == "shape"

        assert pi is not None
        assert pi.get('id') == f"proportionInvariant.s:{p}"
        assert pi.get('spec') == "parameter.RealParameter"
        assert pi.get('name') == "proportionInvariant"
    

def test_convert_treelikelihood_covarion(covarion, partitions):
    covarion._convert_treelikelihood()
    # check all userDataTypes
    #<userDataType id="TwoStateCovarion.1" spec="beast.base.evolution.datatype.TwoStateCovarion"/>
    udt = covarion.root.xpath(".//distribution[@spec='TreeLikelihood']//userDataType")
    assert len(udt) == len(partitions)
    
    for i, u in enumerate(udt):
        assert u.get('id') == f'userDataType:{partitions[i]}'
        assert u.get('spec') == 'beast.base.evolution.datatype.TwoStateCovarion'

    # check that we have useAmbiguities=True
    for p in partitions:
        d = covarion.root.xpath(f".//distribution[@id='treeLikelihood.{p}']")
        assert d[0].get('useAmbiguities') == 'true'

    for p in partitions:
        sm = covarion.root.xpath(f".//distribution/siteModel[@id='SiteModel.s:{p}']")
        assert len(sm), f'siteModel.s:{p} missing'
        assert sm[0].get('substModel') == '@covarion:combined'


def test_convert_treelikelihood_ctmc(ctmc, partitions):
    ctmc._convert_treelikelihood()
    # check that we have useAmbiguities=False
    for p in partitions:
        d = ctmc.root.xpath(f".//distribution[@id='treeLikelihood.{p}']")
        assert d[0].get('useAmbiguities', 'false') == 'false'
    
    # check userDataType
    udt = ctmc.root.xpath(".//distribution[@spec='TreeLikelihood']//userDataType")
    assert len(udt) == len(partitions)
    for i, u in enumerate(udt):
        assert u.get('id') == f'userDataType:{partitions[i]}'
        assert u.get('spec') == 'beast.base.evolution.datatype.Binary'
    
    # check substModel
    for p in partitions:
        sm = ctmc.root.xpath(f".//distribution/siteModel[@id='SiteModel.s:{p}']")
        assert len(sm), f'siteModel.s:{p} missing'
        assert sm[0].get('substModel') == f'@CTMC.s:{p}'


### --------------------------------------------------------------------------------------------------###
### Misc
### --------------------------------------------------------------------------------------------------###
@pytest.mark.parametrize("fixture", ["covarion", "ctmc"])
def test_clock(request, fixture):
    m = request.getfixturevalue(fixture)
    # we should have a clock before we change things
    clock = m.root.xpath(".//branchRateModel")
    assert len(clock) == 1, 'missing clock??'
    
    m.convert()
    
    clock = m.root.xpath(".//branchRateModel")
    assert len(clock) == 1, 'missing clock'

    