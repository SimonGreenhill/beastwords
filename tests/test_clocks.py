from pathlib import Path

import pytest
from lxml import etree

from beastwords.clocks import detect_clock, strip_clock, extract_mean, convert_clock


STRICT_XML = """
<beast>
  <run id="mcmc">
    <state id="state">
      <parameter id="clockRate.c:clock" name="stateNode" lower="0.0" upper="1.0">1.0</parameter>
    </state>
    <distribution id="posterior">
      <distribution id="prior">
        <prior id="ClockPrior.c:clock" name="distribution" x="@clockRate.c:clock">
          <LogNormal id="LogNormal.clock" name="distr" M="-9.21" S="1.25"/>
        </prior>
      </distribution>
      <distribution id="likelihood">
        <distribution id="treeLikelihood.x" spec="TreeLikelihood" tree="@Tree.t:tree">
          <branchRateModel id="StrictClock.c:clock"
            spec="beast.base.evolution.branchratemodel.StrictClockModel"
            clock.rate="@clockRate.c:clock"/>
        </distribution>
      </distribution>
    </distribution>
    <operator id="StrictClockRateScaler.c:clock" spec="ScaleOperator" parameter="@clockRate.c:clock" weight="5.0"/>
    <logger id="tracelog"><log idref="clockRate.c:clock"/></logger>
  </run>
</beast>
"""

ORC_XML = """
<beast>
  <run id="mcmc">
    <state id="state">
      <parameter id="ucldMean.c:clock" name="stateNode" lower="0.0" upper="1.0">1.0</parameter>
      <parameter id="ORCsigma.c:clock" name="stateNode">0.1</parameter>
      <parameter id="ORCRates.c:clock" dimension="6" lower="1.0E-100" name="stateNode">1.0</parameter>
    </state>
    <distribution id="posterior">
      <distribution id="prior">
        <prior id="ucldMeanPrior.c:clock" name="distribution" x="@ucldMean.c:clock">
          <LogNormal id="LogNormal.ucldMean" name="distr" M="-9.21" S="1.25"/>
        </prior>
      </distribution>
      <distribution id="likelihood">
        <distribution id="treeLikelihood.x" spec="TreeLikelihood" tree="@Tree.t:tree">
          <branchRateModel id="OptimisedRelaxedClock.c:clock"
            spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel"
            clock.rate="@ucldMean.c:clock" rates="@ORCRates.c:clock" tree="@Tree.t:tree">
            <LogNormal id="LogNormalDistributionModel.c:clock" S="@ORCsigma.c:clock" meanInRealSpace="true" name="distr">
              <parameter name="M">1.0</parameter>
            </LogNormal>
          </branchRateModel>
        </distribution>
      </distribution>
    </distribution>
  </run>
</beast>
"""


from beastwords.main import Converter

FIXTURE = Path(__file__).parent / "overall-ctmc.xml"  # 3 taxa, fixed strict clock


def test_converter_applies_clock_before_partitioning():
    o = Converter.from_file(FIXTURE)
    o.clock = "orc"
    o.convert()
    assert detect_clock(o.root) == "orc"
    # the ORC branchRateModel survives into the partitioned likelihood
    assert o.root.xpath(".//branchRateModel[@id='OptimisedRelaxedClock.c:clock']")
    ids = o.root.xpath(".//*/@id")
    assert sorted({i for i in ids if ids.count(i) > 1}) == []


def test_converter_default_preserves_clock():
    o = Converter.from_file(FIXTURE)
    o.convert()
    assert detect_clock(o.root) == "strict"


def test_cli_clock_orc(tmp_path, monkeypatch):
    import sys
    from beastwords.main import main
    out = tmp_path / "out.xml"
    monkeypatch.setattr(sys, "argv",
        ["beastwords", "--clock", "orc", str(FIXTURE), str(out)])
    main()
    root = etree.parse(str(out)).getroot()
    assert detect_clock(root) == "orc"


def _load():
    return etree.parse(str(FIXTURE)).getroot()


def _dup_ids(root):
    ids = root.xpath(".//*/@id")
    return sorted({i for i in ids if ids.count(i) > 1})


def test_convert_strict_to_orc():
    root = _load()
    convert_clock(root, "orc")
    assert detect_clock(root) == "orc"

    brm = root.xpath(".//branchRateModel")
    assert len(brm) == 1
    assert brm[0].get("id") == "OptimisedRelaxedClock.c:clock"

    rates = root.xpath(".//parameter[@id='ORCRates.c:clock']")[0]
    assert rates.get("dimension") == "4"  # 2*3 - 2

    assert root.xpath(".//parameter[@id='ucldMean.c:clock']")
    samplers = root.xpath(".//operator[@spec='AdaptableOperatorSampler']")
    assert len(samplers) == 4
    metric = root.xpath(".//metric[@spec='beastlabs.evolution.tree.RNNIMetric']")[0]
    assert metric.get("taxonset") == "@TaxonSet.overall"
    tree_logger = root.xpath(".//log[contains(@spec, 'TreeWithMetaDataLogger')]")[0]
    assert tree_logger.get("branchratemodel") == "@OptimisedRelaxedClock.c:clock"
    assert _dup_ids(root) == []


def test_fixed_source_mean_carried_as_fixed():
    # overall-ctmc fixes clockRate (estimate=false, 1.0): ORC mean stays fixed,
    # so no mean scaler / prior / log, but per-branch ORCRates are still estimated.
    root = _load()
    convert_clock(root, "orc")
    ucld = root.xpath(".//parameter[@id='ucldMean.c:clock']")[0]
    assert ucld.get("estimate") == "false"
    assert root.xpath(".//operator[@id='ucldMeanScaler.c:clock']") == []
    assert root.xpath(".//prior[@id='ucldMeanPrior.c:clock']") == []


def test_fixed_mean_inlined_not_in_state():
    # a fixed mean must not sit in <state> (BEAST warns: no operator); it is
    # inlined as the branchRateModel's clock.rate child instead.
    root = _load()
    convert_clock(root, "orc")
    state = root.xpath(".//state[@id='state']")[0]
    assert state.xpath("./parameter[@id='ucldMean.c:clock']") == []
    brm = root.xpath(".//branchRateModel")[0]
    child = brm.xpath("./parameter[@name='clock.rate']")[0]
    assert child.get("id") == "ucldMean.c:clock"
    assert child.get("estimate") == "false"


def test_mean_always_logged_even_when_fixed():
    # fixed source clock -> mean still appears in the trace log
    root = _load()
    convert_clock(root, "orc")
    tracelog = root.xpath(".//logger[@id='tracelog']")[0]
    assert tracelog.xpath("./log[@idref='ucldMean.c:clock']")

    root = _load()
    convert_clock(root, "strict")
    tracelog = root.xpath(".//logger[@id='tracelog']")[0]
    assert tracelog.xpath("./log[@idref='clockRate.c:clock']")


def test_orc_back_to_strict_roundtrip():
    root = _load()
    convert_clock(root, "orc")
    convert_clock(root, "strict")
    assert detect_clock(root) == "strict"
    assert root.xpath(".//*[starts-with(@id, 'ORC')]") == []
    assert len(root.xpath(".//branchRateModel")) == 1
    assert _dup_ids(root) == []


def test_strict_to_strict_idempotent():
    root = _load()
    convert_clock(root, "strict")
    assert detect_clock(root) == "strict"
    assert len(root.xpath(".//branchRateModel")) == 1
    assert root.xpath(".//parameter[@id='clockRate.c:clock']")
    assert _dup_ids(root) == []


def test_detect_strict():
    root = etree.fromstring(STRICT_XML)
    assert detect_clock(root) == "strict"


def test_detect_orc():
    root = etree.fromstring(ORC_XML)
    assert detect_clock(root) == "orc"


CHILD_FIXED_XML = """
<beast><run id="mcmc">
  <state id="state"/>
  <distribution id="likelihood">
    <distribution id="treeLikelihood.x" spec="TreeLikelihood" tree="@Tree.t:tree">
      <branchRateModel id="StrictClock.c:clock"
        spec="beast.base.evolution.branchratemodel.StrictClockModel">
        <parameter id="clockRate.c:clock" spec="parameter.RealParameter"
          estimate="false" name="clock.rate">1.0</parameter>
      </branchRateModel>
    </distribution>
  </distribution>
</run></beast>
"""


def test_extract_mean_estimated_state_param():
    root = etree.fromstring(STRICT_XML)
    mean = extract_mean(root)
    assert mean.estimate is True
    assert mean.value == "1.0"
    assert mean.lower == "0.0"
    assert mean.upper == "1.0"
    assert mean.distr is not None
    assert mean.distr.tag == "LogNormal"
    assert mean.distr.get("M") == "-9.21"


def test_extract_mean_fixed_child_param():
    root = etree.fromstring(CHILD_FIXED_XML)
    mean = extract_mean(root)
    assert mean.estimate is False
    assert mean.value == "1.0"
    assert mean.distr is None


def test_strip_removes_all_clock_ids():
    root = etree.fromstring(ORC_XML)
    strip_clock(root)
    assert root.xpath(".//*[contains(@id, '.c:clock')]") == []


def test_strip_removes_branchratemodel():
    root = etree.fromstring(STRICT_XML)
    strip_clock(root)
    assert root.xpath(".//branchRateModel") == []


def test_strip_removes_dangling_idref_logs():
    # STRICT_XML logs the estimated clock rate via <log idref="clockRate.c:clock"/>;
    # that reference must go or BEAST fails with "Could not find object".
    root = etree.fromstring(STRICT_XML)
    strip_clock(root)
    assert root.xpath(".//*[contains(@idref, '.c:clock')]") == []


def test_strip_clears_dangling_branchratemodel_attr():
    xml = """
    <beast><run id="mcmc">
      <distribution id="likelihood">
        <distribution id="treeLikelihood.x" spec="TreeLikelihood" tree="@Tree.t:tree">
          <branchRateModel id="StrictClock.c:clock"
            spec="beast.base.evolution.branchratemodel.StrictClockModel"
            clock.rate="@clockRate.c:clock"/>
        </distribution>
      </distribution>
      <logger id="treelog">
        <log id="TreeWithMetaDataLogger.t:tree" spec="TreeWithMetaDataLogger"
             branchratemodel="@StrictClock.c:clock" tree="@Tree.t:tree"/>
      </logger>
    </run></beast>
    """
    root = etree.fromstring(xml)
    strip_clock(root)
    log = root.xpath(".//log[@id='TreeWithMetaDataLogger.t:tree']")[0]
    assert log.get("branchratemodel") is None
    # the surviving anchor element is untouched
    assert root.xpath(".//log[@id='TreeWithMetaDataLogger.t:tree']")
