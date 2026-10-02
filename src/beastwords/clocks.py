"""Clock-model conversion for beastwords.

A BEAST2 clock is a single shared ``branchRateModel`` plus its supporting state
parameters, priors, operators and logs -- all carrying the ``.c:<clock>`` id
suffix. It is orthogonal to the substitution model and to word-partitioning, so
converting it is a separate pass over the whole XML, run before partitioning.
"""

from copy import deepcopy
from dataclasses import dataclass

from lxml import etree

CLOCK = "clock"  # the clock-model suffix: ids end in ``.c:clock``


@dataclass
class Mean:
    """The clock-rate (``clock.rate``) parameter, carried across conversions.

    ``distr`` is the inner prior distribution element (e.g. ``<LogNormal>``) or
    ``None`` when the mean is fixed and has no prior.
    """
    value: str
    estimate: bool
    lower: str | None
    upper: str | None
    distr: object  # lxml element or None


def detect_clock(root):
    """Return the clock type of ``root``: ``"strict"`` or ``"orc"``."""
    brm = root.xpath(".//branchRateModel")
    if not brm:
        raise ValueError("No branchRateModel found")
    spec = brm[0].get("spec", "")
    if spec.endswith("StrictClockModel"):
        return "strict"
    if spec.endswith("UCRelaxedClockModel"):
        return "orc"
    raise ValueError(f"Unrecognised clock spec: {spec}")


def extract_mean(root):
    """Return the clock-rate :class:`Mean` carried from the source clock.

    Handles both shapes: a ``clock.rate`` child ``<parameter>`` of the
    ``branchRateModel`` (fixed-clock inputs) and a ``clock.rate="@id"`` reference
    to a state parameter (calibrated inputs).
    """
    brm = root.xpath(".//branchRateModel")[0]
    ref = brm.get("clock.rate")
    if ref:
        param_id = ref.lstrip("@")
        param = root.xpath(f".//parameter[@id='{param_id}']")[0]
    else:
        param = brm.xpath("./parameter[@name='clock.rate']")[0]
        param_id = param.get("id")

    distr = None
    priors = root.xpath(f".//prior[@x='@{param_id}']")
    if priors:
        distr = deepcopy(priors[0][0])

    return Mean(
        value=(param.text or "").strip(),
        estimate=param.get("estimate") != "false",
        lower=param.get("lower"),
        upper=param.get("upper"),
        distr=distr,
    )


def strip_clock(root):
    """Remove every clock element and dangling clock reference from ``root``.

    Clock elements carry the ``.c:clock`` id suffix; the ``branchRateModel`` is
    removed whether or not it does. Surviving loggers keep their id but lose any
    ``branchratemodel`` attribute pointing at the old clock.
    """
    for el in root.xpath(".//*[contains(@id, '.c:clock')]"):
        el.getparent().remove(el)
    for brm in root.xpath(".//branchRateModel"):
        brm.getparent().remove(brm)
    # drop elements that only *reference* the clock (e.g. a tracelog
    # <log idref="clockRate.c:clock"/>), else BEAST fails on the dangling idref
    for el in root.xpath(".//*[contains(@idref, '.c:clock')]"):
        el.getparent().remove(el)
    for el in root.xpath(".//*[contains(@branchratemodel, '.c:clock')]"):
        del el.attrib["branchratemodel"]


def _graft(anchor, fragment):
    """Parse an XML ``fragment`` (possibly several siblings) into ``anchor``."""
    wrapper = etree.fromstring(f"<_>{fragment}</_>")
    for child in wrapper:
        anchor.append(child)


def _wire_mean(brm, mean, anchors, mean_id, prior_id, scaler, updown):
    """Attach the carried-over clock-rate mean to ``brm``.

    Estimated: add a state parameter + prior + scaler/updown operators, and
    reference it via ``clock.rate="@id"``. Fixed: inline it as the
    ``branchRateModel``'s ``clock.rate`` child so it never enters ``<state>``
    (which would warn about a state node with no operator). Either way the mean is
    added to the trace log. ``scaler``/``updown`` are ``(id, scaleFactor, weight)``
    triples.
    """
    if not mean.estimate:
        p = etree.SubElement(brm, "parameter",
            id=mean_id, spec="parameter.RealParameter",
            estimate="false", name="clock.rate")
        p.text = mean.value or "1.0"
    else:
        p = etree.SubElement(anchors["state"], "parameter",
            id=mean_id, spec="parameter.RealParameter", name="stateNode")
        if mean.lower is not None:
            p.set("lower", mean.lower)
        if mean.upper is not None:
            p.set("upper", mean.upper)
        p.text = mean.value or "1.0"
        brm.set("clock.rate", f"@{mean_id}")

        if mean.distr is not None:
            pr = etree.SubElement(anchors["prior"], "prior",
                id=prior_id, name="distribution")
            pr.set("x", f"@{mean_id}")
            pr.append(deepcopy(mean.distr))

        sid, sfac, sweight = scaler
        etree.SubElement(anchors["run"], "operator",
            id=sid, spec="ScaleOperator", parameter=f"@{mean_id}",
            scaleFactor=sfac, weight=sweight)
        uid, ufac, uweight = updown
        ud = etree.SubElement(anchors["run"], "operator",
            id=uid, spec="UpDownOperator", scaleFactor=ufac, weight=uweight)
        etree.SubElement(ud, "up", idref=mean_id)
        etree.SubElement(ud, "down", idref=anchors["tree"])

    # always log the clock-rate mean, even when fixed
    etree.SubElement(anchors["tracelog"], "log", idref=mean_id)


def _anchors(root):
    tl = root.xpath(".//distribution[@spec='TreeLikelihood']")
    tree_ref = tl[0].get("tree") if tl else root.xpath(".//branchRateModel")[0].get("tree")
    return {
        "state": root.xpath(".//state[@id='state']")[0],
        "prior": root.xpath(".//distribution[@id='prior']")[0],
        "run": root.xpath(".//run[@id='mcmc']")[0],
        "tracelog": root.xpath(".//logger[@id='tracelog']")[0],
        "tree": (tree_ref or "@Tree.t:tree").lstrip("@"),
    }


def build_strict(root, mean, brm_parent, anchors):
    mean_id = f"clockRate.c:{CLOCK}"
    brm = etree.SubElement(brm_parent, "branchRateModel",
        id=f"StrictClock.c:{CLOCK}",
        spec="beast.base.evolution.branchratemodel.StrictClockModel")
    _wire_mean(brm, mean, anchors, mean_id,
        prior_id=f"ClockPrior.c:{CLOCK}",
        scaler=(f"StrictClockRateScaler.c:{CLOCK}", "0.75", "5.0"),
        updown=(f"strictClockUpDownOperator.c:{CLOCK}", "0.75", "3.0"))


def build_orc(root, mean, brm_parent, anchors):
    mean_id = f"ucldMean.c:{CLOCK}"
    tree = anchors["tree"]
    nbranches = 2 * len(root.xpath(".//sequence")) - 2
    taxonset = root.xpath(".//taxonset[@id]")[0].get("id")
    brm_id = f"OptimisedRelaxedClock.c:{CLOCK}"

    brm = etree.fromstring(f"""
      <branchRateModel id="{brm_id}"
        spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel"
        rates="@ORCRates.c:{CLOCK}" tree="@{tree}">
        <LogNormal id="LogNormalDistributionModel.c:{CLOCK}" S="@ORCsigma.c:{CLOCK}" meanInRealSpace="true" name="distr">
          <parameter id="ORCLogNormalM.c:{CLOCK}" spec="parameter.RealParameter" estimate="false" name="M">1.0</parameter>
        </LogNormal>
      </branchRateModel>""")
    _wire_mean(brm, mean, anchors, mean_id,
        prior_id=f"ucldMeanPrior.c:{CLOCK}",
        scaler=(f"ucldMeanScaler.c:{CLOCK}", "0.5", "3.0"),
        updown=(f"relaxedUpDownOperator.c:{CLOCK}", "0.9", "15.0"))
    brm_parent.append(brm)

    _graft(anchors["state"], f"""
      <parameter id="ORCsigma.c:{CLOCK}" spec="parameter.RealParameter"
        lower="0.0" upper="1.0" name="stateNode">0.1</parameter>
      <parameter id="ORCRates.c:{CLOCK}" spec="parameter.RealParameter"
        dimension="{nbranches}" lower="1.0E-100" name="stateNode">1.0</parameter>
    """)

    _graft(anchors["prior"], f"""
      <prior id="ORCRatePriorDistribution.c:{CLOCK}" name="distribution" x="@ORCRates.c:{CLOCK}">
        <LogNormal id="ORCRatesLogNormal.c:{CLOCK}" S="@ORCsigma.c:{CLOCK}" meanInRealSpace="true" name="distr">
          <parameter id="ORCRatesLogNormalM.c:{CLOCK}" spec="parameter.RealParameter" estimate="false" name="M">1.0</parameter>
        </LogNormal>
      </prior>
      <prior id="ORCsigmaPrior.c:{CLOCK}" name="distribution" x="@ORCsigma.c:{CLOCK}">
        <Gamma id="ORCsigmaGamma.c:{CLOCK}" name="distr">
          <parameter id="ORCsigmaGammaAlpha.c:{CLOCK}" spec="parameter.RealParameter" estimate="false" name="alpha">0.5396</parameter>
          <parameter id="ORCsigmaGammaBeta.c:{CLOCK}" spec="parameter.RealParameter" estimate="false" name="beta">0.3819</parameter>
        </Gamma>
      </prior>
    """)

    _graft(anchors["run"], f"""
      <operator id="ORCAdaptableOperatorSampler_sigma.c:{CLOCK}" spec="AdaptableOperatorSampler" weight="1.0">
        <parameter idref="ORCsigma.c:{CLOCK}"/>
        <operator id="ORCucldStdevScaler.c:{CLOCK}" spec="orc.consoperators.UcldScalerOperator" distr="@LogNormalDistributionModel.c:{CLOCK}" rates="@ORCRates.c:{CLOCK}" scaleFactor="0.5" stdev="@ORCsigma.c:{CLOCK}" weight="1.0"/>
        <operator id="ORCUcldStdevRandomWalk.c:{CLOCK}" spec="operator.kernel.BactrianRandomWalkOperator" parameter="@ORCsigma.c:{CLOCK}" scaleFactor="0.1" weight="1.0"/>
        <operator id="ORCUcldStdevScale.c:{CLOCK}" spec="kernel.BactrianScaleOperator" parameter="@ORCsigma.c:{CLOCK}" scaleFactor="0.5" upper="10.0" weight="1.0"/>
        <operator id="ORCSampleFromPriorOperator_sigma.c:{CLOCK}" spec="orc.operators.SampleFromPriorOperator" parameter="@ORCsigma.c:{CLOCK}" prior2="@ORCsigmaPrior.c:{CLOCK}" weight="1.0"/>
      </operator>
      <operator id="ORCAdaptableOperatorSampler_rates_root.c:{CLOCK}" spec="AdaptableOperatorSampler" weight="1.0">
        <parameter idref="ORCRates.c:{CLOCK}"/>
        <tree idref="{tree}"/>
        <operator id="ORCRootOperator1.c:{CLOCK}" spec="orc.consoperators.SimpleDistance" clockModel="@{brm_id}" rates="@ORCRates.c:{CLOCK}" tree="@{tree}" twindowSize="0.005" weight="1.0"/>
        <operator id="ORCRootOperator2.c:{CLOCK}" spec="orc.consoperators.SmallPulley" clockModel="@{brm_id}" dwindowSize="0.005" rates="@ORCRates.c:{CLOCK}" tree="@{tree}" weight="1.0"/>
      </operator>
      <operator id="ORCAdaptableOperatorSampler_rates_internal.c:{CLOCK}" spec="AdaptableOperatorSampler" weight="5.0">
        <parameter idref="ORCRates.c:{CLOCK}"/>
        <tree idref="{tree}"/>
        <operator id="ORCInternalnodesOperator.c:{CLOCK}" spec="orc.consoperators.InConstantDistanceOperator" clockModel="@{brm_id}" rates="@ORCRates.c:{CLOCK}" tree="@{tree}" twindowSize="0.005" weight="1.0"/>
        <operator id="ORCRatesRandomWalk.c:{CLOCK}" spec="operator.kernel.BactrianRandomWalkOperator" parameter="@ORCRates.c:{CLOCK}" scaleFactor="0.1" weight="1.0"/>
        <operator id="ORCRatesScale.c:{CLOCK}" spec="kernel.BactrianScaleOperator" parameter="@ORCRates.c:{CLOCK}" scaleFactor="0.5" upper="10.0" weight="1.0"/>
        <operator id="ORCSampleFromPriorOperator_rates.c:{CLOCK}" spec="orc.operators.SampleFromPriorOperator" parameter="@ORCRates.c:{CLOCK}" prior2="@ORCRatePriorDistribution.c:{CLOCK}" weight="1.0"/>
      </operator>
      <operator id="ORCAdaptableOperatorSampler_NER.c:{CLOCK}" spec="AdaptableOperatorSampler" weight="2.0">
        <tree idref="{tree}"/>
        <operator id="ORCNER_null.c:{CLOCK}" spec="orc.operators.MetaNEROperator" rates="@ORCRates.c:{CLOCK}" tree="@{tree}" weight="1.0"/>
        <operator id="ORCNER_dAE_dBE_dCE.c:{CLOCK}" spec="orc.ner.NEROperator_dAE_dBE_dCE" rates="@ORCRates.c:{CLOCK}" tree="@{tree}" weight="1.0"/>
        <metric id="RNNIMetric.c:{CLOCK}" spec="beastlabs.evolution.tree.RNNIMetric" taxonset="@{taxonset}"/>
      </operator>
    """)

    _graft(anchors["tracelog"], f"""
      <log idref="ORCsigma.c:{CLOCK}"/>
      <log id="ORCRatesStat.c:{CLOCK}" spec="beast.base.evolution.RateStatistic" branchratemodel="@{brm_id}" tree="@{tree}"/>
    """)

    for log in root.xpath(".//log[contains(@spec, 'TreeWithMetaDataLogger')]"):
        log.set("branchratemodel", f"@{brm_id}")


_BUILDERS = {"strict": build_strict, "orc": build_orc}


def convert_clock(root, target):
    """Convert the clock model of ``root`` in place to ``target`` (strict/orc).

    Carries the source clock-rate mean (value, estimate flag, prior) across,
    then rebuilds the rate-variation machinery for the target clock. Run before
    word-partitioning so the single ``branchRateModel`` is cloned per partition
    by the normal pipeline.
    """
    if target not in _BUILDERS:
        raise ValueError(f"Unknown clock target: {target!r}")
    mean = extract_mean(root)
    brm_parent = root.xpath(".//branchRateModel")[0].getparent()
    anchors = _anchors(root)
    strip_clock(root)
    _BUILDERS[target](root, mean, brm_parent, anchors)
