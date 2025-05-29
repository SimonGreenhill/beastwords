from copy import deepcopy
from collections import defaultdict
from pathlib import Path

from lxml import etree

from warnings import warn

from beastwords.utils import repartition_by_size, repartition_by_groupsize


class Converter(object):
    
    userDataType_spec = '?'
    useAmbiguities = 'false'
    
    def __init__(self, xmlfile, tree=None, root=None, model=None):
        if not xmlfile.exists():
            raise IOError(f"File {xmlfile} does not exist")
        self.xmlfile = xmlfile
        self.tree = tree if tree is not None else etree.parse(xmlfile)
        self.root = root if root is not None else self.tree.getroot()
        self.model = model if model is not None else self.root.get("beautitemplate")

        self.words = self.get_words()
        self.partitions, self.ascertainment = self.get_partitions()
    
    @classmethod
    def from_file(cls, xmlfile):
        tree = etree.parse(xmlfile)
        root = tree.getroot()
        model = root.get("beautitemplate")
        
        if model == "BinaryCovarion":
            return CovarionConverter(xmlfile, tree=tree, root=root, model=model)
        elif model == "BinaryCTMC":
            return CTMCConverter(xmlfile, tree=tree, root=root, model=model)
        else:
            warn(f"Unsupported beauti template: {model}")
            return Converter(xmlfile, tree=tree, root=root, model=model)
    
    def get_words(self):
        words = []
        for e in self.root.findall(".//charstatelabels"):
            words.append((e.get('characterName'), e.get('id')))
        return words
    
    def patch(self, element, newattrib={}, update=False):
        """
        Clones `element`, updating it with key:values in `newattrib`
        
        If `update` is False - returns a clone, if true then it alters the original
        """
        new = deepcopy(element) if not update else element
        for key, value in newattrib.items():
            new.set(key, value)
        return new
        
    def patch_child_ids(self, element, partition):
        """Iterates over children ids and adds the partition name to their IDs"""
        new = deepcopy(element)
        # do root
        old_id = new.get("id").split(":")[0]
        new.set("id", f"{old_id}:{partition}")
        # now do children
        for el in new.xpath(".//*[@id]"):
            old_id = el.get("id").split(":")[0]
            el.set("id", f"{old_id}:{partition}")
        return new

    def replace(self, xpath, **kwargs):
        """
        Replaces a single element with one for each partition, setting values to kwargs
        """
        old = self.root.xpath(xpath)
        if len(old) == 0:
            raise ValueError(f"Can't find element: {xpath}")
        elif len(old) > 1:
            raise ValueError(f"Found many elements: {xpath}")
        
        for p in sorted(self.partitions):
            attr = {k: v.format(p) for (k, v) in kwargs.items()}
            new = self.patch(old[0], newattrib=attr, update=False)
            parent = old[0].getparent()
            index = parent.index(old[0])
            parent.insert(index + 1, new)
        old[0].getparent().remove(old[0])  # remove old one
        
    def parse_word(self, w):
        return w.replace("_u_", "_").rsplit("_" ,1)
    
    def set_partitions(self, size):
        try:
            size = int(size)
            self.partitions = repartition_by_size(size, self.partitions)
        except ValueError:
            self.partitions = repartition_by_groupsize(size, self.partitions)
        except:
            raise
    
    def get_partitions(self):
        partitions = defaultdict(list)
        for i, (char, _id) in enumerate(self.words, 0):
            partitions[self.parse_word(char)[0]].append(i)
        ascertainment = partitions.pop("_ascertainment", [])
        return (partitions, ascertainment)

    def get_gamma(self):
        return 1  # assuming we don't want a gamma per partition here.

    def get_partition_range(self, partition):
        """Returns a string showing the range of sites in this partition"""
        sites = sorted(self.partitions.get(partition, []))
        if not len(sites):
            return ""
        runs = []
        start = prev = sites[0]
        for s in sites[1:]:
            if s == prev + 1:
                prev = s
            else:
                # finalize the current run
                runs.append(f"{start}-{prev}" if start != prev else f"{start}")
                start = prev = s
        # add the final run
        runs.append(f"{start}-{prev}" if start != prev else f"{start}")
        return ",".join(runs)

    def _convert_sequences(self):  # i.e. add ascertainment characters into each partition
        # get sequences
        sequences = {s.get('id'): s.get('value') for s in self.root.xpath('.//sequence')}
        # convert sequences in XML to dictionary of {'partition': {'taxon1': '...', 'taxon2': '...'}}
        # n.b. this will ignore the old 'ascertainment' character (effectively deleting it) 
        # as it's not in the list of partitions
        new = defaultdict(lambda: defaultdict(list))
        for taxon in sequences:
            for partition, sites in self.partitions.items():
                for s in sites:
                    new[partition][taxon].append(sequences[taxon][s])
        
        # loop over partitions, create an ascertainment character and insert it at pos 0
        for partition in new:
            for taxon in new[partition]:
                # identify ascertainment character
                if all(c == '?' for c in new[partition][taxon]):
                    ascertainment = '?'
                elif all(c == '-' for c in new[partition][taxon]):
                    ascertainment = '-'
                else:
                    ascertainment = '0'
                new[partition][taxon].insert(0, ascertainment)
        
        # ok, now regenerate sequences and figure out positions
        positions = [] # we have `taxon` initialised above, we'll count that one
        for oldseq in self.root.xpath('.//sequence'):
            seqid = oldseq.get('id')
            newseq = etree.Element("sequence", id=seqid, taxon=oldseq.get('taxon'), spec="Sequence", totalcount="2")
            value = []
            for partition in sorted(new):
                value.append("".join(new[partition][seqid]))
                
                if oldseq.get('id') == taxon:
                    positions.extend([
                        (partition, i) for i, _ in enumerate(new[partition][seqid])
                    ])
            
            newseq.set('value', " ".join(value))
            oldseq.getparent().append(newseq)  # add new sequence
            oldseq.getparent().remove(oldseq)  # remove old seq
        
        # generate userDataType -- find old userDataType, and update
        # while we're here we will update ascertainment/partitions
        self.partitions, self.ascertainment = defaultdict(list), []
        udt = self.root.xpath('./data/userDataType')[0]
        # remove old chars
        for o in udt.getchildren():
            udt.remove(o)
        
        for i, (char, index) in enumerate(positions, 1):
            o = etree.Element("charstatelabels",
                id=f"UserDataType.{i}",
                spec="beast.base.evolution.datatype.UserDataType",
                characterName=f"{char}_{index}",
                codeMap="", states="-1", value="")
            udt.append(o)
            
            self.partitions[char].append(i)
            if index == 0:
                self.ascertainment.append(i)
        
    def _convert_state(self):
        path = ".//state[@id='state']/parameter[starts-with(@id, 'mutationRate.s:')]"
        mr = self.root.xpath(path)
        if len(mr) == 0: # Simon likes to delete these from one partiton runs. Make one up
            #<parameter id="mutationRate.s:overall" spec="parameter.RealParameter" name="stateNode">1.0</parameter>
            p = etree.Element("parameter",
                id="mutationRate.s:dummy", spec="parameter.RealParameter", name="stateNode")
            p.text = "1.0"
            self.root.xpath(".//state[@id='state']")[0].append(p)
        self.replace(path, id="mutationRate.s:{}")

    def _convert_prior(self):
        prior = self.root.xpath(".//distribution[@id='prior']")[0]
        path = ".//prior[starts-with(@id, 'MutationRatePrior.s:')]"
        mrp = prior.xpath(path)
        if len(mrp) == 0: # Simon likes to delete these from one partiton runs. Make one up
            mrp = etree.Element("prior",
                id="MutationRatePrior.s:dummy", name="distribution", x="@mutationRate.s:dummy")
            etree.SubElement(mrp, "OneOnX", id="OneOnX.0", name="distr")
            prior.append(mrp)
        
        self.replace(path, id="MutationRatePrior.s:{}", x="@mutationRate.s:{}")
        # and update internal OneOnX
        for o in self.root.xpath(path):
            p = o.get('id').split(":")[1]
            o.getchildren()[0].set('id', f"OneOnX:{p}")
            
    def _add_substmodel(self, partition, siteModel):
        return siteModel

    def _convert_treelikelihood(self):
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
        
        # find data/sequence
        data = self.root.xpath('.//data[@spec="FilteredAlignment"]')
        if len(data) > 1:
            raise ValueError("I can't handle multiple partitions")
        data = data[0]
        seq = data.get('data')

        # find brm
        brm = self.root.xpath('.//branchRateModel')[0]
        assert brm is not None, "Unable to find branchRateModel"
        brm_id = brm.get('id')
        # -> we will move this into the first partition later

        # find tree
        tree = self.root.xpath('.//init')[0]
        assert tree is not None, "Unable to find tree"
        tree_id = tree.get('initial')
        
        # find substModel
        substModel = self.root.xpath(".//substModel")[0]
        assert substModel is not None, "Unable to find substModel"
        
        # find Lh and treeLh
        likelihood = self.root.xpath(f".//distribution[@id='likelihood']")[0]
        assert likelihood is not None, "Unable to find likelihood"
        treeLh = likelihood.xpath(".//distribution[@spec='TreeLikelihood']")[0]
        
        # add the required substModels and put them after the state
        state = self.root.xpath(f".//state[@id='state']")[0]
        
        for i, p in enumerate(self.partitions):
            # 1. construct <distribution>
            distribution = etree.Element("distribution",
                id=f"treeLikelihood.{p}",
                spec="TreeLikelihood",
                tree=f"{tree_id}",
                useAmbiguities=self.useAmbiguities)
            
            # add the branch rate model to the first partition
            if i == 0:
                distribution.append(self.patch(brm))  # use patch just to clone
                brm.getparent().remove(brm)  # now delete brm
            else:
                distribution.set('branchRateModel', f"@{brm_id}") 
                
            # 2. construct <data>
            chars = self.get_partition_range(p)
            d1 = etree.Element("data",
                id=f"data:{p}", spec="FilteredAlignment", ascertained="true", excludeto="1", 
                filter="-"
            )
            d2 = etree.SubElement(d1, "data", id=f"{p}", spec="FilteredAlignment", data=f"{seq}",
                filter=self.get_partition_range(p)
            )
            udt = etree.SubElement(d2, "userDataType", id=f"userDataType:{p}", spec=self.userDataType_spec)
            
            distribution.append(d1)
            
            # add siteModel
            siteModel = etree.Element("siteModel",
                id=f"SiteModel.s:{p}",
                spec="SiteModel",
                gammaCategoryCount="%d" % self.get_gamma(),
                mutationRate=f"@mutationRate.s:{p}")
            
            # add substModel
            self._add_substmodel(p, siteModel)

            p2 = etree.SubElement(siteModel, "parameter",
                id=f"proportionInvariant.s:{p}", spec="parameter.RealParameter", estimate="false",
                lower="0.0", name="proportionInvariant", upper="1.0"
            )
            p2.text = '0.0'
            
            distribution.append(siteModel)
            likelihood.append(distribution)
        
        # cleanup old stuff.
        data.getparent().remove(data)
        treeLh.getparent().remove(treeLh)
        substModel.getparent().remove(substModel)

    def _convert_operators(self):
        path = ".//operator[starts-with(@id, 'mutationRateScaler.s:')]"
        mrs = self.root.xpath(path)
        if len(mrs) == 0: # Simon likes to delete these from one partiton runs. Make one up
            # <operator id="mutationRateScaler.s:hand" spec="ScaleOperator" parameter="@mutationRate.s:hand" scaleFactor="0.5" weight="0.1"/>
            mrs = etree.Element("operator",
                id="mutationRateScaler.s:dummy", spec="ScaleOperator", parameter="@mutationRate.s:dummy",
                scaleFactor="0.5", weight="0.1")

            # find last operator
            last = self.root.xpath(".//operator")[-1]
            last.getparent().append(mrs)
            parent = last.getparent()
            index = parent.index(last)
            parent.insert(index + 1, mrs)
            
        self.replace(path, id="mutationRateScaler.s:{}", parameter="@mutationRate.s:{}")

    def _convert_log(self):
        # <log idref="treeLikelihood.hand"/>
        self.replace(".//log[starts-with(@idref, 'treeLikelihood.')]", idref="treeLikelihood.{}")
        # <log idref="mutationRate.s:hand"/>
        self.replace(".//log[starts-with(@idref, 'mutationRate.s')]", idref="mutationRate.s:{}")
    
    def convert(self):
        self._convert_sequences() # should go first i think
        self._convert_state()
        self._convert_prior()
        self._convert_treelikelihood()
        self._convert_operators()
        self._convert_log()
        
    def __str__(self):
        return self.write(self.tree)
        
    def write(self, el=None):
        el = el if el is not None else self.tree
        etree.indent(el)  # needed to 'reset' the indentation
        return etree.tostring(el, pretty_print=True, encoding='unicode')
    
    def to_file(self, filename):
        etree.indent(self.tree)  # needed to 'reset' the indentation
        self.tree.write(
            filename,
            xml_declaration=True,
            encoding="UTF-8",
            standalone="no",
            pretty_print=True
        )


class CovarionConverter(Converter):
    
    userDataType_spec = "beast.base.evolution.datatype.TwoStateCovarion"
    useAmbiguities = 'true'
    has_siteModel = False
    
    def _convert_state(self):
        super()._convert_state()
        el = self.root.xpath(".//parameter[starts-with(@id, 'bcov_alpha.s:')]")[0]
        assert el is not None, 'Unable to find original parameter/bcov_alpha.s:<.*>'
        el.set('id', "bcov_alpha.s:combined")
        
        el = self.root.xpath(".//parameter[starts-with(@id, 'bcov_s.s:')]")[0]
        assert el is not None, 'Unable to find original parameter/bcov_s.s:<.*>'
        el.set('id', "bcov_s.s:combined")
        
        el = self.root.xpath(".//parameter[starts-with(@id, 'frequencies.s:')]")[0]
        assert el is not None, 'Unable to find original parameter/frequencies.s:<.*>'
        el.set('id', "frequencies.s:combined")
        return el
    
    def _add_substmodel(self, partition, siteModel):
        if self.has_siteModel: # just add an attribute
            siteModel.set("substModel", "@covarion:combined")
        else:  # add the whole model
            # 1. find old one
            old = self.root.xpath(".//*/substModel")
            assert len(old) == 1, f"Expected 1 substModel, got {old}"
            new = self.patch(old[0], {
                'id': 'covarion:combined',
                'alpha': "@bcov_alpha.s:combined",
                'switchRate': "@bcov_s.s:combined",
                'vfrequencies': "@frequencies.s:combined",
            })
            new = self.patch_child_ids(new, 'combined')
            siteModel.insert(0, new)

            # add gammaShape parameter
            gammaShape = etree.SubElement(siteModel, "parameter",
               id=f"gammaShape.s:{partition}", spec="parameter.RealParameter", estimate="false", name="shape")
            gammaShape.text = '1.0'
            siteModel.append(gammaShape)

            self.has_siteModel = True  # set flag
        return siteModel
    
    def _convert_prior(self):
        super()._convert_prior()
        self.patch(
            self.root.xpath(".//prior[starts-with(@id, 'bcov_alpha_prior.s:')]")[0],
            {'id': 'bcov_alpha_prior.s:combined', 'x': "@bcov_alpha.s:combined"},
            update=True)
        self.patch(
            self.root.xpath(".//prior[starts-with(@id, 'bcov_s_prior.s:')]")[0],
            {'id': 'bcov_s_prior.s:combined', 'x': "@bcov_s.s:combined"},
            update=True)

    def _convert_operators(self):
        super()._convert_operators()
        op = self.patch(
            self.root.xpath(".//operator[starts-with(@id, 'bcovAlphaScaler.s:')]")[0],
            {'id': 'bcovAlphaScaler.s:combined', 'parameter': "@bcov_alpha.s:combined"},
            update=True)
            
        op = self.patch(
            self.root.xpath(".//operator[starts-with(@id, 'bcovSwitchParamScaler.s:')]")[0],
            {'id': 'bcovSwitchParamScaler.s:combined', 'parameter': "@bcov_s.s:combined"},
            update=True)
        
        op = self.patch(
            self.root.xpath(".//operator[starts-with(@id, 'frequenciesDelta.s:')]")[0],
            {'id': 'frequenciesDelta.s:combined'},
            update=True)
        self.patch(op.getchildren()[0], {'idref': "frequencies.s:combined"}, update=True)
    
    def _convert_log(self):
        super()._convert_log()
        self.patch(
            self.root.xpath(".//log[starts-with(@idref, 'bcov_alpha.s:')]")[0],
            {'idref': 'bcov_alpha.s:combined'},
            update=True)

        self.patch(
            self.root.xpath(".//log[starts-with(@idref, 'bcov_s.s:')]")[0],
            {'idref': 'bcov_s.s:combined'},
            update=True)
        
        self.patch(
            self.root.xpath(".//log[starts-with(@idref, 'frequencies.s:')]")[0],
            {'idref': 'frequencies.s:combined'},
            update=True)


class CTMCConverter(Converter):
    
    userDataType_spec = "beast.base.evolution.datatype.Binary"
    useAmbiguities = 'false'
    
    def _convert_state(self):
        super()._convert_state()
        # already have mutationRate.*
        # add gammaShape & freqParameter
        try:
            self.replace(".//parameter[starts-with(@id, 'gammaShape.s:')]", id="gammaShape.s:{}")
        except ValueError:
            pass
            
        self.replace(".//parameter[starts-with(@id, 'freqParameter.s:')]", id="freqParameter.s:{}")


    def _convert_prior(self):
        super()._convert_prior()
        # already have MutationRatePrior.*, add GammaShapePrior
        # <prior id="GammaShapePrior.s:foot" name="distribution" x="@gammaShape.s:foot">
        #     <Exponential id="Exponential.1" name="distr">
        #         <mean id="Function$Constant.1" spec="Function$Constant" value="1.0"/>
        #     </Exponential>
        # </prior>
        path = ".//prior[starts-with(@id, 'GammaShapePrior.s:')]"
        
        try:
            self.replace(path, id="GammaShapePrior.s:{}", x="@gammaShape.s:{}")
        except ValueError:
            return  # no gamma - nothing to do
        
        # and update internal Exponential
        for o in self.root.xpath(path):
            p = o.get('id').split(":")[1]
            exp = o.getchildren()[0]
            old_id = exp.get('id').split(".")[0]
            exp.set('id', f"{old_id}:{p}")
            # and the nested <mean>
            old_id = exp.getchildren()[0].get('id').split(".")[0]
            exp.getchildren()[0].set('id', f"{old_id}:{p}")

    def _add_substmodel(self, partition, siteModel):
        # ctmc gets one substModel per word
        old = self.root.xpath(".//*/substModel")  # don't care which one it is 
        new = self.patch_child_ids(old[0], partition)
        # get freq subelement and change frequencies="@freqParameter.s:.."
        #  <frequencies id="estimatedFreqs.s:eye" spec="Frequencies" frequencies="@freqParameter.s:overall"/>
        f = new.xpath('frequencies')[0]
        f.set('frequencies', f'@freqParameter.s:{partition}')
        siteModel.insert(0, new)
        
        # handle gamma by removing it being estimated
        if 'shape' in siteModel.attrib:
            del(siteModel.attrib['shape'])
        
        gammaShape = etree.SubElement(siteModel, "parameter",
            id=f"gammaShape.s:{partition}", spec="parameter.RealParameter", estimate="false", name="shape")
        gammaShape.text = '1.0'
        siteModel.append(gammaShape)
        return siteModel
        
    def _convert_operators(self):
        super()._convert_operators()
        
        try:
            self.replace(".//operator[starts-with(@id, 'gammaShapeScaler.s:')]", id="gammaShapeScaler.s:{}")
        except:
            pass  # no gamma
            
        self.replace(".//operator[starts-with(@id, 'FrequenciesExchanger.s:')]", id="FrequenciesExchanger.s:{}")

        # patch internal freqParameters and gammaShapeScaler
        for p in self.partitions:
            op = self.root.xpath(f".//operator[@id='FrequenciesExchanger.s:{p}']")[0]
            self.patch(op.getchildren()[0], {'idref': f"freqParameter.s:{p}"}, update=True)
            
            try:
                op = self.root.xpath(f".//operator[@id='gammaShapeScaler.s:{p}']")[0]
                op.set("parameter", f"@gammaShape.s:{p}")
            except:
                pass  # no gamma

    def _convert_log(self):
        super()._convert_log()
        # <log idref="freqParameter.s:hand"/>
        self.replace(".//log[starts-with(@idref, 'freqParameter.s')]", idref="freqParameter.s:{}")
        # <log idref="gammaShape.s:overall"/>
        try:
            self.replace(".//log[starts-with(@idref, 'gammaShape.s')]", idref="gammaShape.s:{}")
        except ValueError:
            pass  # no gamma

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Converts a one partition XML to a partitioned one')
    parser.add_argument("input", help='filename', type=Path)
    parser.add_argument("output", help='filename', type=Path)
    parser.add_argument(
        '-p', "--partitions", dest='partitions', default=None, type=str,
        help="set partition number. If this is None use words", action='store'
    )
    args = parser.parse_args()
    
    xml = Converter.from_file(args.input)
    if args.partitions:
        xml.set_partitions(args.partitions)
    xml.convert()
    xml.to_file(args.output)

if __name__ == "__main__":
    main()


