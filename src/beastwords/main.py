from copy import deepcopy
from collections import defaultdict
from pathlib import Path

from lxml import etree

from warnings import warn

        
class Converter(object):
    
    userDataType_spec = '?'
    
    def __init__(self, xmlfile, tree=None, root=None, model=None):
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
    
    def patch(self, element, newattrib, update=False):
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
        
        for p in self.partitions:
            attr = {k: v.format(p) for (k, v) in kwargs.items()}
            new = self.patch(old[0], newattrib=attr, update=False)
            parent = old[0].getparent()
            index = parent.index(old[0])
            parent.insert(index + 1, new)
        old[0].getparent().remove(old[0])  # remove old one
        
    def parse_word(self, w):
        return w.replace("_u_", "_").rsplit("_" ,1)
    
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

    def _convert_state(self):
        self.replace(".//parameter[starts-with(@id, 'mutationRate.s:')]", id="mutationRate.s:{}")

    def _convert_prior(self):
        self.replace(
            ".//prior[starts-with(@id, 'MutationRatePrior.s:')]", 
            id="MutationRatePrior.s:{}", x="@mutationRate.s:{}")

    def _convert_substmodel(self):
        pass

    def _convert_treelikelihood(self):
        # find data/sequence
        data = self.root.xpath('.//data[@spec="FilteredAlignment"]')
        if len(data) > 1:
            raise ValueError("I can't handle multiple partitions")
        data = data[0]
        seq = data.get('data')

        # find brm
        brm = self.root.xpath('.//branchRateModel')[0].get('id')
        assert brm is not None, "Unable to find branchRateModel"

        # find tree
        tree = self.root.xpath('.//init')[0].get('initial')
        assert tree is not None, "Unable to find tree"
        
        # find substModel
        substModel = self.root.xpath(".//substModel")[0]
        assert substModel is not None, "Unable to find substModel"
        
        # find Lh and treeLh
        likelihood = self.root.xpath(f".//distribution[@id='likelihood']")[0]
        assert likelihood is not None, "Unable to find likelihood"
        treeLh = likelihood.xpath(".//distribution[@spec='TreeLikelihood']")[0]
        
        # add the required substModels and put them after the state
        state = self.root.xpath(f".//state[@id='state']")[0]
        substmodel_ids = []
        for s in self._convert_substmodel():
            parent = state.getparent()
            index = parent.index(state)
            parent.insert(index + 1, s)
            substmodel_ids.append(s.get('id'))
            
        # remove the old one
        substModel.getparent().remove(substModel)
        
        for i, p in enumerate(self.partitions):
            # 1. construct <distribution>
            distribution = etree.Element("distribution",
                id=f"treeLikelihood.{p}", 
                branchRateModel=brm,
                tree=tree,
                useAmbiguities='true')
            
            # 2. construct <data>
            chars = self.get_partition_range(p)
            d1 = etree.Element("data",
                id=f"data:{p}", spec="FilteredAlignment", ascertained="true", excludeto="1", 
                filter="-"
            )
            d2 = etree.SubElement(d1, "data", id=f"{p}", spec="FilteredAlignment", data=f"{seq}",
                filter=self.get_partition_range(p)
            )
            udt = etree.SubElement(d2, "userDataType", id=f"{p}", spec=self.userDataType_spec)
            
            distribution.append(d1)

            # figure out id of substmodel
            if len(substmodel_ids) == 1:    # covarion
                smid = substmodel_ids[0] 
            elif len(substmodel_ids) == len(self.partitions):  # one per partition
                smid = substmodel_ids[i]
            else:
                raise ValueError("Substitution Models != Partition count")
            
            # add siteModel
            siteModel = etree.Element("siteModel",
                id=f"SiteModel.s:{p}",
                spec="SiteModel",
                gammaCategoryCount="%d" % self.get_gamma(),
                mutationRate="@mutationRate.s:{p}",
                substModel=f"@{smid}")
            
            p1 = etree.SubElement(siteModel, "parameter",
                id=f"gammaShape.s:{p}", spec="parameter.RealParameter", estimate="false", name="shape")
            p1.text = '1.0'

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

    def _convert_operators(self):
        self.replace(
            ".//operator[starts-with(@id, 'mutationRateScaler.s:')]", 
            id="mutationRateScaler.s:{}", parameter="@mutationRate.s:{}")

    def _convert_log(self):
        # <log idref="treeLikelihood.hand"/>
        self.replace(".//log[starts-with(@idref, 'treeLikelihood.')]", idref="treeLikelihood.{}")
        # <log idref="mutationRate.s:hand"/>
        self.replace(".//log[starts-with(@idref, 'mutationRate.s')]", idref="mutationRate.s:{}")
    
    def convert(self):
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


class CovarionConverter(Converter):
    
    userDataType_spec = "beast.base.evolution.datatype.TwoStateCovarion"
    
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
    
    def _convert_substmodel(self):
        sm = self.root.xpath(".//substModel")[0]
        sm = self.patch(sm, {
            'id': 'covarion:combined',
            'alpha': '@bcov_alpha.s:combined',
            'switchRate': '@bcov_s.s:combined',
            'vfrequencies': '@frequencies.s:combined'
        })
        sm = self.patch_child_ids(sm, 'combined')
        yield sm

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
        self.patch(op.getchildren()[0], {'idref': "@frequencies.s:combined"}, update=True)
    
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

    def _convert_state(self):
        super()._convert_state()
        # already have mutationRate.*
        # add gammaShape & freqParameter
        self.replace(".//parameter[starts-with(@id, 'gammaShape.s:')]", id="gammaShape.s:{}")
        self.replace(".//parameter[starts-with(@id, 'freqParameter.s:')]", id="freqParameter.s:{}")


    def _convert_prior(self):
        super()._convert_prior()
        # already have MutationRatePrior.*, add GammaShapePrior
        self.replace(".//prior[starts-with(@id, 'GammaShapePrior.s:')]", id="GammaShapePrior.s:{}")

    def _convert_substmodel(self):
        sm = self.root.xpath(".//substModel")[0]
        # ctmc gets one substModel per word
        for p in self.partitions:
            new = self.patch_child_ids(sm, p)
            # get freq subelement and change frequencies="@freqParameter.s:.."
            f = new.xpath('frequencies')[0]
            f.set('frequencies', f'@freqParameter.s:{p}')
            yield new
    
    def _convert_operators(self):
        super()._convert_operators()

        self.replace(".//operator[starts-with(@id, 'FrequenciesExchanger.s:')]", id="FrequenciesExchanger.s:{}")
        # patch internal freqParameters
        for p in self.partitions:
            op = self.root.xpath(f".//operator[@id='FrequenciesExchanger.s:{p}']")[0]
            self.patch(op.getchildren()[0], {'idref': f"freqParameter.s:{p}"}, update=True)

    def _convert_log(self):
        super()._convert_log()
        # <log idref="freqParameter.s:hand"/>
        self.replace(".//log[starts-with(@idref, 'freqParameter.s')]", idref="freqParameter.s:{}")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Converts a one partition XML to a partitioned one')
    parser.add_argument("input", help='filename', type=Path)
    parser.add_argument("output", help='filename', type=Path)
    args = parser.parse_args()
    
    xml = Converter.from_file(args.input)
    xml.convert()
    print(xml)
    print("TODO write")

if __name__ == "__main__":
    main()


