"""Miscellaneous tests for files encountered in wild etc"""

from pathlib import Path

from beastwords.main import Converter


def test_misc1():
    # ascertainment character was not what I expected
    conv = Converter.from_file(Path(__file__).parent / 'misc1.xml')
    words = conv.get_words()  # just a list of (charname, type)
    assert words[0] == ('_ascertainment', 'UserDataType.0')
    partitions, ascertainment = conv.get_partitions()
    assert partitions == {'I': [1]}  # no _ascertainment! 
    assert ascertainment == [0]

