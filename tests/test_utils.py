import pytest

from beastwords.utils import repartition_by_size, repartition_by_groupsize, _split

@pytest.fixture
def data():
    return {                         # cumulative   
        'book': [16],                # 1            
        'elbow': [15],               # 2            
        'hand': [1, 2],              # 4            
        'eye': [7, 8, 9],            # 7            mid
        'foot': [3, 4, 5, 6],        # 11
        'arm': [10, 11, 12, 13, 14], # 16
    }  # n.b. this is sorted by set size to help debug.


def test_repartition_by_size(data):
    
    assert repartition_by_size(1, data) == {
        'p1': data['book'] + data['elbow'] + data['hand'] + data['eye'] + data['foot'] + data['arm']
    }
    
    assert repartition_by_size(2, data) == {
        'p1': data['book'] + data['elbow'] + data['hand'] + data['eye'],
        'p2': data['foot'] + data['arm']
    }
    
    assert repartition_by_size(3, data) == {
        'p1': data['book'] + data['elbow']  + data['hand'],   # 4
        'p2': data['eye'],                                     # 3
        'p3': data['foot'] + data['arm']                      # 9
    }
    
    assert repartition_by_size(4, data) == {
        'p1': data['book'] + data['elbow'] + data['hand'],    # 4
        'p2': data['eye'],                                    # 3
        'p3': data['foot'],                                   # 4
        'p4': data['arm']                                     # 5
    }    

    with pytest.raises(ValueError):
        repartition_by_size(10, data)


def test__split():
    assert _split("1-2") == [[1,2]]
    assert _split("4-9") == [[4,5,6,7,8,9]]
    assert _split("1,3-5,9") == [[1],[3,4,5],[9]]
    assert _split("1-5,6-10,11-15,16-20,25") == [
        [1, 2, 3, 4, 5], [6, 7, 8, 9, 10], [11, 12, 13, 14, 15],
        [16, 17, 18, 19, 20], [25]
    ]


def test_repartition_by_groupsize(data):
    
    with pytest.warns(UserWarning):
        assert repartition_by_groupsize("1,2", data) == {
            'p1': data['book'] + data['elbow'],
            'p2': data['hand']
        }
    
    # no warning
    assert repartition_by_groupsize("1-3,4-9", data) == {
        'p1-3': data['book'] + data['elbow'] + data['hand'] + data['eye'],
        'p4-9': data['foot'] + data['arm']
    }
    
    # no warning
    assert repartition_by_groupsize("1,2,3,4-10", data) == {
        'p1': data['book'] + data['elbow'],
        'p2': data['hand'],
        'p3': data['eye'],
        'p4-10': data['foot'] + data['arm']
    }
    
    # check errors on overlap
    with pytest.raises(ValueError):
        repartition_by_groupsize("1-2,2-9", data)
