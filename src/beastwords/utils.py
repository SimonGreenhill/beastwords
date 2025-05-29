import math
from warnings import warn

def repartition_by_size(partitions, data):
    """
    Repartitions `data` by splitting into npartitions sizes
    
    > repartition_by_size("1-5,6-10", {...})
    """
    
    # 1) Sort groups by (length, key) so that equal‐sized groups
    #    are ordered by their dict key name
    if partitions > len(data.keys()):
        raise ValueError(f"too many partitions for {len(data)} words")
    
    sorted_items = sorted(data.items(), key=lambda kv: (len(kv[1]), kv[0]))
    groups = [values for _, values in sorted_items]

    # 2) Compute total items and ideal max per partition
    total = sum(len(g) for g in groups)
    ideal = math.ceil(total / partitions)

    # 3) Allocate
    result = {f'p{i+1}': [] for i in range(partitions)}
    sizes = [0] * partitions
    j = 0  # current partition index

    for group in groups:
        # if adding this group would exceed ideal, and we still
        # have more partitions left, advance to the next one
        if sizes[j] + len(group) > ideal and j + 1 < partitions:
            j += 1

        result[f'p{j+1}'].extend(group)
        sizes[j] += len(group)

    for k, v in result.items():
        if len(v) == 0:
            warn(f"Warning set {k} is empty")

    return result


def _split(rangestr):
    chunks = []
    for part in rangestr.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            chunks.append(list(range(start, end + 1)))
        else:
            chunks.append([int(part)])
    return chunks


def repartition_by_groupsize(partitions, data):
    """
    Repartitions `data` by splitting into groups
    
    > repartition_by_group("1-5,6-10", {...})
    """
    # 1. collect data by size:
    sizes = [(len(sites), k) for (k, sites) in data.items()]
    #[(1, 'book'), (1, 'elbow'), (2, 'hand'), (3, 'eye'), (4, 'foot'), (5, 'arm')]
    
    out, seen = {}, []
    for window in _split(partitions):
        window = sorted(window)
        label = f'p{window[0]}' if len(window) == 1 else f'p{window[0]}-{window[-1]}'
        out[label] = []
        for size in window:
            for partition in [p for (n,p) in sizes if n == size]:
                out[label].extend(data[partition])
                if dupe := len([site for site in data[partition] if site in seen]):
                    raise ValueError(f"Sites in multiple partitions: {dupe}")
                seen.extend(data[partition])
    
    missing = [s for s in range(1, max(seen)) if s not in seen]
    if len(missing):
        warn(f"Some sites are ignored: {missing}")
    
    return out
        
    
    