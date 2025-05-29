from collections import Counter
from pathlib import Path

from beastwords.main import Converter

def sitedistr(obj, glyph="█"):
    sizes = Counter()
    for p, s in obj.partitions.items():
        sizes[len(s)] += 1
    print(sizes)
    offset = min(sizes) - 5
    for i in range(min(sizes), max(sizes) + 1):
        n = sizes[i]
        bar = glyph * n
        print(f"{i}\t{n}\t{bar}")
    

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Prints a graph of the partition sizes')
    parser.add_argument("input", help='filename', type=Path) 
    parser.add_argument(
        '-p', "--partitions", dest='partitions', default=None, type=int,
        help="set partition number. If this is None use words", action='store'
    )
    args = parser.parse_args()
    
    xml = Converter.from_file(args.input)
    if args.partitions:
        xml.set_partitions(args.partitions)
    sitedistr(xml)

if __name__ == "__main__":
    main()


