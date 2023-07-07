from gff3 import Gff3
def print_3_deep(feat):
    for sf in feat:
        if len(sf['parents']) == 0:
            print(sf['line_raw'].rstrip())
            for ssf in sf['children']:
                print(ssf['line_raw'].rstrip())
                if 'children' in ssf:
                    for sssf in ssf['children']:
                        print(sssf['line_raw'].rstrip())
                        if 'children' in sssf:
                            for ssssf in sssf['children']:
                                print(ssssf['line_raw'].rstrip())


if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-g,--gff',
                        help='gff produced by exonerate',
                        dest='gff_path',
                        required=True,
                        type=str
                        )
    parser.add_argument('-e,--exclude_list',
                        help='list of feature gene IDs to exclude',
                        dest='exclude_list',
                        type=str)
    args = parser.parse_args()

    """
    For some local testing
    class args:
        pass
    args.exclude_list = '/Users/urthstripe/OneDrive - Michigan State University/Conyza_sumatrensis/exclude.list'
    args.gff_path = '/Users/urthstripe/OneDrive - Michigan State University/Conyza_sumatrensis/Cosum_v00_maker.gff'
    """
    
    with open(args.exclude_list) as f:
        elist = set(f.read().strip().splitlines())
    gff = Gff3()
    gff.parse(args.gff_path)
    counter = 0
    y = None
    for feat in gff.features:
        if feat not in elist :
            y = gff.features[feat]
            print_3_deep(feat=y)
            counter += 1




