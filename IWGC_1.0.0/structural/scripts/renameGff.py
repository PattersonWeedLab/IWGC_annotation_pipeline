from gff3 import Gff3

def clean_gff_line(line,new_id,old_id):
    seqid = line['seqid']
    source = line['source']
    feattype = line['type']
    start = line['start']
    end = line['end']
    score = line['score']
    strand = line['strand']
    phase = line['phase']
    id = line['attributes']['ID'].replace(old_id,new_id)
    attrib = f'ID={id};'
    raw_line = line['line_raw']
    if 'Name' in line['attributes']:
        name = id
        attrib += f'Name={name};'
    if 'Parent' in line['attributes']:
        assert len(line['attributes']['Parent']) == 1, f'unpredictable parent in line \n{raw_line}'
        parent = line['attributes']['Parent'][0].replace(old_id,new_id)
        attrib += f'Parent={parent};'

    attrib = attrib.rstrip(';')
    gffline = [seqid, source, feattype,str(start),str(end),str(score),str(strand),str(phase),attrib]
    return gffline

def rename_renunmber(gff_path,org_code):
    gff = Gff3()
    gff.parse(gff_path)
    counter = 0
    org = org_code
    geneid = None
    for line in gff.lines:
        if len(line['parents']) == 0 and line['type'] == 'gene':
            counter += 10
            seqid = line['seqid']
            genenum = '{:06d}'.format(counter)
            geneid = f'{org}{seqid}g{genenum}'
            old_geneid = line['attributes']['ID']
        if geneid is not None and '#' not in line['line_raw'][0:3]:
            print('\t'.join(clean_gff_line(line=line, new_id=geneid, old_id=old_geneid)))
        else:
            print(line['line_raw'].rstrip())


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-g,--gff',
                        help='gff produced by exonerate',
                        dest='gff_path',
                        required=True,
                        type=str
                        )
    parser.add_argument('-t,--tag,--org_code',
                        dest='org_code',
                        required=True,
                        help='Code that designates organism for the genome. Optimally a five letter code with the '
                             'first 2 indicating genus and the last three species for '
                             'example, Conyza sumatrensis -> Cosum')
    args = parser.parse_args()

    rename_renunmber(gff_path=args.gff_path, org_code=args.org_code)

