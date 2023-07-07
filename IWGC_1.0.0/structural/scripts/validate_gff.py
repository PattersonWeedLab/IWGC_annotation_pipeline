import os
import gffutils
from BCBio import GFF


def load_gff(gff, dbpath):
    """

    :param gff: path to gff file
    :param dbpath: path to write gffdb to
    :return: gffutils gff database.
    """

    if os.path.isfile(dbpath):
        db = gffutils.FeatureDB(dbpath, keep_order=True)
    else:
        db = gffutils.create_db(gff, dbfn=dbpath, force=True, keep_order=True,
                                merge_strategy='create_unique', sort_attribute_values=True)
    return db

def create_mini_gffdb(tmp_dir,tmp_file,cache,uniqit = None):
    """
    Notes: This assumes that the subfeatures to be updated are nearly unique in the file. That is they contain the
    unique gene ID as part of their names.
    for example:
    KSGCh10005.gene.cds
    KSGCh10005.gene.cds
    KSGCh10005.gene.cds
    would convert into
    KSGCh10005.gene.cds
    KSGCh10005.gene.cds_1
    KSGCh10005.gene.cds_2
    However,
    CDS
    CDS
    CDS
    would be unacceptable because it would only be created as a unique feature in the context of all the features within
    the cache list. Now if the Cache list contains the entire annotation, then this could be gotten away with ...
    probably

    :param tmp_dir: tempory directry
    :param tmp_file: tmp file to write gff to
    :param cache: list containing gff lines including line returns. for writing to gff
    :param uniqit: list of feature types that are to be updated with unique IDs. Includes CDS, five_prime_utr
                    and three_prime_utr
    :return: prints to stdout an updated gff with unique terms for cds.
    """
    if uniqit is None:
        uniqit = ['CDS', 'five_prime_UTR', 'three_prime_UTR']  # these have no children.
        # if features with childern are added  addtional coding will need to be done to make
        # sure the gff is not broken by mis-matched Parent and ID tags!!!
    with open(f'{tmp_dir}/{tmp_file}', 'w') as f:
        f.write(''.join(cache))
    db = load_gff(f'{tmp_dir}/{tmp_file}', f'{tmp_dir}/{tmp_file}.db')
    for gene_feat in db.all_features(featuretype='gene'):
        print(gene_feat)
        for subfeat in db.children(id=gene_feat.id):
            if subfeat.featuretype in uniqit:
                subfeat.attributes['ID'] = subfeat.id
            print(subfeat)
    os.system(f'rm {tmp_dir}/{tmp_file}.db {tmp_dir}/{tmp_file}')

def main():
        import argparse

        description = """
                   Takes as input a gff and checks that five_primer_utr, three_prime_utr, and CDS features are all unique. 
                   In the event they are not unique, it assigns a unique ID to them. Do not use to modify features whith
                   "Children", since no code is written to update 'Parent' attribute with new ID attribute. 
                   """
        import os

        # command line usage:
        cwd = os.getcwd()
        tmpdir = f'{cwd}/vg_tmp_dir'
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument('--gff',
                            help="gff to be converted to database",
                            dest='gff')
        parser.add_argument('--tmp_dir',
                            help="name for dir tmp databases to written to",
                            dest='tmp_dir',
                            default=tmpdir)
        parser.add_argument('--tmp_file',
                            help="name for tmp database",
                            dest='tmp',
                            default='tmp')
        arg = parser.parse_args()

        # make directory for tmp files if it does not exist already.

        if os.path.isdir(arg.tmp_dir) is False:
            os.system(f' mkdir {arg.tmp_dir}')

        in_file = arg.gff
        tmp_dir = arg.tmp_dir
        tmp_file = arg.tmp
        db_size = 1000
        tmp_dir = tmp_dir.rstrip('/')
        with open(in_file) as in_handle:
            counter = 0
            cache = ['##gff-version 3']
            print(''.join(cache))  # this comment line is not returned by the database so we need to print it now.
            for line in in_handle:
                sline = line.split("\t")
                if len(sline) >= 3:
                    if sline[2] == 'gene' and counter >= db_size:
                        create_mini_gffdb(tmp_dir=tmp_dir, tmp_file=tmp_file, cache=cache, uniqit=None)
                        counter = 0
                        cache = ['##gff-version 3\n']
                    elif sline[2] == 'gene':
                        counter += 1
                cache.append(line)
            """
            Get any remaining genes printed out.
            """
            create_mini_gffdb(tmp_dir=tmp_dir, tmp_file=tmp_file, cache=cache, uniqit=None)

if __name__ == '__main__':
    main()














