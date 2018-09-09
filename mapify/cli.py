import argparse

import products


def validprods():
    prods = ['all']
    prods.extend(list(products.prodmap().keys()))
    return prods


def main(args: argparse.Namespace) -> None:
    print(args)
    pass


if __name__ == '__main__':
    desc = '''
    Map generation tool for CCDC test results.
    This version is currently orientated towards "ncompare" CCD and the annualized classification.
    *Not to be used with production results* They will use a different formatting schema.
    Dates are given in ISO-8601 which is YYYY-MM-DD.
    Valid products are: {prods}
    '''.format(prods=validprods())

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('jdir',
                        type=str,
                        help='Input directory of JSON files to process')
    parser.add_argument('pdir',
                        type=str,
                        help='Input directory of pickle files to process')
    parser.add_argument('outdir',
                        type=str,
                        help='Output directory for GeoTiff maps')
    parser.add_argument('dates',
                        type=str,
                        help='Comma separated list of dates in ISO-8601 to build the requested products for')
    parser.add_argument('prods',
                        type=str,
                        help='Comma separated list of products to build')
    parser.add_argument('--no-fill-begin',
                        dest='fill_begin',
                        action='store_false',
                        default=True,
                        help='Turn off the LC filling technique for the beginning of a timeseries')
    parser.add_argument('--no-fill-end',
                        dest='fill_end',
                        action='store_false',
                        default=True,
                        help='Turn off the LC filling technique for the end of a timeseries')
    parser.add_argument('--no-fill-samelc',
                        dest='fill_samelc',
                        action='store_false',
                        default=True,
                        help='Turn off the LC filling technique between segments that are the same LC')
    parser.add_argument('--no-fill-difflc',
                        dest='fill_difflc',
                        action='store_false',
                        default=True,
                        help='Turn off the LC filling technique between segments that are different LC')
    parser.add_argument('--no-fill-nomodel',
                        dest='fill_nomodel',
                        action='store_false',
                        default=True,
                        help='Turn off the LC filling technique where there isn\'t a model')
    parser.add_argument('--cpu',
                        dest='cpu',
                        action='store',
                        type=int,
                        default=1,
                        help='Number processes to spawn')
    args = parser.parse_args()

    main(args)