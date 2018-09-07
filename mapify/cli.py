import argparse

def main(args: argparse.Namespace) -> None:
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Map generation tool for CCDC test results.')
    parser.add_argument('jdir',
                        type=str,
                        help='Input directory of JSON files to process.')
    parser.add_argument('pdir',
                        type=str,
                        help='Input directory of pickle files to process.')
    parser.add_argument('outdir',
                        type=str,
                        help='Output directory for GeoTif maps.')
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
                        default=1,
                        help='Number CPU\'s to use')
    args = parser.parse_args()

    main(args)