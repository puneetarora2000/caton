from scipy.signal import filtfilt

if __name__ == '__main__':
    usage = """
    Downsample your data by a factor of 16x
    
    %prog your_dat_file.dat [options]
    %prog -h displays help"""
    parser = OptionParser(usage)
    parser.add_options([caton.myopts.probe,caton.myopts.output])
    (opts,args) = parser.parse_args()
