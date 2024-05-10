#!/usr/bin/env python3
from itertools import product, combinations, chain
import sys, argparse, os
import numpy as np


def get_site_patterns(array_in,comb):
    N_alns = array_in.shape[0]
    N_taxa = array_in.shape[1]
    character_alphabet = list(map(str,np.unique(array_in)))
    character_alphabet = [item for item in character_alphabet if item != '-15']
    site_patterns = list(product(character_alphabet,repeat = comb))
    site_patterns = list(map('_'.join,site_patterns))
    print(f"\nTotal site patterns: {len(site_patterns)}")
    taxa_list = list(map(str,list(range(0,N_taxa))))
    taxa_list = list(combinations(taxa_list,comb))
    taxa_list = list(map('_'.join,taxa_list))
    print(f"\nTotal number of taxon combinations: {len(taxa_list)}")
    
    dtype =[tuple([i,"f8"]) for i in list(chain(*[['_'.join([x,y]) for y in site_patterns] for x in taxa_list]))]
    struc_array = np.zeros(N_alns, dtype=dtype)
    
    #dic_patterns = dict.fromkeys(site_patterns,np.repeat(0,N_alns))
 
    taxa_groups = [list(x) for x in list(combinations(list(range(N_taxa)),comb))]
    
    for i in range(N_alns):
        aln = array_in[i,:,:]
        aln = aln[:,(aln != -15).any(axis=0)]
        Aln_length = aln.shape[1] 
        for s in range(Aln_length):     
            site = aln[:,s]
            for group in taxa_groups:
                array_key = '_'.join(map(str,group+site[group].tolist()))
                struc_array[array_key][i]+=1/Aln_length     
    return(struc_array.view((float,len(struc_array.dtype.names))))


def main():
    parser = argparse.ArgumentParser(description='numeric2pattern conversion Ready for Keras')
    parser.add_argument( '--tr', help = "Train dataset in NUMPY",dest='TRAIN',default="TRAIN.npy")
    parser.add_argument( '--te', help = "Test dataset in NUMPY",dest='TEST',default="TEST.npy")
    parser.add_argument( '--co', help = "Taxa combinations",dest='COMB',type=int)
    
    args = parser.parse_args()
    
    print("Reading input")
    test_data = np.load(args.TEST)
    test_data_site = get_site_patterns(test_data,args.COMB)
    np.save("TEST_sitepattern_"+str(args.COMB), test_data_site)
    
    train_data = np.load(args.TRAIN)
    train_data_site = get_site_patterns(train_data,args.COMB)
    np.save("TRAIN_sitepattern_"+str(args.COMB), train_data_site)
    
     
if __name__ == "__main__":
    main()
