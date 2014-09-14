#!/usr/bin/env python

import pdb
import string
import sys

class ConvertToDict:
    def __init__(self):
        self.tmp_dict = {}
        self.return_dict = {}
    def walk_string(self, some_string):
        '''walk given text string and return a dictionary. 
        Maintain state in instance attributes in case we hit an exception'''
        l = string.split(some_string)
        for i in range(len(l)):
            key = str(i)
            self.tmp_dict[key] = int(l[i])
        return_dict = self.tmp_dict
        self.return_dict = self.tmp_dict
        self.reset()
        return return_dict
    def reset(self):
        '''clean up'''
        self.tmp_dict = {}
        self.return_dict = {}
    def get_number_dict(self, some_string):
        '''do super duper exception handling here'''
        try:
            return self.walk_string(some_string)
        except:
			#if we hit an exception, we can rely on tmp_dict 
			being a backup to the point of the exception
            return self.tmp_dict

def main():
    ctd = ConvertToDict()
    for line in file(sys.argv[1]):
        line = line.strip()
        print "*" * 40
        print "line>>", line
        print ctd.get_number_dict(line)
        print "*" * 40
    
if __name__ == "__main__":
    pdb.runcall(main)
    #main()