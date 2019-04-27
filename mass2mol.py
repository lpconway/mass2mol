# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 21:53:09 2019

@author: Louis
"""
import numpy as np

class mass2mol:
    """
    A Python implementation of the algorithm developed by Sebastian Bocker and Zsuzsanna Liptak (A Fast and Simple Algorithm for the Money Changing Problem, Algorithmica, 2007, 48, 413–432)
    The mass2mol object generated the extended residue table for CHNOPS with a blowup factor of 10^5 upon initialisation.
    The combinations of CHNOPS atoms which have the same monoisotopic masses as a given mass can then be calculated using the find_formula method.
    """
    
    def __gcd(self, a, b):
        """Calculate the Greatest Common Divisor of a and b.
    
        Unless b==0, the result will have the same sign as b (so that when
        b is divided by it, the result comes out positive).
        """
        while b:
            a, b = b, a%b
        return a
    
    def __lcm(self, a, b):
        """
        Calculate the least common multiple of a and b.
        """
        return abs(a * b)//self.__gcd(a, b)
    
    def __find_all(self, M, i, c, output, ppm):
        """
        Find all decompositions using the extended residue table
        """
        if i == 0:
            c[0] = M//self.a[0]
            output.append((c[:], ppm))
            return
        LCM = self.__lcm(self.a[0], self.a[i])
        l = LCM//self.a[i]
        for j in range(l):
            c[i] = j
            m = M - j * self.a[i]
            r = m % self.a[0]
            lbound = self.ERT[r, i-1]
            while (m >= lbound) and (lbound != -1) :
                self.__find_all(m, i-1, c, output, ppm)
                m = m - LCM
                c[i] = c[i] + l
        return
    
    def __gen_ERT(self):
        """
        Generate extended residue table
        """
        k = len(self.a)
            
        self.ERT = -np.ones((self.a[0], k), dtype='int64') # Used -1 in place of infinities: requires more computational steps to check this. Better to use unsigned int - 1?
        self.ERT[0, ] = 0
        
        for i in range(1, k):
            self.ERT[:,i] = self.ERT[:, i-1]
            d = self.__gcd(self.a[0], self.a[i])
            for p in range(d):
                n = [self.ERT[q, i-1] for q in range(self.a[0]) if (q % d == p and self.ERT[q, i-1] != -1)]
                if n == []:
                    n = -1
                else:
                    n = min(n)
                    for x in range(1, self.a[0]//d):
                        n = n + self.a[i]
                        r = n % self.a[0]
                        if self.ERT[r, i-1] != -1:
                            if n < self.ERT[r, i-1]:
                                self.ERT[r, i] = n
                            else:
                                n = self.ERT[r, i-1]
                        else:
                            self.ERT[r, i] = n

    def __init__(self, ppm=5):
        """
        Initialise mass2mol object with masses of CHNOPS and blowup factor.
        Generate extended residue table for these masses.
        """
        self.a = [100783, 1200000, 1400307, 1599492, 3097376, 3197207]
        self.blowup = 100000
        self.alphabet = 'HCNOPS'
        self.storage = []
        self.ppm = ppm
        self.__gen_ERT()
    
    def __formula(self, x):
        """
        Function to generate chemical formula string from list of numbers of atoms.
        """
        output = ''
        if x[1] != 0:
            output = output + 'C' + str(x[1])
        if x[0]!= 0:
            output = output + 'H' + str(x[0])
        for i in range(2,len(x)):
            if x[i] != 0:
                output = output + self.alphabet[i] + str(x[i])
        return output
                
    
    def find_formula(self, mass):
        """
        Function to calculate all possible combinations of atoms which have a combined mass equal to the given mass, within a certain error range.
        Returns a list of tuples consisting of the chemical formula followed by the difference from the given mass in parts per million.
        """
        mass = int(mass * self.blowup)
        error = int(mass/1000000 * self.ppm)
        output = []
        for m in range(mass-error, mass + error +1):
            ppm = (m - mass)/mass * 1000000
            self.__find_all(m, len(self.a)-1, [0]*len(self.a), output, ppm)
        formula_list = []
        for x in output:
            formula_list.append((self.__formula(x[0]), x[1]))
        return formula_list