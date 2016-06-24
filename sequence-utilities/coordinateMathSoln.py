#!/usr/bin/env python3
# Name: Louis Chatta(lchatta)
# Group Members: none

import math
TO_DEGREES = 57.2958
class Triad :
    """
 
    Author: David Bernick
    Date: March 21, 2013
    This class calculates angles and distances among a triad of points.
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad.  
            p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,0,0) ).
        """
        self.p = p
        self.q = q
        self.r = r
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b))) 
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r)) 
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) /    math.sqrt(self.d2(self.q,self.p)*self.d2(self.r,self.p))) * TO_DEGREES
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /  math.sqrt(self.d2(self.p,self.q)*self.d2(self.r,self.q))) * TO_DEGREES
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /  math.sqrt(self.d2(self.p,self.r)*self.d2(self.q,self.r))) * TO_DEGREES

"""
This program takes 3 sets of coordinates with each 
element coming before each coordinate. Formatted in 
the following format: C = (39.447, 94.657, 11.824) 
N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465) .
Based on this format, I use code mentioned in class to 
seperate the elements and the values into an array. Based 
on the seperation, the element strings will be placed at 
0 , 4 , and 8. From there I save the element string by 
removing all whitespace and stripping the '=' character.
I then map each of the values to each element string
in a dictionary. Then I load the corresponding values
into a Triad object. I then print the dPQ, dQR, the 
angleQ .
"""
coordinates = input("Enter 3 sets ot coordinates : ")

# format input to be a split list of elements and coordinates
t1 = coordinates.replace ('(', ',') # change ( to ,
t2 = t1.replace (')', ',') # change ) to ,
l = t2.split (',') # split on ,

# elements at positions 0, 4, and 8 are going to 
# element strings

element1 = "".join(l[0].split()).strip('=')
element2 = "".join(l[4].split()).strip('=')
element3 = "".join(l[8].split()).strip('=')

# map each element with its values
elementDict = {element1:(float(l[1]),float(l[2]),float(l[3])), element2:(float(l[5]),float(l[6]),float(l[7])),element3:(float(l[9]),float(l[10]),float(l[11]))}

# create the triad with each respecting element
triad  = Triad(elementDict[element1],elementDict[element2],elementDict[element3])

# print the corresponding lengths
# dPQ is the length between the first two elements
print('%s-%s bond length = %.2f'%(element1,element2,triad.dPQ()))
# dQR is the length between the last two elements
print('%s-%s bond length = %.2f'%(element2,element3,triad.dQR()))
# angleQ gives us the total angle
print('%s-%s-%s bond angle = %.2f'%(element1,element2,element3,triad.angleQ()))


