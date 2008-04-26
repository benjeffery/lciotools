import math

def cosTheta(mom):
    pt2 = (mom[0] * mom[0]) + (mom[1] * mom[1])
    pt = math.sqrt(pt2)
    ptot = math.sqrt(pt2 + (mom[2] * mom[2]))
    return mom[2]/ptot
    
def xyRadius(vec):
    return math.sqrt((vec[0]*vec[0])+(vec[1]*vec[1]))
    
def threeDRadius(vec):
    return math.sqrt((vec[0]*vec[0])+(vec[1]*vec[1])+(vec[2]*vec[2]))
    
def sortByFunc(func):
    def sorter(x,y):
        if func(x)>func(y):
            return 1
        elif func(x)<func(y):
            return -1
        else: # func(x)==func(y):
            return 0
    return sorter
    
def sub(vec1,vec2):
    return (vec1[0]-vec2[0],vec1[1]-vec2[1],vec1[2]-vec2[2])
    


    

