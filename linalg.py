import typing
import itertools
import math

#Implementation of Vector data type as they are not native to Python
class Vector():
    #Hardcoded dimensions as we are working exclusively in 2D
    def __init__(self, a:float, b:float):
        self.x:float = a
        self.y:float = b

    #Calculates Magnitude of vector.
    def Magnitude(self):
        #Calculates Magnitude
        M: float = math.sqrt((self.x)**2 + (self.y)**2)
        return M

    #Scalar multiplication
    def __mul__(self, S:float):
        x1 = self.x * S
        x2 = self.y * S
        return Vector(x1,x2)

    #Returns vector as string
    def __str__(self):
        return str(self.x) + "," + str(self.y)

    #Calculates Dot Product of two vectors
    def DotProduct(self, V):
        DP:float = self.x * V.x + self.y * V.y
        return DP

    #Piecewise addition of vectors
    def __add__(self, V):
        x1 = self.x + V.x
        y1 = self.y + V.y
        return Vector(x1,y1)

    #Calculates clockwise angle of vector
    def CalculateAngle(self):
        return math.atan2(self.y,self.x)

    def Normalise(self):
        return Vector(self.x / self.Magnitude(), self.y / self.Magnitude())

#Euclid's Algorithm for calculating the Greatest Common Divisor of two rational numbers
def EuclidAlgorithm(a,b):
    while a != b:
        if a > b:
            a = a - b
        else:
            b = b - a
    return a

class ExtendedMatrix():

    def __init__(self, Rows:int, Columns:int, Values):
        self.Numbers:list = Values
        self.Rows:int = Rows
        self.Columns:int = Columns
        self.RREF:list = []

    #Performs Gaussian Elimination
    def Gaussian(self):
        #Creates copy of original matrix
        self.RREF = self.Numbers
        #Sets pointer
        pointer:int = 0
        for r in range(self.Rows):
            if pointer >= self.Columns:
                return self.RREF
            i = r
            #Iterates through the row until a nonzero value is found
            while self.RREF[i][pointer] == 0:
                i += 1
                if i == self.Rows:
                    i = r
                    pointer += 1
                    if self.Columns == pointer:
                        
                        return self.RREF
            #Swaps rows
            self.RREF[i],self.RREF[r] = self.RREF[r],self.RREF[i]
            LeftValue = self.RREF[r][pointer]
            #Divides through to make leftmost value equal to 1
            self.RREF[r] = [ x / float(LeftValue) for x in self.RREF[r]]
            for i in range(self.Rows):
                if i != r:
                    LeftValue = self.RREF[i][pointer]
                    
                    self.RREF[i] = [ iv - LeftValue*rv for rv,iv in zip(self.RREF[r],self.RREF[i])]
            pointer += 1

    def AddRow(self, Row):
        self.Rows += 1
        self.Numbers.append(Row)
        

    #Returns the solution of a system of linear equations represented by a matrix in RREF
    def SolveFromRREF(self):
        Values:list = [0] * self.Rows
        if not self.RREF:
            return Values
        
        #Gets the value from the bottom row
        Values[-1] = self.RREF[-1][-1] / self.RREF[-1][-2]
        i = 0
        while i < self.Rows - 1:
            a = 0
            j = 0
            while j < i:
                a += self.RREF[-(i+2)][-(j+2)] * Values[-(j+1)]
                j += 1
            Values[-(i+2)] = (self.RREF[-(i+2)][-1] - a) / self.RREF[-(i+1)][-(j+2)]
            i += 1
            
        #Performs Euclid's algorithm to calculate Greatest Common Divisor
        LCM = EuclidAlgorithm(Values[0],Values[1])
        #Iterates over the list, calculating GCD of all elements
        for i in range(len(Values)-1):
            LCM = EuclidAlgorithm(LCM, Values[i+1])
        #Divides each element by the LCM to yield an integer
        Values = [i / LCM for i in Values]
        return Values

    def __str_(self):
        pass
