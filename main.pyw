#dependencies
import time
import datetime
import typing
from linalg import *
import random
import math
import random
import typing
import os
#hides pygame support prompt
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"
import pygame
from pygame.locals import *

#Named colour constants
BLACK = (0,0,0)
WHITE = (255,255,255)
RED = (192,0,0)
GREEN = (0,192,0)
BLUE = (0,0,192)
YELLOW = (192,192,0)
PURPLE = (192,0,192)
TASTEFULBLUE = (16,48,64)
XMaximum = 800
YMaximum = 600

Colours = {0:WHITE, 1:RED, 2: GREEN, 3:BLUE,4:YELLOW, 5:PURPLE}
#Boltzmann's constant, scaled up
Boltzmann = 1.38




#Samples standard normal distribution for Boltzmann's Dist
def BoxMuller():
    a,b = random.uniform(0,1), random.uniform(0,1)
    c = (a,b) if a > b else (b,a)
    d = math.sqrt(-2 * math.log(c[0]))

    X = d * math.cos(2*math.pi*c[1])
    Y = d * math.sin(2*math.pi*c[1])

    return X,Y

#Samples speed from Boltzmann's Distribution
def BoltzmannDist(Mass):

    #Samples standard normal distribution
    X,Y = BoxMuller()
    Const = math.sqrt(Boltzmann / Mass)
    Speeds = Vector(X,Y)
    #returns component Velocities
    return Speeds * Const



#Moleule class
class Molecule():
    
    def __init__(self,Pos,Bound, Type, Mass):
        #Sets position randomly
        self.Position = Pos if Pos else Vector(random.randint(40 + int(Bound),XMaximum-40 - int(Bound)), random.randint(40 + int(Bound),YMaximum-40 - int(Bound)))
        #Sets velocity as a normalised vector.
        self.Radius = 8
        #Number of frames which the particle can't react for
        #(prevents newly created particles reacting with each other)
        self.CreationFlag = 15
        self.Type = Type
        self.Mass = Mass
        #Assigns colour based off of type
        self.Colour = Colours[self.Type]
        #Gets random velocities
        a = BoltzmannDist(Mass)
        self.Velocity = a.Normalise()
        self.Speed = max(a.Magnitude(),0.03)

    #Calculates position of molecule from the current position + velocity
    def CalculateMotion(self,Bound,Temperature):
        #Checks if particle is touching edge of box
        self.Collision(Bound)
        self.Position = self.Position + self.Velocity * self.Speed * Temperature
        #Moves particle inside box if it ends up outside
        self.Position = Vector(random.randint(40 + int(Bound),XMaximum-40 - int(Bound)), random.randint(40 + int(Bound),YMaximum-40 - int(Bound))) if self.Position.x > XMaximum else self.Position

    #Draws molecule on screen
    def Draw(self):
        GlobalPG.DrawAtom(self.Colour,(self.Position.x,self.Position.y),self.Radius)

        #Changes direction of particle if the particle is touching the edge of the screen
    def Collision(self,Bound):
        self.Velocity.x *= 1 if (self.Position.x - self.Radius > Bound and self.Position.x + self.Radius < XMaximum - Bound) else -1
        self.Velocity.y *= 1 if (self.Position.y - self.Radius> Bound and self.Position.y + self.Radius< YMaximum - Bound) else -1

#Reaction class
class Reaction():

    def __init__(self,S,E,A,P,M):

        self.Particles = S
        self.Enthalpy = E
        self.Activation = A
        self.Products = P
        self.SubstancesCount = len(self.Particles)
        self.Masses = M


#Pygame button class
class PGButton():

    def __init__(self,position,colour,fun,args,dim,text,type):
        self.Position:tuple = position
        self.Colour:tuple = colour
        self.Function:str = fun
        self.Dimensions:tuple = dim
        self.Text:str = text
        self.Args = args
        self.Type = type
        
    #Draws rectangle equivalent to button
    def Draw(self,font,surface):
            A = pygame.Rect(self.Position[0], self.Position[1], self.Dimensions[0], self.Dimensions[1])
            pygame.draw.rect(surface,WHITE,A)
            text_surface = font.render(self.Text, False, self.Colour)
            surface.blit(text_surface, (self.Position[0] + 100, self.Position[1]))
            
            #Draws minus sign or plus sign depending on Type attribute
            if self.Type > 0:
                B = pygame.Rect(self.Position[0]+5, self.Position[1] + self.Dimensions[1]/2 - 5, self.Dimensions[0] - 10, 10)
                pygame.draw.rect(surface, BLACK,B)
                if self.Type > 1:
                    C = pygame.Rect(self.Position[0] + self.Dimensions[0]/2 - 5, self.Position[1] +5, 10, self.Dimensions[1] - 10)
                    pygame.draw.rect(surface, BLACK,C)

    #Checks if coordinates are within the button and returns True if so
    def CheckInButton(self,ClickPos):
        if ClickPos[0] > self.Position[0] and ClickPos[0] < self.Position[0] + self.Dimensions[0]:
            if ClickPos[1] > self.Position[1] and ClickPos[1] < self.Position[1] + self.Dimensions[1]:
                return True
        return False

class Simulation():
    
    def __init__(self):

        #Initialises variables
        
        self.Reac = Reaction(["CH3COOH","CH3CH2OH","CH3COOCH2CH3","H2O"], 0.6, 6, [[0,1],[2,3]], [32, 26, 48, 10])
        self.particlecount = [0] * self.Reac.SubstancesCount
        self.InitTime = time.time()
        
        #Opens file for outputting data
        self.Date = str(datetime.datetime.utcfromtimestamp(self.InitTime).strftime('%d-%m-%Y-%H%M-%S'))
        #File is automatically created
        try:
            self.OutputFile = open("data/simulation/" + self.Date + ".txt", "w")
        #Error handling for if simulation folder does not exist
        except FileNotFoundError:
            self.OutputFile = open(self.Date + ".txt", "w")

        #Initialises pygame
        pygame.init()
        self.screen = pygame.display.set_mode((XMaximum + 480, YMaximum))
        pygame.display.set_caption("Chemical Reaction Simulation")
        try:
            icon = pygame.image.load("assets/icon.png").convert()
            icon = pygame.transform.scale(icon, (32,32))
            pygame.display.set_icon(icon)
        except FileNotFoundError:
            pass
        
        self.C = pygame.time.Clock()
        pygame.font.init()
        self.displayfont = pygame.font.SysFont('Calibri', 20)

        #Sets size of buttons as a function of the dimensions of the window
        #So the UI is scaled
        ButtonDim = (XMaximum // 24, XMaximum // 24)
        ButtonSpacing = XMaximum // 16

        #Instantiates Button objects
        self.Buttons = []
        for i in range(self.Reac.SubstancesCount):
            self.Buttons.append(PGButton((XMaximum+80, 24 + ButtonSpacing * (4+i)),Colours[i],"AddSubstance", i,ButtonDim,"",2))
            self.Buttons.append(PGButton((XMaximum+160, 24 +  ButtonSpacing * (4+i)),Colours[i],"RemoveSubstance", i,ButtonDim,"Add "+self.Reac.Particles[i],1))

        self.Buttons.append(PGButton((XMaximum+80, 24),WHITE,"TogglePause",1,ButtonDim,"",0))
        self.Buttons.append(PGButton((XMaximum+160 ,24),WHITE,"Reset",24,ButtonDim,"Pause/Reset",0))
        self.Buttons.append(PGButton((XMaximum+80, 24 + ButtonSpacing),WHITE,"IncreasePressure",0.1,ButtonDim,"",2))
        self.Buttons.append(PGButton((XMaximum+160, 24+ ButtonSpacing),WHITE,"DecreasePressure",0.1,ButtonDim,"Pressure",1))
        self.Buttons.append(PGButton((XMaximum+80, 24+ 2*ButtonSpacing),WHITE,"IncreaseTemp",1,ButtonDim,"",2))
        self.Buttons.append(PGButton((XMaximum+160, 24+ 2*ButtonSpacing),WHITE,"DecreaseTemp",1,ButtonDim,"Temperature",1))
        self.Buttons.append(PGButton((XMaximum+160, 24+ 3*ButtonSpacing),WHITE,"Clear",1,ButtonDim,"Clear",1))

        self.Reset(24)

    #Increases pressure
    def IncreasePressure(self,n):
        b = (self.Pressure-1) * self.PressureConst
        self.Pressure = self.Pressure + n if (self.Pressure-1+n) * self.PressureConst <= YMaximum / 2 - 50 else self.Pressure

        #Repositions molecules outside the new simulation area
        for i in self.Objects:
            if i.Position.x < (self.Pressure-1) * self.PressureConst:
                i.Position.x += self.PressureConst
            elif i.Position.x > XMaximum - (self.Pressure-1) * self.PressureConst:
                i.Position.x -= self.PressureConst
            if i.Position.y < (self.Pressure-1) * self.PressureConst:
                i.Position.y += self.PressureConst
            elif i.Position.y > YMaximum - (self.Pressure-1) * self.PressureConst:
                i.Position.y -= self.PressureConst

    #Decreases pressure
    def DecreasePressure(self,n):
        self.Pressure = self.Pressure - n if self.Pressure > 1 else self.Pressure

    #Increases pressure
    def IncreaseTemp(self,n):
        self.Temperature+= 1

    #Decreases pressure
    def DecreaseTemp(self,n):
        self.Temperature -= 1 if self.Temperature > n else 0
        

    #Adds specified substance to simulation
    def AddSubstance(self, Type):
        #Adds some molecules
        for i in range(3):
            self.Objects.append(Molecule((),(self.Pressure-1) * self.PressureConst, Type,self.Reac.Masses[Type]))

    #Removes specified substance from simulation
    def RemoveSubstance(self, Type):
        count = 0
        for i in self.Objects:
            if i.Type == Type:
                self.Objects.remove(i)
                count += 1
            if count > 2:
                break

    #Toggles the Pause flag, which determines whether the simulation runs or not
    def TogglePause(self,n):
        self.Pause = 0 if self.Pause else 1

    #Resets simulation to initial status
    def Reset(self,n):
        self.Pause = 1
        self.Objects = []

        self.displaycount = 0
        self.displaycycle = 0
        self.Temperature = 20
        self.Pressure = 1.0
        self.count = 0
        for i in range(n):
            for j in range(self.Reac.SubstancesCount):
                self.Objects.append(Molecule(0,0, j, self.Reac.Masses[j]))

    #Removes all particles from screen
    def Clear(self,n):
        self.Objects = []
        
        
    #Draws a circle representing an atom
    def DrawAtom(self,Colour,Centre,Radius):
        pygame.draw.circle(self.screen,color=Colour,center=Centre,radius=Radius)

    #Checks if two particles have collided
    def CollisionCheck(self,a,b):
        #Gets distance between centres of particles as vector
        Distance = Vector(b.Position.x - a.Position.x, b.Position.y - a.Position.y)
        #Calculates magnitude of distance
        Magn = Distance.Magnitude()
        if Magn <= a.Radius + b.Radius:
            return True
        return False

    #Calculates the velocity of the particles formed when two particles collide
    def CalculateVelocityAfterCollision(self,Particles):
        
        #Calculates angle that each particle travels at
        Angles = []
        Velocities = []
        for i in Particles:
            Angles.append(i.Velocity.CalculateAngle())
        #Sets the angle of each new particle as the average of the two particles before but with variance of pi/3
        NewAngle = sum(Angles) / len(Angles)
        Diff = 2 * math.pi / len(Particles) + 1
        for i in range(len(Particles)):
            Velocities.append(Vector(math.cos(NewAngle + 2*math.pi/Diff * i / 3), math.sin(NewAngle + Diff * i / 3)))
        return Velocities

    #Routine if two particles collide
    def OnCollision(self,Particles):
        e = self.Reac.Enthalpy if Particles[0].Type < len(self.Reac.Products[0]) else 1/self.Reac.Enthalpy
        #Calculates velocity of new particles
        Velocities = self.CalculateVelocityAfterCollision(Particles)
        #Generates new particles

        Molecules = []

        for i in range(len(Particles)):

            a = (Molecule(Particles[i].Position,(self.Pressure-1) * self.PressureConst,self.Reac.SubstancesCount - Particles[i].Type - 1, self.Reac.Masses[Particles[i].Type - 1]))
            a.Velocity = Velocities[i]
            a.CalculateMotion((self.Pressure-1)*self.PressureConst,self.Temperature / 20)
            Molecules.append(a)
        
        return Molecules


    #Renders UI
    def RenderButtonsAndText(self):
        #Draws Buttons
        for i in self.Buttons:
            i.Draw(self.displayfont, self.screen)


            
        #Renders the number of reactions
        text_surface = self.displayfont.render("Temperature: " + str(self.Temperature) + "C, Pressure: " + str(round(self.Pressure,1)) + "atm" , False, WHITE)
        self.screen.blit(text_surface, (XMaximum + 80, YMaximum - 80))
        text_surface = self.displayfont.render("No. of reactions in last second: " + str(self.displaycount) , False, WHITE)
        self.screen.blit(text_surface, (XMaximum + 80, YMaximum - 40))

        #Renders the number of each species 
        for i in range(self.Reac.SubstancesCount):
                text_surface = self.displayfont.render(str(self.particlecount[i]) , False, Colours[i])
                self.screen.blit(text_surface, (XMaximum + 80 + 40 * i, YMaximum - 120))
        
        #Draws rectangles to represent area that the reaction takes place in
        pygame.draw.rect(self.screen,BLACK,(0,0,XMaximum, (self.Pressure-1)*self.PressureConst))
        pygame.draw.rect(self.screen,BLACK,(0,YMaximum-(self.Pressure-1) * self.PressureConst,XMaximum,self.Pressure*self.PressureConst))
        pygame.draw.rect(self.screen,BLACK,(0,(self.Pressure-1) * self.PressureConst,(self.Pressure-1) * self.PressureConst,YMaximum - 2*(self.Pressure-1)*self.PressureConst))
        pygame.draw.rect(self.screen,BLACK,(XMaximum-(self.Pressure-1) * self.PressureConst,(self.Pressure-1) * self.PressureConst,(self.Pressure-1) * self.PressureConst,YMaximum - 2*(self.Pressure-1)*self.PressureConst))

            



        

    #Main simulation loop
    def MainLoop(self):

        #Clears Screen of stuff
        self.screen.fill(TASTEFULBLUE)
        #Lists of particles to be added
        self.ToAdd = []
        self.ToRemove = []
        if not self.Pause:
            #Iterates through each particle
            for i in range(len(self.Objects)-1):
                flag = False
                for j in range(i+1,len(self.Objects)):
                    #Checks if particles have collided
                    a = self.CollisionCheck(self.Objects[i],self.Objects[j])
                    #Collision does not count if a particle has already collided or has recently been created
                    if (a and (not flag) and (not i==j) and (not i in self.ToRemove) and (not j in self.ToRemove) and not(self.Objects[i].CreationFlag) and not(self.Objects[j].CreationFlag)):
                        #Checks if the two particles can react with each other
                        e = [self.Objects[i].Type, self.Objects[j].Type]
                        if (e in self.Reac.Products or e[::-1] in self.Reac.Products):
                            #Checks if the activation energy is present
                            if (self.Objects[i].Speed + self.Objects[j].Speed) * self.Temperature > self.Reac.Activation:
                                NewMolecules = self.OnCollision([self.Objects[i], self.Objects[j]])
                                self.ToAdd += NewMolecules
                                #Sets flag so particles cannot react again
                                flag = True
                                self.ToRemove += [i,j]
                                #Increments counter for number of reactions
                                self.count += 1
                                
                #Calculates position of particle if particle does not react
                if not flag:
                    self.Objects[i].CalculateMotion((self.Pressure-1)*self.PressureConst,self.Temperature)
                

            for i in range(len(self.Objects)):
                if not i in self.ToRemove:
                    self.ToAdd.append(self.Objects[i])
                #Decrements crration 
                if self.Objects[i].CreationFlag:
                    self.Objects[i].CreationFlag -= 1
            self.Objects = self.ToAdd

        #Updates frame counter
        self.displaycycle += self.C.tick(30)
        #Routine every 1 second
        if self.displaycycle > 1000:
            if not self.Pause:
                self.displaycount = self.count
                self.displaycycle = 0
                
                #Writes data to file
                self.OutputFile.write(str(self.particlecount)[1:-1]+ ", " + str(self.count) + "\n")
                self.count = 0
                self.particlecount = [0] * self.Reac.SubstancesCount

                #Counts number of particles
                for i in self.Objects:
                    self.particlecount[i.Type] += 1

        #Draws all objects
        for i in self.Objects:
            i.Draw()

        #Renders buttons and text
        self.RenderButtonsAndText()
        #Updates display
        pygame.display.flip()
            


        #Event handling
        #Closes the program if the quit button or esc key is pressed.
        for event in pygame.event.get():
           if event.type == QUIT:
               self.OutputFile.close()
               pygame.quit()
               sys.exit()
              
           #If ESC key pressed
           elif event.type == KEYDOWN:
              if event.key == K_ESCAPE:
                  self.OutputFile.close()
                  pygame.quit()
                  sys.exit()
                 
           #On click
           elif event.type == MOUSEBUTTONDOWN:
                 pos = pygame.mouse.get_pos()
                 for i in self.Buttons:
                    if i.CheckInButton(pos):
                        #Executes function associated with button
                        a = getattr(self, i.Function)
                        a(i.Args)

            
                
        
#Initialises simulation
GlobalPG = Simulation()

#Execute main loop of simulation
while True:
    GlobalPG.MainLoop()

