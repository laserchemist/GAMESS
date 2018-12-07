# JMS 04 December 2018
#     Python script to read GAMESS MD trajectory from output file and save as xmol
#     and process data. Beta version
import re

class gamess(object): # Class to extract trajectory from GAMESS file
    """docstring forgamess."""
    def __init__(self,file):
        self.description="Class to parse GAMESS trajectory output"
        self.file=file

    def xyz(self):
        #print("In calling function",self.file)
        with open( self.file, "r" ) as source:
            coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")
        col2=[[]] # Key to creating an empty 2D List to fill with x,y,z coordinate
        index=-1
        for i in coord:
            #print(i[2])
            index+=1
            col2.append([])
            for j in i:
               col2[index].append([float(i) for i in j[2:5]]) # Convert to float element by element
        coordxyz=col2
        self.xyz=coordxyz
        return coordxyz

    def atom(self):
        #print("In calling function",self.file)
        with open( self.file, "r" ) as source:
            coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")
        col2=[[]] # Key to creating an empty 2D List to fill with x,y,z coordinate
        index=-1
        for i in coord:
            #print(i[2])
            index+=1
            col2.append([])
            for j in i:
               col2[index].append(j[0]) # Grab element label
        self.atom=col2
        return col2
    def mass(self):
        #print("In calling function",self.file)
        with open( self.file, "r" ) as source:
            coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")
            col2=[[]] # Key to creating an empty 2D List to fill with x,y,z coordinate
            index=-1
            for i in coord:
                #print(i[2])
                index+=1
                col2.append([])
                for j in i:
                    col2[index].append(float(j[1])) # Grab element label
        self.mass=col2
        return col2




# Need function to return array of COORDINATES
def coordinate(source_file,head,foot):
    buffer= []
    coord=[]
    index=-1
    atom=-1
    for line in source_file:
        if line.startswith( head ):
            index+=1
            coord.append([])
            for line in source_file:
                 if line.find( foot)>0:
                     break
                 atom+=1
                 coord[index].append(line.rstrip('\n').split())
    return coord

path='/Users/jms/local/molecules/imidazole/'
file='output65.txt'
outlines=[]
numline=[]
substr='fsec'
pattern = re.compile(r"(\+\d{1,2})?[\s.-]?\d{3}[\s.-]?\d{4}") # Compile regex
filepath=path+file
with open(filepath, "rt") as in_file: # open file given in file for reading text data
	for linenum, line in enumerate(in_file):
		if line.lower().find(substr) != -1:
			outlines.append((linenum,line.rstrip('\n')))
		if pattern.search(line) != None: # If pattern search finds a match,
                        numline.append((linenum, line.rstrip('\n')))

for linenum, line in outlines:
  print("Line ",linenum, ": ", line, sep=' ')
for linenum, line in numline:
  print("Line ",linenum, ": ", line, sep=' ')

coord=[]
step=0
atom=0
with open( filepath, "r" ) as source:
    coord=coordinate(source," QM ATOM COORDINATES","CARTESIAN COORDINATES (ANG)")

print("COORDINATES:")
calc1=gamess(filepath)
print(calc1.description)
calc2=gamess(filepath)
print(calc1.file)
print("Testing class\n||||||||||||||||||")
coords=calc1.xyz()
print(coords[0])
print("atoms: ",calc1.atom())
print("mass: ",calc1.mass())
print("Length:",len(coords[0]))
print("Length overall:",len(coords))
# Write xmol file
natoms=len(coords[0])
steps=len(coords)
mass=calc1.mass
filename=filepath+".xmol"
print("Writing file to:",filename)
with open(filename, 'w') as f:
    for item in coords:
        f.write(str(natoms)+"\n")
        f.write("Output from GAMESS\n")
        atom=0
        for xyz in item:
            strxyz=" ".join([str(i) for i in xyz])
            textline=str(int(mass[0][atom]))+" "+strxyz
            print(textline)
            f.write(textline)
            f.write("\n")
            atom+=1
