
# 2020/12/07 11:22:38  zt
# Convert Fluent pre-integrated table to foam format
import numpy as np
import sys,os
import cantera as ct

def lineParser(line):
    line = line.replace('(','')
    line = line.replace(')','')
    value = float(line.split()[-1])
    return value

def readBlock(f,blockTag):
    print('Reading ',blockTag)
    data = []
    line = f.readline() # start with sth like (100 (102 nSpecies)(
    nData = int(lineParser(line))
    for i in range(nData):
        line = f.readline()
        try:
            # numbers
            temp = float(line.split()[0])
        except ValueError:
            # string
            temp = line.split()[0]
        data.append(temp)
    line = f.readline() # ))
    line = f.readline() # blank line
    return data

def readMultiDimsBlock(f,blockTag):
    print('Reading ',blockTag)

    line = f.readline() # start with sth like (100 (102 nSpecies)(
    line = line.replace('(','')
    line = line.replace(')','')
    dims = [int(temp) for temp in line.split()[2:]]

    n1 = dims[0] # Z
    n2 = dims[1] # Zvar
    n3 = dims[2] # Cvar
    n4 = dims[3] # C
    if (len(dims) == 5):
        nData = dims[4]
    else:
        nData = 1
    
    data = np.zeros([n1,n2,n3,n4,nData],dtype=float)
    
    # The data structure has been checked and the read values are 
    #   same as those from `Display PDF table`
    for i5 in range(nData):
        for i4 in range(n4):
            for i3 in range(n3):
                for i2 in range(n2):
                    for i1 in range(n1):
                        line = f.readline()
                        temp = lineParser(line)
                        data[i1,i2,i3,i4,i5] = temp
    line = f.readline() # ))
    line = f.readline() # blank line

    return data

def calcHs(Temperature,Y):
    
    gas = ct.Solution("gri30.cti")

    dims = Y.shape
    n1 = dims[0]
    n2 = dims[1]
    n3 = dims[2]
    n4 = dims[3]
    nSpecies = dims[4]-1

    hs = np.zeros([n1,n2,n3,n4,1],dtype=float)


    for i4 in range(n4):
        for i3 in range(n3):
            for i2 in range(n2):
                for i1 in range(n1):
                    T = Temperature[i1,i2,i3,i4,0]
                    Yi = Y[i1,i2,i3,i4,0:nSpecies]

                    gas.set_unnormalized_mole_fractions(Yi)
                    gas.TP = T,1e5

                    ht = gas.enthalpy_mass
                    gas.TP = 298.15,1e5
                    hf = gas.enthalpy_mass
                    hs[i1,i2,i3,i4,0] = ht-hf
                    
                    # Convert mole fraction to mass fraction
                    Y[i1,i2,i3,i4,0:nSpecies] = gas.Y

    return hs

def calcYcCsrc(species,Y):

    dims = Y.shape
    n1 = dims[0]
    n2 = dims[1]
    n3 = dims[2]
    n4 = dims[3]
    nCsrc = dims[4]

    YcCsrc = np.zeros([n1,n2,n3,n4,1],dtype=float)
    Yc = np.zeros([n1,n2,n3,n4,1],dtype=float)
    
    for n in range(nCsrc-1):
        if ( species[n] == 'CO2' or species[n] == 'CO' or species[n] == 'H2'
                or species[n] == 'H2O' ):
            Yc[:,:,:,:,0] += Y[:,:,:,:,n]

    for i4 in range(n4):
        for i3 in range(n3):
            for i2 in range(n2):
                for i1 in range(n1):
                    YcCsrc[i1,i2,i3,i4,0] = Yc[i1,i2,i3,i4,0]*Y[i1,i2,i3,i4,nCsrc-1]
    return YcCsrc

def toFoamTable(data,varName,fileName):
    
    print('Writing:', fileName)
    dims = data.shape
    n1 = dims[0] # Z
    n2 = dims[1] # Zvar
    n3 = dims[2] # Cvar
    n4 = dims[3] # C

    header = "FoamFile\n" \
           + "{\n" \
           + "    version     2.0;\n" \
           + "    format      ascii;\n" \
           + "    class       dictionary;\n" \
           + '    location    "constant";\n' \
           + "    object   "+varName+";\n" \
           + "}\n" \
           + varName+"_table\n"

    with open(fileName,'w') as f:
        f.write(header)
        f.write(str(n4)+'\n')
        f.write('(\n')
        for i4 in range(n4):
            f.write(str(n2)+'\n')
            f.write('(\n')
            for i2 in range(n2):
                f.write(str(n1)+'\n')
                f.write('(\n')
                for i1 in range(n1):
                    f.write(str(n3)+'\n')
                    f.write('(\n')
                    for i3 in range(n3):
                        f.write(str(data[i1,i2,i3,i4])+'\n')
                    f.write(')\n')
                f.write(')\n')
            f.write(')\n')
        f.write(');\n')

    return

def toFoamTableProps(Z,Zvar,C,Cvar,fileName):

    print('Writing:',fileName)

    with open(fileName,'w') as f:
        header = "FoamFile\n" \
               + "{\n" \
               + "    version     2.0;\n" \
               + "    format      ascii;\n" \
               + "    class       dictionary;\n" \
               + '    location    "constant";\n' \
               + "    object      tableProperties;\n" \
               + "}\n" 
        f.write(header)
    
        line = 'mixtureFractionDefinition "readFromTable";\n'
        f.write(line)
    
        line = 'interpolationType       linear4DInterpolation;\n'
        f.write(line)
    
        line = 'operatingPressure       1e5;\n'
        f.write(line)
    
        line = 'canteraFileName         "tables/Table";\n'
        f.write(line)
        
        n1 = len(Z)
        n2 = len(Zvar)
        n3 = len(Cvar)
        n4 = len(C) 
    
        line = 'C_param\n'+str(n4)+"\n(\n"
        f.write(line)
        for n in range(n4):
            f.write(str(C[n])+'\n')
        f.write(');\n')
        
        line = 'Zeta_param\n'+str(n2)+"\n(\n"
        f.write(line)
        for n in range(n2):
            f.write(str(Zvar[n])+'\n')
        f.write(');\n')
    
        line = 'Z_param\n'+str(n1)+"\n(\n"
        f.write(line)
        for n in range(n1):
            f.write(str(Z[n])+'\n')
        f.write(');\n')
    
        line = 'Ceta_param\n'+str(n3)+"\n(\n"
        f.write(line)
        for n in range(n3):
            f.write(str(Cvar[n])+'\n')
        f.write(');\n')
    
    return

def tableReader(file):

    with open(file,'r') as f:
        line = ''
        while ('))' not in line):
            line = f.readline()
            if ( 'pdf/max-species' in line ):
                max_species = lineParser(line)
                max_species = int(max_species)
            elif ( 'pdf/numf1point' in line ): 
                Znum = lineParser(line)
            elif ('pdf/numf2point' in line):
                Cnum = lineParser(line)
            elif ('pdf/numfvpoint' in line):
                Zvarnum = lineParser(line)
            elif ('pdf/numfv2point' in line):
                Cvarnum = lineParser(line)

        line = f.readline()
        line = f.readline() # Species input data

        species = readBlock(f,"Species name")

        dummy = readBlock(f,"Fuel boundary")
        dummy = readBlock(f,"Oxygen boundary")
        dummy = readBlock(f,"Another Oxygen boundary")

        dummy = readBlock(f,"Properties")

        Z = readBlock(f,"Mixture fraction grid")
        Zvar = readBlock(f,"Mixture fraction variance grid")
        Cvar = readBlock(f,"Progress variable variance grid")
        C = readBlock(f,"Progress Variable grid")

        T = readMultiDimsBlock(f,"Temperature data")

        rho = readMultiDimsBlock(f,"Density data")

        dummy = readMultiDimsBlock(f,"Unknown data")
        
        # Note: 
        #   1. Y is mole fraction
        #   2. Y contains nSpecies+1 data
        #   3. The last dim of Y is for OmegaC(Progress variable source term [1/s])
        Y = readMultiDimsBlock(f,'Species Mole fraction data')

        dummy = readMultiDimsBlock(f,"Unknown data")

    # Calculate sensible enthalpy and convert Mole fraction to Mass fraction
    # Also, calculate YcCsrc for the closure purpose of Cvar equation
    hs = calcHs(T,Y)
    
    YcCsrc = calcYcCsrc(species,Y)

    # Output in foam format
    toFoamTable(T[:,:,:,:,0],'T','Tables/T_table')

    toFoamTable(hs[:,:,:,:,0],'he','Tables/he_table')

    for n in range(int(max_species)):
        toFoamTable(Y[:,:,:,:,n],species[n],'Tables/'+species[n]+"_table")

    toFoamTable(Y[:,:,:,:,max_species],"Csrc","Tables/Csrc_table")
    
    toFoamTable(YcCsrc[:,:,:,:,0],"YcCsrc","Tables/YcCsrc_table")

    toFoamTableProps(Z,Zvar,C,Cvar,'tableProperties')

    print("Done writing to Foam files")
    

def main():
 
    pdfFile = sys.argv[-1]
    tableReader(pdfFile)

if __name__ == '__main__':
    main()
