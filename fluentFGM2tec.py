#!/usr/bin/python

# 2019/01/24 22:48:48  zt
# convert FGM laminar flamelets to Tecplot plt

Z = []
data = []
var_name = []
#var_name.append("Z\n")
n = 0
n_grid = 0
def read_data_block(f):
    global Z, data, var_name, n, n_grid

    block = []

    line = f.readline() # header
    if(line == ""):
        # it's the end of file
        return False
    else:
        line = f.readline() # premix_chi
        
        line = f.readline() # Z
        z_temp = float(line.split()[1])
        Z.append(z_temp)
        
        line = f.readline() # species
        n_species = int(line.split()[1])
    
        line = f.readline() # gridpoints
        n_grid = int(line.split()[1])
    
        line = f.readline() # pressure
        line = f.readline() # body
        
        n_var = 2+n_species+1 # Reaction_progress+Temperature+N_species+Reaction_progress_sourceterm

        for n_temp in range(n_var):
            temp_block = []
            line = f.readline() # variable header
            if ( n == 0):
                var_name.append(line)
            for i in range(n_grid/5+1):
                line = f.readline()
                for temp in line.split():
                    temp_block.append(float(temp))
            block.append(temp_block)
    
        data.append(block)

        f.readline() # blank line
        n += 1

        return True

def readFlamelet(file):

    f = open(file,'r')
    
    flag = True
    while( flag ):
        flag = read_data_block(f )
    f.close()
    
def output_tec(file):
    f = open(file, 'w')

    header = 'Variables = '
    for var in var_name:
        header = header + '\"' +var[:len(var)-1]+'\"' +','
    f.write(header[:len(header)-1]+'\n')
    
    for i in range(len(Z)):
        header = "Zone T = \"Z=" +str(Z[i]) +"\" \n"
        print header
        f.write(header)
        for j in range(n_grid):
            line = ""
            for k in range(len(var_name)):
                line += str(data[i][k][j]) + '\t'
            f.write(line+'\n')
    f.close()
    
    

readFlamelet("fluent-fgm.fla")

output_tec("FGM-flamelet.plt")
