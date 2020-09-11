import numpy as np
import os
import matplotlib.pyplot as plt

import dassflow2d
print("==========================================")
print("	  Import dassflow2d successful!		 ")
print("==========================================")

from dassflow2d import m_mesh as m_mesh
from dassflow2d import m_common as m_common
from dassflow2d import m_model as m_model
from dassflow2d import m_linear_solver as m_linear_solver
from dassflow2d import call_model as call_model
from dassflow2d import m_numeric as m_numeric

def test_wrappers():
    #============================================
    #		SIMPLE TESTS
    #============================================
    typlim = m_mesh.NodeTypeLim(23,5)
    print(typlim)
    typ = m_mesh.NodeType(2,2)
    print(typ)
    print(typ.coord)
    typ.coord.x = -5
    typ.coord.y = 7
    typ.cell[0] = 1
    typ.cell[1] = 2
    typ.edge[0] = 3
    typ.edge[1] = 4	
    print(typ)
    print(typ.coord)
    print("==========================================")
    print("	SIMPLE TESTS SUCCESSFUL!	  	     ")
    print("==========================================")

    #============================================
    #		SIMULATION PART
    #============================================
    
    ##### Read input ######
    # =====================

    filename = "input.txt"
    filepython = "input_for_python.txt"

    # Create a specific data file to read data from Python
    outF = open(filepython,"w")
    with open(filename, 'r') as filehandle:
        for line in filehandle:
            write_value = False
            if (line[0] not in ['!','=','&','/']) and ('=' in line):
                for index in range(len(line)):
                    if line[index] == ',':
                        break
                    if line[index] not in ['\t',' '] and write_value:
                        if line[index] is 'd' and line[index-1] is '.':
                            outF.write('e')  
                        else:
                            outF.write(line[index])
                    if line[index] == '=':
                        write_value = True
                outF.write("\n")
    outF.close()

    # Convert variable type from the specific data file
    fname = open(filepython,"r")
    data = fname.readlines()
    data = data[:-1]
    fname.close()
    input_data = []
    for i in range(len(data)):
        try:
            input_data += [int(data[i])]
        except ValueError:
            try:
                input_data += [float(data[i])]
            except ValueError:
                input_data += [str(data[i][1:-2])]

    # Set variable as a type of inputdata
    inData = m_model.inputdata()
    inData.mesh_type_ = input_data[0]
    inData.mesh_name_ = input_data[1]
    inData.lx_ = input_data[2]
    inData.ly_ = input_data[3]
    inData.nx_ = input_data[4]
    inData.ny_ = input_data[5]
    inData.bc_n_ = input_data[6]
    inData.bc_s_ = input_data[7]
    inData.bc_w_ = input_data[8]
    inData.bc_e_ = input_data[9]
    inData.ts_ = input_data[10]
    inData.dtw_ = input_data[11]
    inData.dtp_ = input_data[12]
    inData.dta_ = input_data[13]
    inData.temp_scheme_ = input_data[14]
    inData.spatial_scheme_ = input_data[15]
    inData.adapt_dt_ = input_data[16]
    inData.cfl_ = input_data[17]
    inData.heps_ = input_data[18]
    inData.friction_ = input_data[19]
    inData.g_ = input_data[20]
    inData.w_tecplot_ = input_data[21]
    inData.w_obs_ = input_data[22]
    inData.use_obs_ = input_data[23]
    inData.eps_min_ = input_data[24]
    inData.c_manning_ = input_data[25]
    inData.eps_manning_ = input_data[26]
    inData.regul_manning_ = input_data[27]
    inData.c_bathy_ = input_data[28]
    inData.eps_bathy_ = input_data[29]
    inData.regul_bathy_ = input_data[30]
    inData.c_hydrograph_ = input_data[31]
    inData.eps_hydrograph_ = input_data[32]
    inData.regul_hydrograph_ = input_data[33]


    ##### Read boudary conditions data ######
    # =========================================

    size_hydrograph = 241
    size_rat = 101

    # hydrograph.txt 
    a_file = open("hydrograph.txt", "r")
    lines = a_file.readlines()
    a_file.close()
    for i in range(len(lines)-size_hydrograph):
        del lines[0]
    new_file = open("hydrograph.txt", "w+")
    for line in lines:
        new_file.write(line)
    new_file.close()

    fname1 = os.path.join("hydrograph.txt")
    hyd = np.loadtxt(fname1) # Import average monthly precip to numpy array
    t1 = hyd[:,0]
    q1 = hyd[:,1]

    # rating_curve.txt
    a_file = open("rating_curve.txt", "r")
    lines = a_file.readlines()
    a_file.close()
    for i in range(len(lines)-size_rat):
        del lines[0]
    new_file = open("rating_curve.txt", "w+")
    for line in lines:
        new_file.write(line)
    new_file.close()

    fname2 = os.path.join("rating_curve.txt")
    rat = np.loadtxt(fname2) # Import average monthly precip to numpy array
    h2 = rat[:,0]
    q2 = rat[:,1]

    ##### Call model #####
    # ====================

    mdl = call_model.Model()
    mdl.mesh = m_mesh.msh()
    mdl.dof = m_model.unk(mdl.mesh)
    mdl.dof0 = m_model.unk(mdl.mesh)

    result = mdl.run_direct(t1, q1, h2, q2, inData) # Run direct

    ##### Set output ##### 
    # ====================

    size_ = result[0].size_series  

    h = []
    u = []
    v = []
    t = []
    for i in range(size_):
        h += [result[0].dof_[i].h]
        u += [result[0].dof_[i].u]
        v += [result[0].dof_[i].v]
        t += [result[0].dof_[i].t_display]
        print("Deep water h at time t = ",int(t[i]), " : ", h[i])
    cost_func = result[1]

    ##### Create a data file for Plotting data #####
    # ==============================================

    x=result[2].x_space
    y=result[2].y_space

    for j in range(size_):
        name = "t="+str(int(t[j]))
        filename = "%s.plt" % name
        filepath = os.path.join('./output_data', filename)
        if not os.path.exists('./output_data'):
            os.makedirs('./output_data')
        f=open(filepath,"w+")
        f.write('## VARIABLES  =  x  y  h  u  v ' + '\n')
        for i in range(len(x)):
            f.write( str(x[i]) + ' ' + str(y[i]) + ' ' + str(h[j][i]) + '  ' + str(u[j][i]) + '  ' + str(v[j][i]) + '\n' )
        f.close()

    print("==========================================")
    print("	SIMULATION PART FINISHED!	  	     ")
    print("==========================================")
    
# To run tests without pytest (debug)
if __name__ == "__main__":
    test_wrappers()
