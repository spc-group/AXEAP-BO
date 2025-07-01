from . import cowan_config
import os
import numpy as np
import shutil
from datetime import datetime
import glob
import sys
import matplotlib.pyplot as plt
import pandas as pd

# Function to generate sticks 

def gen_sticks(el, ox, fdd=100.0, fpd=100.0, gpd=100.0, soc=100.0, sov=100.0, tdq_g=0.0, tdq_e=0.0, dt_g=0.0, dt_e=0.0, ds_g=0.0, ds_e=0.0,
                m_g=0.0, m_e=0.0, em_start=None, em_end=None, G1=1.0, L1=1.2, L2=None, s_pt=None, norm=False):
    
    # Collect configuration parameters
    # ---------------------
    ion = el+ox
    AP = np.array([fdd, fpd, gpd, soc, sov])
    fname = f"{ion}_{fdd}_{fpd}_{gpd}_{soc}_{sov}"

    # Find proper configuration
    RunFileLocation = "C:/cowan"
    config_path = os.path.join(RunFileLocation, "Resources", "Configurations.txt")
    config_data = cowan_config.ConfigMaker(config_path)
    config_els = np.array([x.decode() for x in config_data[:,0]])
    det_ids = np.where(config_els == ion)
    detnum = det_ids[0]
    curr_config = np.squeeze(config_data[detnum])

    # Set cowan input state
    d34 = '3D'; dGr = '0'+ str(int(curr_config[3]))
    Pel = '1S'; eGr = '01'; Pel2 = '3P'
    init_state = Pel + eGr + ' ' + d34 + dGr # Initial electronic state of atom
    final_state = Pel2 + '05' + ' ' + d34 + dGr # Final electronic state of atom
    bindE = curr_config[4] - curr_config[8]
    atom_number = int(curr_config[1])

    # Set energy range values
    if not em_start:
        em_start = np.round(bindE-50, decimals=-1)
    if not em_end:
        em_end = np.round(bindE+50, decimals=-1)
    split = True if L2 else False
    
    # Calculate crystal field parameters
    # Indices: Dq_g = 0, Dt_g = 1, Ds_g = 2, Dq_e = 3, Dt_e = 4, Ds_e = 5, spin_g = 6, spin_e = 7
    CTVal = np.array([tdq_g, dt_g, ds_g, tdq_e, dt_e, ds_e, m_g, m_e]) 
    crystal_bf, crystal_ct, fano = cowan_config.CFG(CTVal) 
    
    # Run rcn programs
    # ---------------------

    # Write first ftn10 input
    val1 = -1;
    f1 = open('ftn10','w', encoding='utf-8')
    f1.write('22 -9    2   10  1.0    5.E-06    1.E-09-2   130   1.0  0.65  0.0 0.50 0.0  0.70 \n')
    f1.write('%s %u %s %s %s %s %s \n' % ('  ',atom_number,'   ', ion, init_state,'       ',init_state))
    f1.write('%s %u %s %s %s %s %s \n' % ('  ',atom_number,'   ', ion, final_state,'       ',final_state))
    f1.write('%s %d \n' % ('  ',-1))
    f1.close()

    # Run rcn31 and rcn2 programs
    os.system(os.path.join(RunFileLocation,'bin','rcn31.exe'))
    shutil.copyfile(os.path.join(RunFileLocation,'bin','RCN2.inp'),'ftn10')
    os.system(os.path.join(RunFileLocation,'bin','rcn2.exe'))
    shutil.copyfile('ftn11','ftn10')

    # Run rcg program
    # ---------------------

    # Update rcg file
    crystal = crystal_bf
    f_out = open('temp.rcg','w',encoding = 'utf-8')
    f1 = open('ftn10','r')
    f2 = open(os.path.join( RunFileLocation,'Resources','TemplateXES.rcg'),'r')
    for ii in range(0,3):
        #readline method is the same function as fgets in MATLAB
        f_out.write('%s' % (f2.readline()))
        f1.readline()
    tline1 = f1.readline()
    tline2 = f1.readline()
    d_req = len(tline2)

    if d_req>11:
        temp_1 = tline1[5:8] 
        temp_2 = tline1[10:13]
        tline1 = cowan_config.changeSTR(tline1,temp_2,5)
        tline1 = cowan_config.changeSTR(tline1,temp_1,10)
        temp_1 = tline2[5:8] 
        temp_2 = tline2[10:13]
        tline2 = cowan_config.changeSTR(tline2,temp_2,5)
        tline2 = cowan_config.changeSTR(tline2,temp_1,10)
    swap_SO = 0; #swapping SO coupling for F element
    pars_ini, line11, line12 = cowan_config.update_line_slater_out(AP, f1, 0, 0, swap_SO)
    pars_fin, line21, line22 = cowan_config.update_line_slater_out(AP, f1, 0, 1, swap_SO)

    f_out.write('%s' % tline1)
    f_out.write('%s' % tline2)
    f_out.write('%s' % line11)

    if (line12 != ''):
        f_out.write('%s' % line12)
    f_out.write('%s' % line21)

    if (line22 !=''):
        f_out.write('%s' % line22)    
    while 1:
        tline = f1.readline()
        if (len(tline) == 0):
            break
        f_out.write('%s' % tline)

    f_out.close(); f1.close(); f2.close()
    
    shutil.copyfile('temp.rcg','ftn10')
    shutil.copyfile('ftn10',fname+'.rcg')
    os.system(os.path.join(RunFileLocation,'bin','rcg9.exe'))
    shutil.move('ftn14','temp.m14')

    f1 = open(os.path.join(RunFileLocation,'Resources','template.rac'),'r')
    f2 = open('temp.rac','w',encoding = 'utf-8')
    crystal[3] = crystal[3]/1000
    crystal[8] = crystal[8]/1000

    num_lines = 39   
    str_nums = np.array([12, 13, 14, 16, 18, 23, 24, 25, 27, 29])-1
    for jj in range(0, num_lines):
        dd = f1.readline()
        curr_str = np.argwhere(str_nums == jj)  
        if not curr_str.size>0:
            f2.write('%s' % dd)
        else:
            f2.write('%s' % dd[:35])
            f2.write('% 3.3f\n' % crystal[curr_str])
            
    f1.close()
    f2.close()    

    os.system(os.path.join(RunFileLocation,'bin','racer.exe temp.m14 temp.ora <temp.rac'))
    shutil.copyfile('temp.ora',fname+'.oea')
    shutil.copyfile('temp.rac',fname+'.rac')

    # Create_nfo_file
    # fname_nfo = File1+'.nfo'
    # fnfo = open(fname_nfo,'w',encoding = 'utf-8')
    # fnfo.write('%s\t %s\n' % ('Filename:',File1))
    # fnfo.write('%s\t\t %s\n' % ('Ion:', config))
    # fnfo.write('%s\t %s\n' % ('Spectrum:','1s3p XES')) 
    # fnfo.write('%s\t %5.1f\n' % ('Binding energy (eV):',bindE))    
    numHoles = int(10 - curr_config[3])
    step_p_ini = np.shape(pars_ini); sw2_4 = 0
    F2dd_ini = 0; F4dd_ini = 0; LS3d_ini = 0

    for jj in range(0, step_p_ini[0]):
        detnum = pars_ini[jj][0]
        if (detnum == 1):
            if (sw2_4 == 0):
                F2dd_ini = pars_ini[jj][1]*0.8
                sw2_4 = 1
            else:
                F4dd_ini = pars_ini[jj][1]*0.8
        elif(detnum == 2):
            LS3d_ini = pars_ini[jj][1]

    F2dd_fin = 0; F4dd_fin = 0; LS2p = 0; LS3d_fin = 0; F2pd = 0
    G1pd = 0; G3pd = 0
    step_p_fin = np.shape(pars_fin); sw2_4 =0; sw1_3 = 0;

    for jj in range(0, step_p_fin[0]):
        detnum = pars_fin[jj][0]
        if(detnum == 1):
            if (sw2_4 == 0):
                F2dd_fin = pars_fin[jj][1]*0.8
                sw2_4 = 1
            else:
                F4dd_fin = pars_fin[jj][1]*0.8
        elif(detnum == 2):
            if (pars_fin[jj][1]>0.4):
                LS2p = pars_fin[jj][1]
            else:
                LS3d_fin = pars_fin[jj][1]
        elif(detnum == 3):
            F2pd = pars_fin[jj][1]*0.8
        elif(detnum == 4):
            if (sw1_3 == 0):
                G1pd = pars_fin[jj][1]*0.8
                sw1_3 = 1;
            else:
                G3pd = pars_fin[jj][1]*0.8
                    
    # fnfo.write('%s\n' % ' ')
    # sl_Fdd = '(' + str(AP[0]) + '%)' 
    # sl_Fpd = '(' + str(AP[1]) + '%)'
    # sl_Gpd = '(' + str(AP[2]) + '%)'
    # so_2p = '(' + str(AP[3]) + '%)'
    # so_3d = '(' + str(AP[4]) + '%)'
    # fnfo.write('%s\t\t %5.3f %s\t  %5.3f %s\n' % ('F2dd:', F2dd_ini,sl_Fdd,F2dd_fin,sl_Fdd))
    # fnfo.write('%s\t\t %5.3f %s\t  %5.3f %s\n' % ('F4dd:', F4dd_ini, sl_Fdd,F4dd_fin, sl_Fdd))
    # fnfo.write('%s\t\t %5.3f %s\t  %5.3f %s\n' % ('LS3d:',LS3d_ini,so_3d, LS3d_fin,so_3d))
    # fnfo.write('%s\t\t\t\t  %5.3f %s\n' % ('LS2p:',LS2p,so_2p))
    # fnfo.write('%s\t\t\t\t  %5.3f %s\n' % ('F2pd:',F2pd,sl_Fpd))
    # fnfo.write('%s\t\t\t\t  %5.3f %s\n' % ('G1pd:',G1pd,sl_Gpd))
    # fnfo.write('%s\t\t\t\t  %5.3f %s\n' % ('G3pd:',G3pd,sl_Gpd))
    # fnfo.write('%s\n' % ' ')
    # fnfo.write('%s\t\t %5.3f \t\t  %5.3f\n' % ('10Dq:',crystal_ct[0], crystal_ct[5]))
    # fnfo.write('%s\t\t %5.3f \t\t  %5.3f\n' % ('Dt:', crystal_ct[1], crystal_ct[6]))
    # fnfo.write('%s\t\t %5.3f \t\t  %5.3f\n' % ('Ds:', crystal_ct[2], crystal_ct[7]))
    # fnfo.write('%s\t\t %5.3f \t\t  %5.3f\n' % ('M(meV):', crystal_ct[3], crystal_ct[8]))
    # fnfo.write('%s\t\t %5.3f\n' % ('Lorenz1:', L1))
    # if SplitTF:
    #     fnfo.write('%s\t\t %5.3f\n' % ('Lorenz2:', L2))
    #     fnfo.write('%s\t\t %5.3f\n' % ('Splitting Point:', SplitPoint))
    # fnfo.write('%s\t\t %5.3f\n' % ('Gauss:', G1))
    # currentSecond= str(datetime.now().second) 
    # currentMinute = str(datetime.now().minute) + ':'
    # currentHour = str(datetime.now().hour) + ':'
    # currentDay = str(datetime.now().day) + ','
    # currentMonth = str(datetime.now().month) + '.'
    # currentYear = str(datetime.now().year) + '.'
    # fnfo.write('%s\t %s\n' % ('Num Holes:', str(numHoles)))
    # fnfo.write('%s %s %s %s %s %s' % (currentYear, currentMonth, currentDay, currentHour, currentMinute, currentSecond))
    # fnfo.close()
        
    dlist1 = glob.glob('ftn*') 
    for jj in dlist1:
        os.remove(jj)
    dlist2 = glob.glob('temp*')
    for jj in dlist2:
        os.remove(jj)
    
    # Stick Generation   
    # ---------------------


    files2plot = fname+'.oea'
    if not os.path.exists(files2plot):
        print('Failed to calculate XES.')
        sys.exit()

    # print('Let us keep going!')
    shutil.copyfile(files2plot,'temp.ora')

    # fid_NFO = open(File1+'.nfo','r')
    BEflag = 1;
    # for jj in range(0,8):
    #     bindln = fid_NFO.readline()
    #     if (bindln[0:2] == 'Bi'):
    #         break
    bindln = '%s\t %5.1f\n' % ('Binding energy (eV):',bindE)
    bindEs = bindln[20:]
    bindEs = str(float(bindEs))
    # fid_NFO.close()

    if split:
        fid1 = open(os.path.join(RunFileLocation,'Resources','template_split.plo'),'r')
    else:
        fid1 = open(os.path.join(RunFileLocation,'Resources','template.plo'),'r')
    fid2 = open('temp.plo','w',encoding = 'utf-8')

    if split:
        dd = fid1.readline()
        fid2.write('%s' % dd[:11])
        fid2.write('%3.2f %3.2f' % (L1, fano))
        fid2.write('%s %4.0f\n' % (' range 0', s_pt))
        dd = fid1.readline()
        fid2.write('%s' % dd[:11])
        fid2.write('%3.2f %3.2f' % (L2, fano))
        fid2.write('%s %4.0f %s\n' % (' range', s_pt,' 9999'))
    else:
        dd = fid1.readline()
        fid2.write('%s' % dd[:11])
        fid2.write('%3.2f %3.2f\n' % (L1, fano))
        
    dd = fid1.readline()
    fid2.write('%s' % dd[:9])
    fid2.write('%3.2f\n' % G1)

    #Force range
    fid2.write('%s' % 'energy_range ')
    fid2.write('%3.2f %3.2f\n' % (em_start, em_end))

    #Account for the temperature
    spectra = 3
    cmd_plot = 'old_racah'
    temp = ''
    fid2.write('%s %s' % (cmd_plot,' '))
    fid2.write('%s\n' % 'temp.ora')
    dd = fid1.readline()
    fid2.write('%s' % dd)
    len_spec_line = 29;

    for ii in range (0, spectra):
        dd = fid1.readline()
        fid2.write('%s' % dd)
        dd = fid1.readline()
        line_st = dd[:len_spec_line]+bindEs+temp
        fid2.write('%s\n' % line_st)
        
    dd = fid1.readline()
    fid2.write('%s' % dd)
    fid1.close(); fid2.close()

    os.system(os.path.join(RunFileLocation,'bin','plo1.exe <temp.plo'))
    shutil.copyfile('temp.plo', fname+'.plo')
        
    fexist = 'temp_left.dat'
    if not os.path.exists(fexist):
        print('Failed to create sticks.')
        sys.exit()

    # print('Let us extract stick and spectrum')

    fid = open(fexist, 'r')
    for line in range(0,1000):
        dd = fid.readline()
        specstr = dd.split()
        if (line == 0):
            x1 = float(specstr[0])
            y1 = float(specstr[1])
        else:
            x1 = np.append(x1, float(specstr[0]))
            y1 = np.append(y1, float(specstr[1]))

    for line in range (0,2):
        dd = fid.readline()
        
    k = 0
    while 1:
        k = k+1
        dd = fid.readline()
        stickstr = dd.split()
        if (len(dd)==0):
            break
        else:
            if (k == 1):
                x2 = float(stickstr[0])
                y2 = float(stickstr[1])
            else:
                x2 = np.append(x2, float(stickstr[0]))
                y2 = np.append(y2, float(stickstr[1]))
        if k>1000:
            print('Can not find sticks')
            break
    # x1, y1 = spectrum data :::: x2, y2 = stick data 
    fid.close()

    fexist = 'temp_right.dat'
    if not os.path.exists(fexist):
        print('Failed to create sticks.')
        sys.exit()

    # print('Let us extract stick and spectrum')

    fid = open(fexist, 'r')
    for line in range(0,1000):
        dd = fid.readline()
        specstr = dd.split()
        if (line == 0):
            x1 = float(specstr[0])
            y1 = float(specstr[1])
        else:
            x1 = np.append(x1, float(specstr[0]))
            y1 = np.append(y1, float(specstr[1]))

    for line in range (0,2):
        dd = fid.readline()
        
    k = 0
    while 1:
        k = k+1
        dd = fid.readline()
        stickstr = dd.split()
        if (len(dd)==0):
            break
        else:
            ener = float(stickstr[0])
            inten = float(stickstr[1])
            if np.count_nonzero(x2 == ener):
                y2[np.argwhere(x2 == ener)] += inten
            else:
                x2 = np.append(x2, ener)
                y2 = np.append(y2, inten)
        if k>1000:
            print('Can not find sticks')
            break
    # x1, y1 = spectrum data :::: x2, y2 = stick data 
    fid.close()

    fexist = 'temp_zero.dat'
    if not os.path.exists(fexist):
        print('Failed to create sticks.')
        sys.exit()

    # print('Let us extract stick and spectrum')

    fid = open(fexist, 'r')
    for line in range(0,1000):
        dd = fid.readline()
        specstr = dd.split()
        if (line == 0):
            x1 = float(specstr[0])
            y1 = float(specstr[1])
        else:
            x1 = np.append(x1, float(specstr[0]))
            y1 = np.append(y1, float(specstr[1]))

    for line in range (0,2):
        dd = fid.readline()
        
    k = 0
    while 1:
        k = k+1
        dd = fid.readline()
        stickstr = dd.split()
        if (len(dd)==0):
            break
        else:
            ener = float(stickstr[0])
            inten = float(stickstr[1])
            if np.count_nonzero(x2 == ener):
                y2[np.argwhere(x2 == ener)] += inten
            else:
                x2 = np.append(x2, ener)
                y2 = np.append(y2, inten)
        if k>1000:
            print('Can not find sticks')
            break
    # x1, y1 = spectrum data :::: x2, y2 = stick data 
    fid.close()


    bindE = float(bindEs)
    #y1 = np.flipud(y1); Was commented out originally
    x1 = -x1 +2*bindE
    if norm:
        y2 = y2/np.trapz(y2, x = x2) # normalization
    x2 = -x2 +2*bindE
    y2 = np.flipud(y2); 
    x2 = np.flipud(x2); 
    y2 = y2/np.sum(y2)

    dlist = glob.glob('temp*')
    for jj in dlist:
        os.remove(jj)
    dlist = glob.glob('*.rac')
    for jj in dlist:
        os.remove(jj)
    dlist = glob.glob('*.plo') 
    for jj in dlist:
        os.remove(jj)
    dlist = glob.glob('*.rcg') 
    for jj in dlist:
        os.remove(jj)
    dlist = glob.glob('*.oea') 
    for jj in dlist:
        os.remove(jj)

    # df1 = pd.DataFrame({'Energy': x1, 'Intensity': y1})
    # df2 = pd.DataFrame({'Energy': x2, 'Intensity': y2})

    # writer = pd.ExcelWriter(File1+'.xlsx', engine = 'xlsxwriter')
    # df1.to_excel(writer,sheet_name='Spectrum')
    # df2.to_excel(writer,sheet_name='Sticks')
    # writer.close()

    # fig, ax = plt.subplots()
    # ax.plot(x1,y1)
    # ax.bar(x2,y2/5, color = 'maroon', width = 0.2)
    # plt.ylim(0,0.2)
    # plt.show()
    return x2, y2, x1, y1

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

