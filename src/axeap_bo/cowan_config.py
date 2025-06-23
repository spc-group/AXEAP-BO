import numpy as np
import math


def ConfigMaker(config_path):
    f = open(config_path,'r')
    # Storing all the columns in configurations.txt. What do these columns represent? 
    C = np.loadtxt(f, delimiter ='\t', dtype = {'names': ('col1', 'col2', 'col3', 'col4', 'col5',
                                                          'col6', 'col7', 'col8', 'col9','col10',
                                                          'col11', 'col12', 'col13', 'col14','col15',
                                                          'col16', 'col17', 'col18', 'col19','col20',
                                                          'col21', 'col22', 'col23', 'col24'),
                                                'formats': ('S4', 'float','float','float','float', 
                                                           'float', 'float','float','float','float',
                                                           'float', 'float','float','float','float',
                                                           'float', 'float','float','float','float',
                                                           'float', 'float','float','float')})
    f.close()
    #if you want some specific column in C, use C['colname']
    config_temp = np.empty(shape = (C.size, 16), dtype = object)
    # Modifying the config data to return... why is this necessary?
    for jj in range(0,C.size):
        temp = C['col1']; config_temp[jj,0] = temp[jj]
        temp = C['col2']; config_temp[jj,1] = temp[jj]
        temp = C['col3']; config_temp[jj,2] = temp[jj]
        temp = C['col4']; config_temp[jj,3] = temp[jj]
        temp = C['col5']; valence = temp[jj]; corr_fac = (valence-2)*1.5
        #1s
        temp = C['col6']; config_temp[jj,4] = temp[jj]+corr_fac
        #2s
        temp = C['col7']; config_temp[jj,5] = temp[jj]+corr_fac
        #2p
        temp = C['col8']; temp1 = C['col9']
        config_temp[jj,6] = (temp[jj]+2*temp1[jj])/3+corr_fac
        #3s
        temp = C['col10']; config_temp[jj,7] = temp[jj]+corr_fac
        #3p
        temp = C['col11']; temp1 = C['col12']
        config_temp[jj,8] = (temp[jj]+2*temp1[jj])/3+corr_fac
        #3d
        temp = C['col13']; temp1 = C['col14']
        config_temp[jj,9] = (2*temp[jj]+3*temp1[jj])/5+corr_fac
        #4s
        temp = C['col15']; config_temp[jj,10] = temp[jj]+corr_fac
        #4p
        temp = C['col16']; temp1 = C['col17']
        config_temp[jj,11] = (temp[jj]+2*temp1[jj])/3+corr_fac
        #4d
        temp = C['col18']; temp1 = C['col19']
        config_temp[jj,12] = (2*temp[jj]+3*temp1[jj])/5+corr_fac
        #5s
        temp = C['col20']; config_temp[jj,13] = temp[jj]+corr_fac
        #5p
        temp = C['col21']; temp1 = C['col22']
        config_temp[jj,14] = (temp[jj]+2*temp1[jj])/3+corr_fac
        #5d
        temp = C['col23']; temp1 = C['col24'];
        config_temp[jj,15] = (2*temp[jj]+3*temp1[jj])/5+corr_fac       
    return config_temp

# First part of the XES Calculation function - looks good
def CFG(input):
    crystal_ct = np.empty(shape = (10), dtype = float)
    crystal = np.empty(shape = (10), dtype = float)
    Dq_g = input[0]; Dt_g = input[1]; Ds_g = input[2]; Dq_e = input[3]
    Dt_e = input[4]; Ds_e = input[5]; spin_g = input[6]; spin_e = input[7]
    
    # Same array assignment as function
    crystal_ct[0]=Dq_g; crystal_ct[1]=Dt_g; crystal_ct[2]=Ds_g 
    crystal_ct[3]=spin_g; crystal_ct[4]=0; crystal_ct[5]=Dq_e 
    crystal_ct[6]=Dt_e; crystal_ct[7]=Ds_e; crystal_ct[8]=spin_e 
    crystal_ct[9]=0
    
    # Looks like we are assuming that symmetry is less than 4
    crystal[0]=6*math.sqrt(30)*Dq_g/10-3.5*math.sqrt(30)*Dt_g
    crystal[1]=-2.5*math.sqrt(42)*Dt_g
    crystal[2]=-math.sqrt(70)*Ds_g
    crystal[3]=spin_g
    crystal[4] = 0
    crystal[5]=6*math.sqrt(30)*Dq_e/10-3.5*math.sqrt(30)*Dt_e
    crystal[6]=-2.5*math.sqrt(42)*Dt_e
    crystal[7]=-math.sqrt(70)*Ds_e
    crystal[8]=spin_e
    crystal[9] = 0
    fano=99999
    return [crystal, crystal_ct, fano]

def changeSTR(fulSTR, chanSTR,idx):
    A = list(fulSTR)
    Alen = len(chanSTR)
    A[idx:idx+Alen] = chanSTR
    modiSTR = ''.join(A)
    return modiSTR

def update_line_slater_out(AtomicParam, fileID, bindE, fin_state, swap_SO):
# MATLAB Code = [pars_ini, line11, line12]=update_line_slater_out(app,fid1, bindE, fin_state, swap_SO,AP);
    line_in1= fileID.readline() # Reading in the appropriate state of the file (ftn10)
    num_param = int(line_in1[18:20]) # Second value in the line
    SO1=-1;SO2=-1; numSO=0; SOcounter=0;
    for jj in range(4, 8):
        indx = jj * 10;
        reduced_val = float(line_in1[indx-8:indx-1]) # 3-7 values in the line
        if (int(line_in1[indx-1]) == 2): # If the last digit is a 2
            numSO = numSO+1;
            if (SO1 == -1):
                SO1 = reduced_val
            else:
                SO2 = reduced_val
    if (num_param>5): # Seems to apply to only the final state (at least for Mn2+)
        line_in2 = fileID.readline()
        for jj in range(1,num_param-4):
            indx = jj*10;
            reduced_val = float(line_in2[(indx-8):(indx-1)])
            if (int(line_in2[indx-1]) == 2):
                numSO = numSO+1;
                if (SO1 == -1):
                    SO1 = reduced_val
                else:
                    SO2 = reduced_val
                    
    if (SO2 == -1):
        if (SO1 != -1):
            SO1 = SO1*AtomicParam[4]/100 # Multiply by valence percentage, which is the fifth parameter
    else:
        if (SO1>SO2):
            SO1 = SO1*AtomicParam[3]/100
            SO2 = SO2*AtomicParam[4]/100
        else:
            SO2 = SO2*AtomicParam[3]/100
            SO1 = SO1*AtomicParam[4]/100
    
    if (swap_SO and numSO > 1 ):
        tempSO = SO1; SO1=SO2; SO2=tempSO
        
    pars_out = np.empty(shape=(4,2), dtype = float) #pre-allocation
    
    for jj in range(4, 8):
        indx = jj*10;
        #Every 10 step has importan information (Warinig: this is not an atomic value)
        #e.g. 12.3451   7.9791 where atomic values are only 12.345 and 7.979
        reduced_val = float(line_in1[(indx-8):(indx-1)])
        detnum = line_in1[indx-1]
        if (detnum == '0'):
            reduced_string = '0';
        elif (detnum == '1'):
            line_in1 = changeSTR(line_in1,'         ',indx-9)
            reduced_val = reduced_val*AtomicParam[0]/100
            reduced_string = "{:5.3f}".format(reduced_val)
            l_led_str = len(reduced_string)
            line_in1 = changeSTR(line_in1,reduced_string+detnum,indx-l_led_str-1)
        elif (detnum == '2'):
            line_in1 = changeSTR(line_in1,'         ',indx-9)
            SOcounter = SOcounter +1
            if (SOcounter == 1):
                reduced_val = SO1
            else:
                reduced_val = SO2
            reduced_string = "{:5.3f}".format(reduced_val)
            l_led_str = len(reduced_string)
            line_in1 = changeSTR(line_in1,reduced_string+detnum,indx-l_led_str-1)
        elif (detnum == '3'):
            line_in1 = changeSTR(line_in1,'         ',indx-9)
            reduced_val = reduced_val*AtomicParam[1]/100
            reduced_string = "{:5.3f}".format(reduced_val)
            l_led_str = len(reduced_string)
            line_in1 = changeSTR(line_in1,reduced_string+detnum,indx-l_led_str-1)
        elif (detnum == '4'):
            line_in1 = changeSTR(line_in1,'         ',indx-9)
            reduced_val = reduced_val*AtomicParam[2]/100
            reduced_string = "{:5.3f}".format(reduced_val)
            l_led_str = len(reduced_string)
            line_in1 = changeSTR(line_in1,reduced_string+detnum,indx-l_led_str-1)
        
        pars_out[jj-4][0] = float(detnum)
        pars_out[jj-4][1] = float(reduced_string)
    
    if bindE> -1:
        bindEs = "{:9.3f}".format(bindE)
        lBE = len(bindEs)
        detnum = line_in1[29]
        line_in1 = changeSTR(line_in1,bindEs+detnum,28-lBE+1)
    
    if (fin_state == 1):
        temp_1 = line_in1[50:60] 
        temp_2 = line_in1[60:70]
        line_in1 = changeSTR(line_in1,temp_2,50)
        line_in1 = changeSTR(line_in1,temp_1,60)
    else:
        line_in1 = changeSTR(line_in1,'     0.000',60)
    
    if num_param>5:
        pars_out2 = np.empty(shape=(num_param-5,2), dtype = float)
        for jj in range(1, num_param-4):
            indx = jj*10;
            reduced_val = float(line_in2[(indx-8):(indx-1)])
            detnum = line_in2[indx-1]
            if (detnum == '1'):
                line_in2 = changeSTR(line_in2,'         ',indx-9)
                reduced_val = reduced_val*AtomicParam[0]/100
                reduced_string = "{:5.3f}".format(reduced_val)
                l_led_str = len(reduced_string)
                line_in2 = changeSTR(line_in2,reduced_string+detnum,indx-l_led_str-1)
            elif (detnum == '2'):
                line_in2 = changeSTR(line_in2,'         ',indx-9)
                SOcounter = SOcounter +1
                if (SOcounter == 1):
                    reduced_val = SO1
                else:
                    reduced_val = SO2
                reduced_string = "{:5.3f}".format(reduced_val)
                l_led_str = len(reduced_string)
                line_in2 = changeSTR(line_in2,reduced_string+detnum,indx-l_led_str-1)
            elif (detnum == '3'):
                line_in2 = changeSTR(line_in2,'         ',indx-9)
                reduced_val = reduced_val*AtomicParam[1]/100
                reduced_string = "{:5.3f}".format(reduced_val)
                l_led_str = len(reduced_string)
                line_in2 = changeSTR(line_in2,reduced_string+detnum,indx-l_led_str-1)
            elif (detnum == '4'):
                line_in2 = changeSTR(line_in2,'         ',indx-9)
                reduced_val = reduced_val*AtomicParam[2]/100
                reduced_string = "{:5.3f}".format(reduced_val)
                l_led_str = len(reduced_string)
                line_in2 = changeSTR(line_in2,reduced_string+detnum,indx-l_led_str-1)
            pars_out2[jj-1][0] = float(detnum)
            pars_out2[jj-1][1] = float(reduced_string)
        pars_out = np.append(pars_out,pars_out2,axis = 0)
    else:
        line_in2 = ''
    return pars_out, line_in1, line_in2