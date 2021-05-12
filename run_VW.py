#14:40:04, Wed, 12-May-2021, By Haibo
# python script to run and get VW of 1D tissue

from subprocess import call
import math
import sys
import numpy as np


def run_SingleCell_For_IC(Mode, Para, AF_model, block_para, Drug_ID, f):

	P=['./NCZ_Model_for_IC', str(BCL), '14', str(BCL), 'WT', Mode, '0', '1']
	for i in range(len(Para)):
		P.append(str(Para[i]))
	for i in range(len(AF_model)):
		P.append(str(AF_model[i]))
	for i in range(len(block_para)):
		P.append(str(block_para[i]))

	P.append(str(Drug_ID));
	call(P, stdout=f)



def run_Normal(Para, AF_model, block_para,f,Ini_cond_File, ID=0):

	call(['cp ICs/1Hz_ISO.0.1.CKII_ih.0.CKII_db.1/ICs.bin.%d ICs/ICs.bin.1000.AF.0.ISO.0.1.CKII_ih.0.CKII_db.1'%ID], shell=True)

	P=['./ONE_D_MPI_Ghost', 'ICs', 'Specific', 'S1_number', '6', 'S2', str(360)] + sim_control

	P.append('Popl_SF_Number')
	P.append('25')
	P.append('Popul_scaling_factors')
	for i in range(len(Para)):
		P.append(str(Para[i]))
	print(P)

	call(P)

	call(['mv  OneD_output.dat.0 cv_out_file.dat OneD_output_diff_coef.dat.0 %s'%sub_folder], shell=True)
	call(['cp -r  Tissue_ICs %s'%sub_folder], shell=True)



def run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0):
                                            # 'CaMKII', 'CaMKII_db',
	P=['./ONE_D_MPI_Ghost_VW', 'ICs', 'Restart', 'Time_Start', '5250', 'S1_number', '6', 'S2', str(S2)] + sim_control

	P.append('Popl_SF_Number')
	P.append('25')
	P.append('Popul_scaling_factors')
	for i in range(len(Para)):
		P.append(str(Para[i]))
	# print(P)


	P.append('OneD_OutFile')
	P.append(f)
	P.append('Sim_ID')
	P.append(str(ID))
	print( P)
	call(P)
	with open('conduction.log',"r") as ff:
		line1 = ff.readline().strip();
		line2 = ff.readline().strip();
	if int(line1) not in [0,1,2] or int(line2) not in [0,1,2]:
		print ('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> error in reading ''conduction.log'', ID = ', ID, 'line1 = ', line1, 'line2 = ', line2,  file=sys.stderr)
		sys.exit(0);


	return int(line1), int(line2);



def check_if_bilateral_cond(c):
	a=c[0]
	b=c[1]
	if (a==2-1 and b == 2-1):
		go=True
	# elif (a>=2-1 and b >= 2-1):  # if last beat: check !!!
	# 	go= True
	else:
		go = False
	return go

def check_if_bilateral_block(c):
	a=c[0]
	b=c[1]
	if (a==1-1 and b == 1-1):
		go=True
	else:
		go=False
	return go


def get_upper_bound_binary_search(para_set, AF_set, block_set,OneD_files, Cell_ID, S2):
	attempts = {}

	BCL = 1000
	d_S2 = 80
	f = open('workfile_upper', 'a')
	
	go = check_if_bilateral_cond(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
	lastgo = go
	print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f, flush=True)

	attempts[S2] = go  # keep record
	compute = True

	up_bound = -1
	low_bound = -1

	while compute:
		if(len(list(set(list(attempts.values())))) == 1):  # if all values are the same
			S2 = S2 - d_S2*go + d_S2*(not go)

		else:
			up_bound = np.min([k for k, v in attempts.items() if v == True])
			low_bound = np.max([k for k, v in attempts.items() if v == False])
			S2 = (up_bound + low_bound) / 2.0
			if (up_bound - low_bound < 0.1):
				# print( S2, up_bound, low_bound )
				print ( "S1$%^= ", BCL, "S2= ", S2, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID, 'number attempt = ', len(attempts), file=f, flush=True)
				break

		go = check_if_bilateral_cond(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
		attempts[S2] = go  # keep record	

		if len(attempts) > 30:
			compute = False
			S2 = np.nan
		print ( "S1= ", BCL, "S2= ", S2, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID,'number attempt = ', len(attempts),  file=f, flush=True)

	call(['mv OneD_output.dat.0 OneD_output.dat.0.S2_upper.%f'%S2], shell=True)
	call(['mv OneD_output.dat.0.S2_upper.%f %s'%(S2,sub_folder)], shell=True)
	
	return S2, up_bound, low_bound 



def get_upper_bound(para_set, AF_set, block_set,OneD_files,Cell_ID, S2):

	f = open('workfile_upper', 'a')
	
	go = check_if_bilateral_cond(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
	lastgo = go
	# go = reduce_S2(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))

	# run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0)
	# call('mv OneD_output.dat S2.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
	compute = True
	d_S2 = 40;
	while compute:
			if(lastgo == (not go) ):
				d_S2=d_S2/2.0

			if go:
				S2 =S2-d_S2 +0.01
			else:
				S2 =S2+d_S2 -0.01
			lastgo=go
			go=check_if_bilateral_cond(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
			print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f)
			if d_S2<0.1 and go:
				compute = False
			if S2> BCL or S2 < 50:
				compute= False
				S2 = np.nan
			if np.isnan(S2):
				compute=False
	call(['mv OneD_output.dat.0 OneD_output.dat.0.S2_upper.%f'%S2], shell=True)
	call(['mv OneD_output.dat.0.S2_upper.%f %s'%(S2,sub_folder)], shell=True)

	return S2,d_S2


def get_lower_bound_binary_search(para_set, AF_set, block_set,OneD_files,Cell_ID, S2):
	attempts = {}

	BCL = 1000
	d_S2 = 5
	f = open('workfile_lower', 'a')
	
	go = check_if_bilateral_block(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
	lastgo = go
	print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f,  flush=True)

	attempts[S2] = go  # keep record
	compute = True

	up_bound = 0
	low_bound = 0
	while compute:
		if(len(list(set(list(attempts.values())))) == 1):  # if all values are the same
			S2 = S2 + d_S2*go - d_S2*(not go)

		else:
			up_bound = np.min([k for k, v in attempts.items() if v == False])
			low_bound = np.max([k for k, v in attempts.items() if v == True])
			S2 = (up_bound + low_bound) / 2.0
			if (up_bound - low_bound < 0.1):
				# print( S2, up_bound, low_bound )
				print ( "S1$%^= ", BCL, "S2= ", S2,"go = ", go, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID, 'number attempt = ', len(attempts), file=f,  flush=True)
				break

		go = check_if_bilateral_block(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
		attempts[S2] = go  # keep record	

		if len(attempts) > 30:
			compute = False
			S2 = np.nan
		print ( "S1= ", BCL, "S2= ", S2, "go =", go,"low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID,'number attempt = ', len(attempts),  file=f,  flush=True)

	call(['mv OneD_output.dat.0 OneD_output.dat.0.S2_lowwer.%f'%S2], shell=True)
	call(['mv OneD_output.dat.0.S2_lowwer.%f %s'%(S2,sub_folder)], shell=True)
	
	return S2, up_bound, low_bound 




def get_lowwer_bound(para_set, AF_set, block_set,OneD_files,Cell_ID, S2):

	f = open('workfile_lowwer', 'a')
	
	go = check_if_bilateral_block(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
	lastgo = go
	# go = reduce_S2(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))

	# run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0)
	# call('mv OneD_output.dat S2.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
	compute = True
	d_S2 = 5;
	while compute:
			if(lastgo == (not go) ):
				d_S2=d_S2/2.0
			if go:
				S2 =S2+d_S2 -0.01
			else:
				S2 =S2-d_S2 +0.01
			lastgo=go
			go=check_if_bilateral_block(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))
			print("S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID,file=f)
			if d_S2<0.1 and go:
				compute = False
			if S2> BCL or S2 < 50:
				compute= False
				S2 = np.nan
			if np.isnan(S2):
				compute=False

	call(['mv OneD_output.dat.0 OneD_output.dat.0.S2_lowwer.%f'%S2], shell=True)
	call(['mv OneD_output.dat.0.S2_lowwer.%f %s'%(S2,sub_folder)], shell=True)

	return S2, d_S2


files=[]
singlecell_file=[]
OneD_files=[]
OneD_files.append("OneD.log.dat."+str(0))
OneD_files.append("OneD.log.dat."+str(1))



BCL=int(sys.argv[1]);
ISO_con = float(sys.argv[2]);
CaMKII = sys.argv[3]#'CaMKII_db'

sim_control = ['BCL', str(BCL), 'ISO_con', str(ISO_con), 'Total_time', '6550'] #+ ['CaMKII_K1_NaV_GapG', 'No_CaMKII_GapG']
if(CaMKII):
	sim_control = sim_control + ['CaMKII', CaMKII]
print(sim_control)
print ('Simlation control: ', sim_control)



S2 = 500
# go = True
para_set = [ [1] * 25 ]  #np.loadtxt('para.log')  #[ [1] * 25 ] #
AF_set = [0,0]
block_set = [0,0]
Drug_ID=0
AF_ID=0

folder = 'BCL.%.0f.ISO_con.%.1f.%s'%(BCL, ISO_con,CaMKII)
call(['mkdir' , folder])



for  Cell_ID in range(1):
	sub_folder = folder+'/ID_%d'%Cell_ID;
	call(['mkdir' , sub_folder])
	# pass
	# 319.51125,0.078125# 
	# run_Normal(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[0],'None',0)

	# run_VW(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[0], 'Specific', S2, Cell_ID )
	#### run_Normal(Para, AF_model, block_para,f,Ini_cond_File, ID=0):
	S2_upper, S2_upper_up_bound, S2_upper_low_bound = get_upper_bound_binary_search(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[0],Cell_ID, S2)
	S2_lowwer, S2_lower_up_bound, S2_lower_low_bound = get_lower_bound_binary_search(para_set[Cell_ID], AF_set[AF_ID], block_set[Drug_ID],OneD_files[0],Cell_ID, S2_upper-10)

	call(['mv workfile_lowwer workfile_upper %s'%sub_folder], shell=True)

	vw_results = open(folder+"/vw.results.ID.%d"%Cell_ID, "w+")
	print (Cell_ID, S2_upper,S2_upper_up_bound, S2_upper_low_bound, S2_lowwer, S2_lower_up_bound, S2_lower_low_bound, file=vw_results)
	vw_results.close()
