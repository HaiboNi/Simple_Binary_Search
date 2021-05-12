import numpy as np

from subprocess import call

def check_if_bilateral_cond(S2):
	target = 578.888
	if S2 - target >= 0.001:
		return True
	else:
		return False


def check_if_bilateral_block(S2):
	target = 566.888
	if S2 - target <= 0.001:
		return True
	else:
		return False



def get_upper_bound(Cell_ID, S2):

	f = open('workfile_upper', 'a')
	
	go = check_if_bilateral_cond(S2)
	lastgo = go
	BCL = 1000
	# go = reduce_S2(run_VW(para_set, AF_set, block_set,OneD_files, 'Specific', S2, Cell_ID ))

	# run_VW(Para, AF_model, block_para,f,Ini_cond_File, S2, ID=0)
	# call('mv OneD_output.dat S2.BCL.3Hz.Popul.Index.'+str(Cell_ID)+'.AF_Model.' + str(AF_ID) + '.block.'+str(Drug_ID), shell=True);
	compute = True
	d_S2 = 80;
	while compute:
			if(lastgo == (not go) ):
				d_S2=d_S2/2.0

			if go:
				S2 =S2-d_S2 +0.01
			else:
				S2 =S2+d_S2 -0.01
			lastgo=go
			go=check_if_bilateral_cond(S2)
			print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f)
			if d_S2<0.1 and go:
				compute = False
			if S2> BCL or S2 < 50:
				compute= False
				S2 = np.nan
			if np.isnan(S2):
				compute=False
	# call(['mv OneD_output.dat.0 OneD_output.dat.0.S2_upper.%f'%S2], shell=True)
	# call(['mv OneD_output.dat.0.S2_upper.%f %s'%(S2,sub_folder)], shell=True)

	return S2,d_S2


def get_upper_bound_binary_search(Cell_ID, S2):
	attempts = {}

	BCL = 1000
	d_S2 = 80
	f = open('workfile_upper', 'a')
	
	go = check_if_bilateral_cond(S2)
	lastgo = go
	print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f)

	attempts[S2] = go  # keep record
	compute = True

	up_bound = 0
	low_bound = 0
	while compute:
		if(len(list(set(list(attempts.values())))) == 1):  # if all values are the same
			S2 = S2 - d_S2*go + d_S2*(not go)

		else:
			up_bound = np.min([k for k, v in attempts.items() if v == True])
			low_bound = np.max([k for k, v in attempts.items() if v == False])
			S2 = (up_bound + low_bound) / 2.0
			if (up_bound - low_bound < 0.1):
				# print( S2, up_bound, low_bound )
				print ( "S1$%^= ", BCL, "S2= ", S2, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID, 'number attempt = ', len(attempts), file=f)
				break

		go = check_if_bilateral_cond(S2)
		attempts[S2] = go  # keep record	

		if len(attempts) > 30:
			compute = False
			S2 = np.nan
		print ( "S1= ", BCL, "S2= ", S2, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID,'number attempt = ', len(attempts),  file=f)
	return S2, up_bound, low_bound 





def get_lowwer_bound_binary_search(Cell_ID, S2):
	attempts = {}

	BCL = 1000
	d_S2 = 5
	f = open('workfile_lower', 'a')
	
	go = check_if_bilateral_block(S2)
	lastgo = go
	print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f)

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
				print ( "S1$%^= ", BCL, "S2= ", S2,"go = ", go, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID, 'number attempt = ', len(attempts), file=f)
				break

		go = check_if_bilateral_block(S2)
		attempts[S2] = go  # keep record	

		if len(attempts) > 30:
			compute = False
			S2 = np.nan
		print ( "S1= ", BCL, "S2= ", S2, "go =", go,"low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID,'number attempt = ', len(attempts),  file=f)
	return S2, up_bound, low_bound 


def get_lower_bound_binary_search_after_upper(Cell_ID, S2):
	attempts = {}

	attempts[S2] = False # based on S2 from upper limit as input

	BCL = 1000
	d_S2 = 5
	f = open('workfile_lower', 'a')

	S2 = S2 - 15
	
	go = check_if_bilateral_block(S2)
	lastgo = go
	print ( "S1= ", BCL, "S2= ", S2, "d_S2 = ", d_S2, "go = ", go, ' Cell_ID = ', Cell_ID, file=f)

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
				print ( "S1$%^= ", BCL, "S2= ", S2,"go = ", go, "low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID, 'number attempt = ', len(attempts), file=f)
				break

		go = check_if_bilateral_block(S2)
		attempts[S2] = go  # keep record	

		if len(attempts) > 30:
			compute = False
			S2 = np.nan
		print ( "S1= ", BCL, "S2= ", S2, "go =", go,"low_bound = ", low_bound, "up_bound = ", up_bound, ' Cell_ID = ', Cell_ID,'number attempt = ', len(attempts),  file=f)
	f.close()
	return S2, up_bound, low_bound 




# get_upper_bound(0, 400)
get_upper_bound_binary_search(0, 578)
get_lower_bound_binary_search_after_upper(0, 578)