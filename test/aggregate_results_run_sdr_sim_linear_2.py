#!./python277/bin/python

import os
import shutil
import os.path
import platform
import sys
import time
import numpy
import csv
import shlex
from scipy import sqrt, pi

print >> sys.stderr, __doc__
print "Version :", platform.python_version()
print "Program :", sys.executable
print 'Script  :', os.path.abspath(__file__)

numpy.set_printoptions(precision=3, suppress = True)


## NEED TO CHANGE THESE VALUES
ngrid = 24
total_runs_per_sim = 100
n_sims = ngrid * total_runs_per_sim
move_files = True
main_folder = 'run_sdr_sim_linear_1'

completed_per_grid = numpy.zeros(ngrid,)


if not os.path.exists(main_folder + "/" + "results"):
	os.mkdir(main_folder + "/" + "results")
else:
	# make sure that these directories exist
	dir_src = main_folder + "/" + "results"
	dir_dst = main_folder
	for file in os.listdir(dir_src):
	    #print file  # testing
	    src_file = os.path.join(dir_src, file)
	    dst_file = os.path.join(dir_dst, file)
	    shutil.move(src_file, dst_file)

totalFailed = 0
results_set_up = False

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def converter(x):
    if x == 'NA':
        return numpy.nan
    else:
        return float(x)

for s in range(n_sims):

	#filename_sim_params    = main_folder + "/" + str(s) + "_sim_params.csv"
	filename_results_angle = main_folder + "/" + str(s) + "_results_angle.csv"
	filename_results_corr  = main_folder + "/" + str(s) + "_results_correlation.csv"

	#filename_sim_params_mv        = main_folder + "/results/" + str(s) + "_sim_params.csv"
	filename_results_angle_mv = main_folder + "/results/" + str(s) + "_results_angle.csv"
	filename_results_corr_mv  = main_folder + "/results/" + str(s) + "_results_correlation.csv"

	notFailed = os.path.exists(filename_results_corr) or os.path.exists(filename_results_corr_mv)

	if notFailed:

		if move_files:
			shutil.move(filename_results_angle, filename_results_angle_mv)
			shutil.move(filename_results_corr, filename_results_corr_mv)

		if not results_set_up:
			results_set_up = True

			if move_files:
				reader = csv.reader(open(filename_results_corr_mv, 'rb'), delimiter=',')
			else:
				reader = csv.reader(open(filename_results_corr, 'rb'), delimiter=',')
			headers = reader.next()
			headers = headers[1:]
			ncols3 = len(headers)
			print(headers)
			if move_files:
				currowResultsCorr  = numpy.genfromtxt(open(filename_results_corr_mv,"rb"), delimiter=",",dtype=None, usecols=range(1,ncols3+1))
			else:
				currowResultsCorr  = numpy.genfromtxt(open(filename_results_corr,"rb"), delimiter=",",dtype=None, usecols=range(1,ncols3+1))
			types = [x.dtype for x in currowResultsCorr[1]]
			print(types)
			results_corr  = numpy.zeros((0, len(headers)),)
			results_angle = numpy.zeros((0, len(headers)),)
			#results_corr[0,:]  = headers
			#results_angle[0,:] = headers

		if move_files:
			#print filename_results_mv
			currowResultsCorr  = numpy.genfromtxt(open(filename_results_corr_mv,"rb"), delimiter=",",dtype=None, usecols=range(1,ncols3+1)) # , usecols=range(1,ncols3+1)
			currowResultsAngle = numpy.genfromtxt(open(filename_results_angle_mv,"rb"),delimiter=",",dtype=None, usecols=range(1,ncols3+1))
		else:
			#print filename_results_mv
			currowResultsCorr  = numpy.genfromtxt(open(filename_results_corr,"rb"), delimiter=",",dtype=None, usecols=range(1,ncols3+1))
			currowResultsAngle = numpy.genfromtxt(open(filename_results_angle,"rb"),delimiter=",",dtype=None, usecols=range(1,ncols3+1))
		#values2enter = numpy.concatenate([gridCur,currowResults])
		#results_corr[s,:]  = currowResultsCorr[1]
		crres  = [x.replace('"', '').strip() for x in currowResultsCorr[1]]
		crresa = [x.replace('"', '').strip() for x in currowResultsAngle[1]]

		results_corr = numpy.vstack((results_corr, crres))
		#results_angle[s,:] = currowResultsAngle[1]
		results_angle = numpy.vstack((results_angle, crresa))

	else:
		totalFailed = totalFailed + 1
		if totalFailed == 1:
			whichFailed = numpy.array([s])
		else:
			numpy.append(whichFailed, s)


if totalFailed > 0:
	if not os.path.exists(main_folder + "/" + "failedrun"):
		os.mkdir(main_folder + "/" + "failedrun")

print "total failed:"
print totalFailed



resultsheader=",".join( headers )+'\n'
#numpy.savetxt("results.csv", results, delimiter=',', header=resultsheader)

curDate = time.strftime("%m_%d_%Y")

#results_corr.tofile(main_folder + "/" + 'results_correlation_' + curDate + '.csv', sep = ",")

# open a file for writing.
csv_out = open(main_folder + "/" + 'results_correlation_' + curDate + '.csv', 'wb')
# create the csv writer object.
mywriter = csv.writer(csv_out, delimiter = ",")
mywriter.writerows([headers])
# all rows at once.
mywriter.writerows(results_corr)

# always make sure that you close the file.
# otherwise you might find that it is empty.
csv_out.close()

#with open(main_folder + "/" + 'results_correlation_' + curDate + '.csv', 'wb') as f:
#	#f.write(resultsheader)
#	for row in results_corr:
#		f.write(repr(row)+' ')

#with open(main_folder + "/" + 'results_angle_' + curDate + '.csv', 'wb') as f:
#	f.write(resultsheader)
#	numpy.savetxt(f, results_angle, delimiter=",")


# open a file for writing.
csv_out = open(main_folder + "/" + 'results_angle_' + curDate + '.csv', 'wb')
# create the csv writer object.
mywriter = csv.writer(csv_out, delimiter = ",")
mywriter.writerows([headers])
# all rows at once.
mywriter.writerows(results_angle)

# always make sure that you close the file.
# otherwise you might find that it is empty.
csv_out.close()


print "All done."
sys.exit(0)
