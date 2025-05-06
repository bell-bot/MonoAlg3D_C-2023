import sys
import numpy as np
import matplotlib.pyplot as plt

def read_ecg_readings (input_file):
	data = np.genfromtxt(input_file)
	nlin, ncol = np.shape(data)
	timesteps = data[:,0]
	currents = data[:,1:]
	num_leads = ncol-1
	
	return timesteps, currents, num_leads
	
def plot_ecg_readings (t, currents, nleads):
	fig, axs = plt.subplots(nleads, figsize=(7, 5*nleads))
	for i in range(nleads):
		ax = axs[i]
		ax.plot(t, currents[:,i], label="lead_%d" % (i), c="black", linewidth=3.0)
		ax.set_xlabel("t (ms)",fontsize=15)
		ax.set_ylabel("Current (mA)",fontsize=15)
		ax.set_title("ECG reading - Lead %d" % (i),fontsize=14)
		ax.legend(loc=0,fontsize=14)
		#plt.show(block=True)
	plt.savefig("ecg_leads.pdf")

def main():
	
	if len(sys.argv) != 2:
		print("-------------------------------------------------------------------------")
		print("Usage:> python %s <input_file>" % sys.argv[0])
		print("-------------------------------------------------------------------------")
		print("<input_file> = Input file with the ECG reading from each timestep")
		print("-------------------------------------------------------------------------")
		return 1

	input_file = sys.argv[1]

	t, currents, num_leads = read_ecg_readings(input_file)

	plot_ecg_readings(t,currents,num_leads)

if __name__ == "__main__":
	main()
