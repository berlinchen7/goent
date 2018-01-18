import sys, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def main():
	x_min, x_max, y_min, y_max, chi, increment, InputFilePath = float(sys.argv[1]), \
											float(sys.argv[2]), float(sys.argv[3]), \
											float(sys.argv[4]), float(sys.argv[5]), \
											float(sys.argv[6]), sys.argv[7]
	InputFile = open("outputs/" + InputFilePath + ".txt")
	data = []
	for line in InputFile:
		dataPoints = line.split()
		data.append(dataPoints)
	data = np.array(data, dtype=float)

	plot = plt.imshow(data)
	plt.pcolor(data, vmin=0, vmax=1)
	plt.xlim([x_min, (x_max-x_min)/increment])
	plt.ylim([y_min, (y_max-y_min)/increment])
	plt.yticks( (0, 10, 20, 30, 40, 50), (0, 1, 2, 3, 4, 5) )
	plt.xticks( (0, 10, 20, 30, 40, 50), (0, 1, 2, 3, 4, 5) )
	plt.title("MC_SY: Chi = " + str(chi))
	plt.colorbar()
	plt.savefig("plots/" + InputFilePath + ".png")
	plt.close()



if __name__ == "__main__":
	main()