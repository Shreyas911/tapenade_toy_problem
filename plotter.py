import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("results_forward_run.txt",skiprows = 2)

length = data.shape[0]
plt.plot(data[:,0],data[:,2], label="glacier profile")
plt.plot(data[:,0],data[:,3], label="inclined base")
plt.legend()
plt.savefig("forward_run_results.png")

plt.close()

plt.plot(data[:,0],data[:,1], label="glacier profile")

data = np.loadtxt("results.txt",skiprows = 2)

plt.plot(data[:,0],data[:,1]/100, label="0.01*(Adjoint sensitivities)")
plt.legend()
plt.savefig("combined_results.png")

