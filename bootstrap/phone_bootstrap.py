# construct a population pickups for our lab
np.random.seed(42)
pickups = np.random.randint(0,500 , size=100)
print(pickups)

# population mean
print(pickups.mean())

# population standard deviation
print(pickups.std())

# draw a sample from population
sample = np.random.choice(pickups, size=30)
print(sample)

# our first sample mean
sample_mean = sample.mean()
print(sample_mean)

# standard deiveation for this sample
sample_std = np.std(sample, ddof=1)
print(sample_std)

# estimated standard error for the sapmle mann
print(sample_std/(30 ** 0.5))


# construct the simulated sampling distribution
sample_props = []
for _ in range(100000):
    sample = np.random.choice(pickups, size=30)
    sample_props.append(sample.mean())

# the simulated mean of the sampling distribution
simulated_mean = np.mean(sample_props)

# the simulated standard deviation of the sampling distribution
simulated_std = np.std(sample_props)

# plot the simulated sampling distribution,
# under the Central Limit Theorem,
# it is expected normal
plt.hist(sample_props)


# the theorical mean and simulated mean
print(pickups.mean(), simulated_mean)

# the theorical standard error and simulated standard error
print(pickups.std()/(30 ** 0.5), simulated_std)


######################################################

sample = np.random.choice(pickups, size=30)
print(sample)

# bootstrap for mean
boot_means = []
for _ in range(10000):
    bootsample = np.random.choice(sample,size=30, replace=True)
    boot_means.append(bootsample.mean())

# simulated mean of mean
bootmean = np.mean(boot_means)

# simulated standard deviation of mean
bootmean_std = np.std(boot_means)

# simulated mean VS true mean
print(pickups.mean(), bootmean)

# the theorical standard error and simulated standard error
print(pickups.std()/(30 ** 0.5), bootmean_std)
