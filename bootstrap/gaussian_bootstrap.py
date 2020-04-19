# mass bootstrap

import numpy as np
import matplotlib.pyplot as plt


# mean and standard deviation
mu_2695    = 2695.0
sigma_2695 = 2.0

mu_2770    = 2766.0
sigma_2770 = 2.0

mu_3000    = 3000.4
sigma_3000 = 0.3742

mu_3050    = 3050.2
sigma_3050 = 0.3317

mu_3066    = 3065.6
sigma_3066 = 0.4359

mu_3090    = 3090.2
sigma_3090 = 0.6557

mu_3188    = 3188.0
sigma_3188 = 13.9228

s_2695 = np.random.normal(mu_2695, sigma_2695, 1000)
s_2770 = np.random.normal(mu_2770, sigma_2770, 1000)
s_3000 = np.random.normal(mu_3000, sigma_3000, 1000)
s_3050 = np.random.normal(mu_3050, sigma_3050, 1000)
s_3066 = np.random.normal(mu_3066, sigma_3066, 1000)
s_3090 = np.random.normal(mu_3090, sigma_3090, 1000)
s_3188 = np.random.normal(mu_3188, sigma_3188, 1000)


print(mu_2695/s_2695.mean(), sigma_2695/s_2695.std())
print(mu_2770/s_2770.mean(), sigma_2770/s_2770.std())
print(mu_3000/s_3000.mean(), sigma_3000/s_3000.std())
print(mu_3050/s_3050.mean(), sigma_3050/s_3050.std())
print(mu_3066/s_3066.mean(), sigma_3066/s_3066.std())
print(mu_3090/s_3090.mean(), sigma_3090/s_3090.std())
print(mu_3188/s_3188.mean(), sigma_3188/s_3188.std())



# print(abs(mu - np.mean(s)) < 0.01)


# # Display the histogram of the samples, along with the probability density function
# count, bins, ignored = plt.hist(s, 30, density=True)
# plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *np.exp( - (bins - mu)**2 / (2 * sigma**2) ),linewidth=2, color='g')

# plt.savefig('gauss_boot.png')