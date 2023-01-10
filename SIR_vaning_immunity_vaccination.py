import matplotlib.pyplot as plt


# definition of SIR model
def sir_model(S, I, R, beta, gamma, v, N, mu):
    dSdt = -beta * S * I - v * S + mu * R
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I + v * S - mu * R
    return dSdt, dIdt, dRdt

# Euler's integration to solve differential equations of SIR model
def euler_integration(S, I, R, beta, gamma, v, N, dt):
    dSdt, dIdt, dRdt = sir_model(S, I, R, beta, gamma, v, N, mu)
    S += dt * dSdt
    I += dt * dIdt
    R += dt * dRdt
    return S, I, R

# initial values of S, I, and R
S = 0.9
I = 0.1
R = 0

# parameters beta, gamma, v, mu and N
# beta: infection rate
# gamma: recovery rate
# v: vaccination rate
# mu: rate at which recovered individuals become susceptible again
beta = 1.5
gamma = 0.5
v = 0.67
N = S + I + R
mu = 0.2

# Set the time step dt
dt = 0.01

# Set the number of time steps to simulate
num_steps = 3*365

# Initialize lists to store the values of S, I, and R at each time step
S_list = [S]
I_list = [I]
R_list = [R]

# Use Euler's method to integrate the SIR model
for i in range(num_steps):
    S, I, R = euler_integration(S, I, R, beta, gamma, v, N, dt)
    S_list.append(S)
    I_list.append(I)
    R_list.append(R)

# Plot the results
plt.plot(S_list, label='Susceptible')
plt.plot(I_list, label='Infected')
plt.plot(R_list, label='Recovered / Removed')
plt.title(f'SIR model with vaccination rate {v} vaning immunity')
plt.text(250, 0.87, fr'infection rate $\beta$: {beta}', fontsize=10)
plt.text(250, 0.82, fr'recovery rate $\gamma$: {gamma}', fontsize=10)
plt.text(250, 0.77, fr'vaning immunity rate $\mu$: {mu}', fontsize=10)
plt.xlabel('Time')
plt.ylabel('Prevalence')
plt.legend()
plt.savefig(f'/Users/leabroennimann/Desktop/computational_epidemiology_essay/essay_tex/vac_vaning_{mu}.png')
plt.show()