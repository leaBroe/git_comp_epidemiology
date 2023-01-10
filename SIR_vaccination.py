import matplotlib.pyplot as plt


# definition of SIR model
def sir_model_vaccination(S, I, R, beta, gamma, v):
    dSdt = -beta * S * I  - v * S
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I + v * S
    return dSdt, dIdt, dRdt

# Euler's integration to solve differential equations of SIR model
def euler_integration_vaccination(S, I, R, beta, gamma, v, dt):
    dSdt, dIdt, dRdt = sir_model_vaccination(S, I, R, beta, gamma, v)
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
beta = 1.5
gamma = 0.5
v = 0.67

# Set the time step dt
dt = 0.01

# Set the number of time steps to simulate
num_steps = 365*3

# Initialize lists to store the values of S, I, and R at each time step
S_list_vaccinated = [S]
I_list_vaccinated = [I]
R_list_vaccinated = [R]

# Use Euler's method to integrate the SIR model
for i in range(num_steps):
    S, I, R = euler_integration_vaccination(S, I, R, beta, gamma, v, dt)
    S_list_vaccinated.append(S)
    I_list_vaccinated.append(I)
    R_list_vaccinated.append(R)

# Plot the results
plt.plot(S_list_vaccinated, label='Susceptible')
plt.plot(I_list_vaccinated, label='Infected')
plt.plot(R_list_vaccinated, label='Recovered / Removed')
plt.text(650, 0.25, fr'infection rate $\beta$: {beta}', fontsize=10)
plt.text(650, 0.2, fr'recovery rate $\gamma$: {gamma}', fontsize=10)
plt.title('SIR model with vaccination rate ' + str(v))
plt.xlabel('Time [days]')
plt.ylabel('Prevalence of the disease')
plt.legend()
plt.savefig(f'/Users/leabroennimann/Desktop/computational_epidemiology_essay/essay_tex/vac_rate_{v}.png')
plt.show()



