
def sir_model(S, I, R, beta, gamma):
    dSdt = -beta * S * I
    dIdt = beta * S * I - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

def euler_integration(S, I, R, beta, gamma, dt):
    dSdt, dIdt, dRdt = sir_model(S, I, R, beta, gamma)
    S += dt * dSdt
    I += dt * dIdt
    R += dt * dRdt
    return S, I, R

# Set the initial values of S, I, and R
S = 0.9
I = 0.1
R = 0

# Set the parameters beta, gamma,
beta = 1.5
gamma = 0.5

# Set the time step dt
dt = 0.01

# Set the number of time steps to simulate
num_steps = 365*3

# Initialize lists to store the values of S, I, and R at each time step
S_list = [S]
I_list = [I]
R_list = [R]

# Use Euler's method to integrate the SIR model
for i in range(num_steps):
    S, I, R = euler_integration(S, I, R, beta, gamma, dt)
    S_list.append(S)
    I_list.append(I)
    R_list.append(R)

# Plot the results
import matplotlib.pyplot as plt

plt.plot(S_list, label='Susceptible')
plt.plot(I_list, label='Infected')
plt.plot(R_list, label='Recovered / Removed')
plt.text(650, 0.25, fr'infection rate $\beta$: {beta}', fontsize=10)
plt.text(650, 0.2, fr'recovery rate $\gamma$: {gamma}', fontsize=10)

plt.title('Basic SIR model')
plt.xlabel('Time [days]')
plt.ylabel('Prevalence of the disease')
plt.legend()
plt.savefig('/Users/leabroennimann/Desktop/computational_epidemiology_essay/essay_tex/basic_SIR.png')
plt.show()