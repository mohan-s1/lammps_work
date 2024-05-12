#/* -------------------------------------------------------------------------------------------------
#   Egorov Group
#   Mohan Shankar, mjs7eek@virginia.edu

#   maxwell_construction.py
#   "This file creates a Maxwell construction in order to find conditions for phase coexistence"
#------------------------------------------------------------------------------------------------- */

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


rho_one = np.array([1.2, 1.3,	1.4,	1.45,	1.5,	1.55,	1.6]) # range of density values for phase 1 

helm_one = np.array([-21.5724,	-22.1546,	-22.5722,	-22.7215,	-22.836,	-22.91365,	-22.96]) # range of free energy values for phase 1

rho_two = np.array([1.6,	1.7,	1.8,	1.9,	2,	2.1,	2.2]) # range of density values for phase 2

helm_two = np.array([-22.96,	-23.0469,	-23.013,	-22.8475,	-22.544,	-22.1004,	-21.5226]) # range of free energy values for phase 2

def quadratic_poly(x, a, b, c): # ax^2 + bx + c
    return  a*x**2 + b*x + c

def deriv_quadratic(x, a, b): # 2ax + b; derive of above 
    return 2*a*x + b

p1_weights, p1_cov = curve_fit(quadratic_poly, rho_one, helm_one) # returns highest order first; ax^2 + bx + c

a_p1, b_p1, c_p1 = p1_weights[0], p1_weights[1], p1_weights[2]

p2_weights, p2_cov = curve_fit(quadratic_poly, rho_two, helm_two) # returns highest order first; ax^2 + bx + c

a_p2, b_p2, c_p2 = p2_weights[0], p2_weights[1], p2_weights[2]

print("Phase One Coexistence Curve Weights: {}, Associated Liquid Errors: {}".format(p1_weights, np.sqrt(np.diag(p1_cov))))
print("Phase Two Coexistence Curve Weights: {}, Associated BCC Errors: {}".format(p2_weights, np.sqrt(np.diag(p2_cov))))

#-------------------------------------------------------------------------------- */
# Define range of r/sigma values to plot phase 1 and 2 curves

x_p1 = np.linspace(1.2, 2.2, 100) 

x_p2 = np.linspace(1.6, 2.2, 100)
#-------------------------------------------------------------------------------- */
# create the fits over the domain specified above
fx_p1 = quadratic_poly(x_p1, a_p1, b_p1, c_p1) 

fx_p1_deriv = deriv_quadratic(x_p1, a_p1, b_p1)
#-------------------------------------------------------------------------------- */
# create the fits over the domain specified above
fx_p2 = quadratic_poly(x_p2,  a_p2, b_p2, c_p2)

fx_p2_deriv = deriv_quadratic(x_p2,  a_p2, b_p2)
#-------------------------------------------------------------------------------- */
# plot data points and fit
fig, ax = plt.subplots()

ax.plot(x_p1, fx_p1, color = 'k')
ax.scatter(rho_one, helm_one, color = 'black', label = 'Phase 1')

ax.plot(x_p2, fx_p2, color = 'RoyalBlue')
ax.scatter(rho_two, helm_two, color = 'RoyalBlue', label = 'Phase 2')

ax.set_xlabel(r'$ r / \sigma $')
ax.set_ylabel(r'$\sigma^3 \beta F/V - 26 \rho \sigma^3$')

ax.legend()
ax.set_title(r"Free Energy Curves for Different Phases, $\frac{K_B T}{\varepsilon} \, = \, 0.2$")
plt.show();
#-------------------------------------------------------------------------------- */
# Solve system of equations to find common tangent 

def equations(p): # establish two equations that we find x1 and x2 for 
    '''
    change the a, b, c values here if you want to check co-existence between different phases
    '''
    x1, x2 = p
    E1 = deriv_quadratic(x1, a_p1, b_p1) - deriv_quadratic(x2, a_p2, b_p2) # equation 1
    E2 = ( (quadratic_poly(x1, a_p1, b_p1, c_p1) - quadratic_poly(x2, a_p2, b_p2, c_p2)) / (x1 - x2) ) - deriv_quadratic(x1, a_p1, b_p1) # equation 2
    return (E1, E2)

from scipy.optimize import least_squares

lb = (1.2, 1.5)   # lower bounds on x1, x2
ub = (2.2, 2.2)   # upper bounds on x1, x2

result = least_squares(equations, [1.2, 1.5], ftol=1e-12, xtol=1e-12, gtol=1e-12, bounds=(lb, ub)) # regress equations to find optimal solution

print(result) # results of least squares optimization 

#-------------------------------------------------------------------------------- */
# Print out relevant data

point_one = result.x[0] # x-value 1 (r/sigma) for first co-existence point
point_two = result.x[1] # x-value 2 (r/sigma) for second co-existence point

print("Coexistence Point 1 (r/sigma): {}, Coexistence Point 2 (r/sigma): {}".format(point_one, point_two)) # MOST IMPORTANT THINGS; R/SIGMA WHERE COEXISTENCE OCCURS 

slope_common_tangent = deriv_quadratic(point_one, a_p1, b_p1)

y1 = quadratic_poly(point_one, a_p1, b_p1, c_p1)
y2 = quadratic_poly(point_two, a_p2, b_p2, c_p2)

print("y-value Tangent on Curve One:", y1)
print("y-value Tangent on Curve Two:", y2)

def line(x, m, b):
    return m*x + b

tangent_weights, tangent_cov = curve_fit(line, (point_one, point_two), (y1, y2)) # returns tangent line, y = mx + b

print("Coexistence Pressure: {}, Coexistence Chemical Potential: {}".format(tangent_weights[0], tangent_weights[1]))

m_tan, b_tan = tangent_weights[0], tangent_weights[1]

final_x = np.linspace(1.2, 2.2, 100) # change this to match the denisty range as needed 

real_tangent = line(final_x, m_tan, b_tan)

#-------------------------------------------------------------------------------- */
# Plot Results
fig, ax = plt.subplots()

ax.plot(x_p1, fx_p1, color = 'k') # plot fit for liq curve
ax.scatter(rho_one, helm_one, color = 'black', label = 'Phase 1') # plot points for liq curve

ax.plot(x_p2, fx_p2, color = 'RoyalBlue') # plot fit for bcc curve
ax.scatter(rho_two, helm_two, color = 'RoyalBlue', label = 'Phase 2') # plot points for bcc curve

ax.plot(final_x, real_tangent, '--', color = 'Pink') # plot tangent line

ax.set_xlabel(r'$ r / \sigma $')

plt.legend()
plt.show();