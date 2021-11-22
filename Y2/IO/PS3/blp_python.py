# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 13:33:40 2021

@author: micha
Code draws strongly from the pyblp docs written by Jeff Gortmaker:
    https://pyblp.readthedocs.io/en/stable/_notebooks/tutorial/nevo.html
I actually know Jeff (former coworker) and he gave me lots of advice when I was applying to
grad schools. Small world!
"""
import pyblp
import numpy as np
import pandas as pd

pyblp.options.digits = 2
pyblp.options.verbose = False


product_data = pd.read_csv(pyblp.data.NEVO_PRODUCTS_LOCATION)

logit_formulation = pyblp.Formulation('prices', absorb='C(product_ids)')
problem = pyblp.Problem(logit_formulation, product_data)

logit_results = problem.solve()
#logit_results

def solve_nl(df):
    groups = df.groupby(['market_ids', 'nesting_ids'])
    df['demand_instruments20'] = groups['shares'].transform(np.size)
    nl_formulation = pyblp.Formulation('0 + prices')
    problem = pyblp.Problem(nl_formulation, df)
    return problem.solve(rho=0.7)


df1 = product_data.copy()
df1['nesting_ids'] = 1
nl_results1 = solve_nl(df1)
#nl_results1

df2 = product_data.copy()
df2['nesting_ids'] = df2['mushy']
nl_results2 = solve_nl(df2)
#nl_results2

nl_results1.beta[0] / (1 - nl_results1.rho) # result on prices
nl_results2.beta[0] / (1 - nl_results2.rho)

# nonlinear parts
X1_formulation = pyblp.Formulation('0 + prices', absorb='C(product_ids)')
X2_formulation = pyblp.Formulation('1 + prices + sugar + mushy')
product_formulations = (X1_formulation, X2_formulation)

# This part specifies integration technique.
# This version specifies MC over 50 consumers
mc_integration = pyblp.Integration('monte_carlo', size=50, specification_options={'seed': 0})
mc_problem = pyblp.Problem(product_formulations, product_data, integration=mc_integration)
bfgs = pyblp.Optimization('bfgs', {'gtol': 1e-4})
results1 = mc_problem.solve(sigma=np.ones((4, 4)), optimization=bfgs)
# Restricted Sigma to be diagonal
results3 = mc_problem.solve(sigma=np.eye(4), optimization=bfgs)

agent_data = pd.read_csv(pyblp.data.NEVO_AGENTS_LOCATION)
agent_formulation = pyblp.Formulation('0 + income + income_squared + age + child')
nevo_problem = pyblp.Problem(
    product_formulations,
    product_data,
    agent_formulation,
    agent_data
)


initial_sigma = np.diag([0.3302, 2.4526, 0.0163, 0.2441])
initial_pi = np.array([
  [ 5.4819,  0,      0.2037,  0     ],
  [15.8935, -1.2000, 0,       2.6342],
  [-0.2506,  0,      0.0511,  0     ],
  [ 1.2650,  0,     -0.8091,  0     ]
])
tighter_bfgs = pyblp.Optimization('bfgs', {'gtol': 1e-5})
nevo_results = nevo_problem.solve(
    initial_sigma,
    initial_pi,
    optimization=tighter_bfgs,
    method='1s'
)

# This nevo_results object is our solved model.



