import math
import numpy as np

def bound_check(tau, N):
	return max(min(tau, N), 1)

def greedy_tau(theta, N):
	return math.ceil(N - np.log(theta)/np.log(1-theta))

def optimal_tau(theta, N):
	return round(N - math.log(theta ** 2/((theta - 1)*math.log(1 - theta)))/math.log(1 - theta))

def get_expec_C(theta, N, tau):
	if tau == 1:
		return N - (N-1) * theta
	else:
		return N - (N-2)*theta + (tau - 2)*theta**2 - (1 - theta)**(N - tau + 1)

def get_sq_expec_C(theta, N, tau):
	return get_expec_C(theta, N, tau)**2

def polynomial_part(theta, N, T):
	return N**2 + (1 - 2*N)*(N - 2)*theta + ((N-2)*(N-3) + (2*N-3)*(T - 2))*(theta**2) + 2 *((3*T - 7) - (T - 2)*(N - 2))*(theta**3) + (T - 3)*(T - 4)*(theta**4)

def exponential_part(theta, N, T):
	return ((1 - 2*N) + 2*(T - 3) * theta + 2*(3 - T)* (theta**2))*((1 - theta)**(N - T + 1))

def get_expec_C_sq(theta, N, T):
	if T == 1:
		return N**2 + (N - 1)*(N - 2)*theta**2 + (1 - 2*N)*(N - 1)*theta 
	elif T == 2:
		return N**2 + (N-2)*(N-3)*theta**2 + (1 - 2*N)*((N-2)*theta + (1 - theta)**(N-1))
	else:
		return polynomial_part(theta, N, T) + exponential_part(theta, N, T)

def get_var_C(theta, N, tau):
	return get_expec_C_sq(theta, N, tau) - get_sq_expec_C(theta, N, tau)

def test(theta, N, tau):
	return N**2 - 2*N + 1 + (N-3)*(N-4) * (theta ** 2) + (N-3)*theta + 2*((N-3)*(theta**3) + ((N-2)*(N-3)/2 - (N - 3))*(theta**4)) + (N - 2)*(theta**2) + 2 * ((1 - N)*(N-3)*theta + (N-1)*(N-2)*(theta**2) - (2*N - 6)*theta**2 - ((N-3)*(N-2) - (2*N - 6))*(theta**3))

def test_2(theta, N, tau):
	return (N - 1 - (N-3)*theta + (N-2) * (theta ** 2))**2

def main():
	N=52
	thetas = [0.01 + 0.01*x for x in range(0, 99)]
	for theta in thetas:
		g_tau = greedy_tau(theta, N)
		print '$' + str(theta) + '$ & $ ' + str(get_expec_C(theta, N, g_tau)) + '$ & $' + str(get_var_C(theta, N, g_tau)) +  '$\\\\'
	correct = 0
	for theta in range(1, 1000):
		theta = float(theta)/1000
		g_tau = bound_check(greedy_tau(theta, N), N)
		o_tau = bound_check(optimal_tau(theta, N), N)
		top_expec = 0
		stopping_card = 0
		for i in range(1, N + 1):
			expec = get_expec_C(theta, N, i)
			if expec > top_expec:
				top_expec = expec
				stopping_card = i
		if stopping_card != g_tau:
			print 'oh no', theta
		else:
			correct += 1
	print correct
if __name__ == "__main__":
    main()

