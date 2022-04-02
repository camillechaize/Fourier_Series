import scipy.integrate as integral
import matplotlib.pyplot as plt
import numpy as np
import math

# Function To Sample
period = 1
pulsation = 2 * math.pi / period
d = 0.7  # cyclic

# Precision
n = 50  # decomposition precision
p = 500  # graph precision
number_periods = 2


# Square function
def periodic_function(t):
    if (t % period) / period < d:
        return 1
    else:
        return 0


def periodic_function_2(t):
    if (t % period) / period < d:
        return math.cos(pulsation*t)
    else:
        return -math.cos(pulsation*t)


def coefficients_calculation(function_to_sample):
    m_0 = integral.quad(lambda x: function_to_sample(x), 0, period)[0]
    a_n = np.empty(n)
    b_n = np.empty(n)
    for k in range(1, n):
        a_n[k] = (2/period)*(integral.quad(lambda x: function_to_sample(x) * math.cos(k * pulsation * x), 0, period))[0]
        b_n[k] = (2/period)*(integral.quad(lambda x: function_to_sample(x) * math.sin(k * pulsation * x), 0, period))[0]

    return m_0, a_n, b_n


def fourier_series_decomposition(function_to_sample):
    x_n = coefficients_calculation(function_to_sample)

    # Evaluate
    x_values = [number_periods * period * t / p for t in range(0, p, 1)]
    y_values = [0.0] * p

    for i, x in enumerate(x_values):
        y = x_n[0]
        for k in range(x_n[1].shape[0]):
            y += x_n[1][k] * math.cos(k * pulsation * x) + x_n[2][k] * math.sin(k * pulsation * x)
        y_values[i] = y

    return x_values, y_values


def plot_estimation(function_to_sample):
    x_values, y_values = fourier_series_decomposition(function_to_sample)
    x_ref, y_ref = [number_periods * period * t / p for t in range(0, p, 1)], [function_to_sample(number_periods * period * t / p) for t in range(0, p, 1)]
    plt.plot(x_values, y_values, "r")
    plt.plot(x_ref, y_ref, "g")
    plt.show()


plot_estimation(periodic_function)
