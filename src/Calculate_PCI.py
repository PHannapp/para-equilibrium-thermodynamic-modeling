import numpy as np
from scipy.optimize import minimize_scalar, brentq
from CONSTANTS import *
import database

gm, x_H = [], []


# defining Gibbs function of gas phase
def gm_gas(T, p):
    return (
        database.GHSERH2(T)
        + T * R_CONST * np.log(0.0000098692327 * p)
        - 4.29e-06 * database.FF1(p)
        - 6.35e-06 * database.FF2(p)
        - 4.25e-06 * database.FF3(p)
        + 1.5e-06 * database.FF4(p)
        + 1.63e-06 * database.FF5(p)
        + 2.479e-06 * p
        + 198428
    ) / 2


# Create function for searching IP
def search_inflection_point():
    # Calculate approximate first derivative
    dy = np.diff(gm)
    # Calculate second derivative (difference of differences)
    d2y = np.diff(dy)
    # Find points where the sign of the second derivative changes
    inflection_points_indices = np.where(np.diff(np.sign(d2y)))[0] + 1
    # Check and remove the first index if it's an inflection point
    if inflection_points_indices.size > 0 and inflection_points_indices[0] == 0:
        inflection_points_indices = inflection_points_indices[1:]
    # Check and remove the last index if it's an inflection point
    if (
        inflection_points_indices.size > 0
        and inflection_points_indices[-1] == len(gm) - 1
    ):
        inflection_points_indices = inflection_points_indices[:-1]
    return inflection_points_indices


# calculates the y value of a fully defined tangent
def tangent_line(indice, slope, intercept):
    return slope * x_H[indice] + intercept


# calculate the minimium offset from slope, intercept, and min and max boundaries
def calculate_tangent_offset(slope, intercept, boundary_min, boundary_max):
    # Iterate over the range of indices and find the maximum negative difference
    max_negative_difference = float("-inf")
    max_negative_indice = 0
    for indice in range(boundary_min, boundary_max + 1):
        diff = (
            tangent_line(indice, slope, intercept) - gm[indice]
        )  # Calculate the difference between tangent and function values
        if diff > max_negative_difference:
            max_negative_difference = diff
            max_negative_indice = indice
    return max_negative_indice, abs(
        max_negative_difference
    )  # Return the maximum offset and its x-coordinate


# make the minimization function with fixed x boundaries
def make_calculate_offset_given_point(indice, boundary_min, boundary_max):
    def calculate_offset_given_point(slope):
        intercept = (
            gm[indice] - slope * x_H[indice]
        )  # Calculate the intercept for the given slope
        _, max_offset = calculate_tangent_offset(
            slope, intercept, boundary_min, boundary_max
        )  # Calculate the maximum offset for the given slope and intercept
        return max_offset  # Return the absolute maximum offset

    return calculate_offset_given_point


# create function for optimizing the slope to tangent the ab5 curve with fixed point x,y
def optimize_slope_given_point(indice, boundary_min, boundary_max, bounds_slope):
    calculate_offset_given_point = make_calculate_offset_given_point(
        indice, boundary_min, boundary_max
    )
    result = minimize_scalar(
        calculate_offset_given_point, bounds=bounds_slope, method="bounded"
    )
    optimal_slope = result.x  # The optimal slope that minimizes the maximum offset
    optimal_intercept = (
        gm[indice] - optimal_slope * x_H[indice]
    )  # Calculate the intercept for the optimal slope
    result_x = calculate_tangent_offset(
        optimal_slope, optimal_intercept, boundary_min, boundary_max
    )
    return result_x[0], optimal_slope, optimal_intercept


# create function for searching common tangents below
def search_common_tangents(IPs, min=0, max=1, CTs=None):
    # print(min, max, CTs)
    if CTs is None:  # if first call --> initialize CT
        CTs = []
    if max > len(IPs) - 1:  # if arrived at last IP -> end recursive function
        return CTs
    if min == 0:
        boundary_min = 0
    else:
        boundary_min = IPs[min - 1]
    if max == len(IPs) - 1:
        boundary_max = len(gm) - 1
    else:
        boundary_max = IPs[max + 1]
    x1, x2 = IPs[min], IPs[max]
    # make line through IPs[min] and IPs[max] --> get slope
    # set IP[max] and optimize slope in boundaries IP[min] and IP[min-1]
    while True:
        initial_slope = (gm[x2] - gm[x1]) / (x_H[x2] - x_H[x1])
        if initial_slope > SLOPE_MAX - ACC:
            CTs = []
            return CTs
        new_x1, optimal_slope, optimal_intercept = optimize_slope_given_point(
            x2, boundary_min, x1, (initial_slope, SLOPE_MAX)
        )
        initial_slope = (gm[x2] - gm[new_x1]) / (x_H[x2] - x_H[new_x1])
        new_x2, optimal_slope, optimal_intercept = optimize_slope_given_point(
            new_x1, x2, boundary_max, (SLOPE_MIN, initial_slope)
        )
        if abs(new_x1 - x1) < ACC and abs(new_x2 - x2) < ACC:
            break
        x1, x2 = new_x1, new_x2
    if x1 < 1 or abs(x2 - BOUND_MAX) < 1:
        return CTs
    # look left of x1 if it cuts
    for j in range(SLICES):
        x1_part = int(x1 / (SLICES - 1) * j)
        if (
            gm[x1_part] - tangent_line(x1_part, optimal_slope, optimal_intercept) < -1
        ):  # or max_offset >
            # if yes --> search_common_tangents(IPs, min+2, max+2, firsttry=False, lasttry=False, CTs=None)
            CTs = search_common_tangents(IPs, max + 1, max + 2, CTs)
            return CTs
    # look right of x2 if it cuts (x3)
    for j in range(SLICES):
        x2_part = int((BOUND_MAX - x2) / (SLICES - 1) * (j + 1) + x2)
        if (
            gm[x2_part] - tangent_line(x2_part, optimal_slope, optimal_intercept) < -1
        ):  # or max_offset >
            # if yes --> take point left of x3:
            for m in range(len(IPs)):
                if IPs[m] > x2_part:
                    new_x2 = m - 1
                    # if its indice difference to x2 is uneven --> take next one
                    if (new_x2 - max) & 1:  #
                        new_x2 = new_x2 + 1
                    CTs = search_common_tangents(IPs, min, new_x2, CTs)
                    return CTs
            # it cuts right of the last
            return CTs
    # if everywhere no --> append to CTs AND search_common_tangents(IPs, min+2, max+2, firsttry=False, lasttry=False, CTs=None) AND return CTs
    CTs.append([x1, x2])
    CTs = search_common_tangents(IPs, max + 1, max + 2, CTs)
    return CTs


# finds p of gm_gas of corresponding Gibbs energy
def find_p(T, y):
    lower_bound = 1e-9  # Adjust this as needed.
    upper_bound = 1e19  # Adjust this as needed.
    try:
        return brentq(
            lambda p: gm_gas(T, p) - y, lower_bound, upper_bound, maxiter=10000
        )
    except ValueError:
        return upper_bound


# finds plateau information
def plateau(T):
    IPs = search_inflection_point()
    CTs = search_common_tangents(IPs)
    # calculate plateau pressure: one point plus tangent slope --> point at x = 1 --> gm_gas
    plateau_slopes = []
    plateau_intercepts = []
    plateau_pressures = []
    for i in range(len(CTs)):
        plateau_slopes.append(
            (gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]])
        )
        plateau_intercepts.append(gm[CTs[i][0]] - plateau_slopes[i] * x_H[CTs[i][0]])
        g_gas = plateau_slopes[i] * G_X + plateau_intercepts[i]
        p = find_p(T, g_gas)
        plateau_pressures.append(p)
    return CTs, plateau_slopes, plateau_intercepts, plateau_pressures


# make the minimization function with fixed x boundaries
def make_calculate_offset_given_point_gas(x, y, boundary_min, boundary_max):
    def calculate_offset_given_point_gas(slope):
        intercept = y - slope * x  # Calculate the intercept for the given slope
        x1, max_offset = calculate_tangent_offset(
            slope, intercept, boundary_min, boundary_max
        )  # Calculate the maximum offset for the given slope and intercept
        return max_offset  # Return the absolute maximum offset

    return calculate_offset_given_point_gas


# create function for optimizing the slope to tangent the ab5 curve with fixed point x,y
def optimize_slope_given_point_gas(x, y, boundary_min, boundary_max, bounds_slope):
    calculate_offset_given_point = make_calculate_offset_given_point_gas(
        x, y, boundary_min, boundary_max
    )
    result = minimize_scalar(
        calculate_offset_given_point, bounds=bounds_slope, method="bounded"
    )
    optimal_slope = result.x  # The optimal slope that minimizes the maximum offset
    optimal_intercept = (
        y - optimal_slope * x
    )  # Calculate the intercept for the optimal slope
    result_x = calculate_tangent_offset(
        optimal_slope, optimal_intercept, boundary_min, boundary_max
    )
    return result_x[0], optimal_slope, optimal_intercept


# optimizer for calculate_x_from_p
def optimize_x_given_p_gas(p, MINIMUM_X, MAXIMUM_X, MINIMUM_SLOPE, MAXIMUM_SLOPE, T):
    x, slope, intercept = optimize_slope_given_point_gas(
        G_X, gm_gas(T, p), MINIMUM_X, MAXIMUM_X, (MINIMUM_SLOPE, MAXIMUM_SLOPE)
    )
    return x, slope, intercept


# calculates the x value to a corresponding p value from the pop file
def calculate_x_from_p(p, plateau_pressures, CTs, T):
    # find where p is out of plateau pressures
    # calculate x from p by SLOPE_MIN or SLOPE_MAX
    # also X_Min and X_max from CTs
    for i in range(len(plateau_pressures)):
        if plateau_pressures[i] > p:  # if p is below the ith plateau-pressure
            if i != 0:  # if it is not the first
                MINIMUM_X, MAXIMUM_X = CTs[i][0], CTs[i][1]
                # MINIMUM_PRESSURE, MAXIMUM_PRESSURE = plateau_pressures[i-1], plateau_pressures[i]
                MINIMUM_SLOPE, MAXIMUM_SLOPE = (
                    (gm[CTs[i - 1][1]] - gm[CTs[i - 1][0]])
                    / (x_H[CTs[i - 1][1]] - x_H[CTs[i - 1][0]]),
                    (gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]),
                )
            else:  # otherwise the lower bounds are zero
                MINIMUM_X, MAXIMUM_X = 0, CTs[i][0]
                # MINIMUM_PRESSURE, MAXIMUM_PRESSURE = 0, plateau_pressures[i]
                MINIMUM_SLOPE, MAXIMUM_SLOPE = (
                    SLOPE_MIN,
                    (gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]),
                )
        elif i == len(plateau_pressures) - 1:  # if no plateau-pressure is above:
            MINIMUM_X, MAXIMUM_X = CTs[i][1], len(x_H) - 1
            MINIMUM_PRESSURE, MAXIMUM_PRESSURE = plateau_pressures[i], P_MAX
            MINIMUM_SLOPE, MAXIMUM_SLOPE = (
                (gm[CTs[i][1]] - gm[CTs[i][0]]) / (x_H[CTs[i][1]] - x_H[CTs[i][0]]),
                SLOPE_MAX,
            )
    return optimize_x_given_p_gas(
        p, MINIMUM_X, MAXIMUM_X, MINIMUM_SLOPE, MAXIMUM_SLOPE, T
    )
