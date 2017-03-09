import wake
from math import sqrt, log, tan, cos, floor, ceil, pi, sin, radians
from Jensen_modular.Jensen_windfarm import distance_to_front, Ct, power
import sys


def jensen_angle(layout_x, layout_y, wind_speed, angle):
    U0 = wind_speed  # Free stream wind speed

    nt = len(layout_x)  # Number of turbines ## Length of layout list

    k = 0.04  # Decay constant
    r0 = 40.0  # Turbine rotor radius
    # angle2 = - 270.0 - angle  # To read windroses where N is 0 and E is 90.
    angle3 = angle + 180.0
    deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
    proportion = [[0.0 for _ in range(nt)] for _ in range(nt)]
    distance = [[0.0 for _ in range(2)] for _ in range(nt)]
    U = [U0 for _ in range(nt)]
    total_deficit = [0.0 for _ in range(nt)]

    for tur in range(nt):
        distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle), tur]
    distance.sort()

    for turbine in range(nt):
        for num in range(turbine):
            total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0

        total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
        U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])

        for i in range(turbine + 1, nt):

            determ = wake.determine_if_in_wake(layout_x[distance[turbine][1]],
                                                 layout_y[distance[turbine][1]],
                                                 layout_x[distance[i][1]], layout_y[distance[i][1]], k,
                                                 r0, angle3)
            proportion[distance[turbine][1]][distance[i][1]] = determ[0]

            if proportion[distance[turbine][1]][distance[i][1]] != 0.0:
                deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[turbine][1]][
                                                                           distance[i][1]] * wake.wake_deficit(
                    Ct(U[distance[turbine][1]]), k, determ[1], r0)
            else:
                deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

    # Farm efficiency
    profit = 0.0

    for l in range(nt):
        profit += power(U[l])

    for g in layout_x:
        print '%f\t' % g
    print
    # print profit
    return - profit

if __name__ == '__main__':
    with open('power3.dat', 'w') as out:
        for i in range(40, 900, 10):
            for j in range(900, 1120, 10):
                out.write('%d\t%d\t%f\n' % (i, j, jensen_angle([0., float(i), float(j)], 8.5, 180.0)))
