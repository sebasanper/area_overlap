__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Jensen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
import wake
from math import sqrt, log, tan, cos, floor, ceil, pi, sin, radians
import time
import wake_geometry as wake_larsen
from eddy_viscosity_integrate import ainslie


def power(U0):
    # return power_bladed(U0)
    return power_v90(U0)


def Ct(U0):
    # return ct_bladed(U0)
    return Ct_v90(U0)


def power_v90(U0):
    if U0 < 4.0:
        return 0.0
    elif U0 <= 25.0:
        return 3.234808e-4 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 ** 5.0 - 30.3162345004 * U0 ** 4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
    else:
        return 0.0


def Ct_v90(U0):
    if U0 < 4.0:
        return 0.1
    elif U0 <= 25.0:
        return 7.3139922126945e-7 * U0 ** 6.0 - 6.68905596915255e-5 * U0 ** 5.0 + 2.3937885e-3 * U0 ** 4.0 + - 0.0420283143 * U0 ** 3.0 + 0.3716111285 * U0 ** 2.0 - 1.5686969749 * U0 + 3.2991094727
    else:
        return 0.0


def power_bladed(U0):
    if U0 < 4.0:
        return 0.0
    else:
        return 0.5 * 61.0 ** 2.0 * pi * 1.225 * 0.485 * U0 ** 3.0


def ct_bladed(U0):
    if U0 < 4.0:
        return 0.1
    else:
        return 0.781


def ct_table_LLT(U0):
    v = U0
    if v == 7: return 0.977
    if v == 8: return 0.943
    if v == 9: return 0.899
    if v == 10: return 0.852
    if v == 11: return 0.804


def ct_LLT(U0):
    if ceil(U0) == floor(U0):
        return ct_table_LLT(U0)
    else:
        return interpolate(floor(U0), ct_table_LLT(floor(U0)), ceil(U0), ct_table_LLT(ceil(U0)), U0)


def power_table_LLT(U0):
    v = U0
    if v == 7: return 970.0
    if v == 8: return 1780.0
    if v == 9: return 2770.0
    if v == 10: return 3910.0
    if v == 11: return 5190.0


def power_LLT(U0):
    if ceil(U0) == floor(U0):
        return power_table_LLT(U0)
    else:
        return interpolate(floor(U0), ct_table_LLT(floor(U0)), ceil(U0), ct_table_LLT(ceil(U0)), U0)


def interpolate(minx, miny, maxx, maxy, valx):
    print maxx, minx
    return miny + (maxy - miny) * ((valx - minx) / (maxx - minx))


def distance_to_front(x, y, theta):
    theta = radians(theta)
    return abs(x + tan(theta) * y - 10000000000.0 / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)


def determine_front(wind_angle, x_t1, y_t1, x_t2, y_t2):
    wind_angle = radians(wind_angle)
    a = (x_t2 - x_t1) * cos(wind_angle) + (y_t2 - y_t1) * sin(wind_angle)
    if a > 0.0:
        return a
    else:
        return 0.0


class Layout:
    def __init__(self, layout_file):

        self.layout = open(layout_file, 'r')
        self.layout_x = []
        self.layout_y = []
        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0
        # Turbine specification
        self.radius = 63.0
        self.hub_height = 90.0

        self.wind_speed = []
        self.wind_direction = []
        self.wind_frequency = []

        for line in self.layout:
            columns = line.split()
            self.layout_x.append(float(columns[0]))
            self.layout_y.append(float(columns[1]))

        self.nt = len(self.layout_x)

    def read_windrose(self, windrose_file):

        windrose = open(windrose_file, 'r')
        for line in windrose:
            columns = line.split()
            self.wind_direction.append(float(columns[0]))
            self.wind_speed.append(float(columns[1]))
            self.wind_frequency.append(float(columns[2]))

    def jensen_angle(self, wind_speed, angle):
        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0
        U0 = wind_speed  # Free stream wind speed

        nt = len(self.layout_y)  # Number of turbines ## Length of layout list

        k = 0.04  # Decay constant
        r0 = self.radius  # Turbine rotor radius
        # angle2 = - 270.0 - angle  # To read windroses where N is 0 and E is 90.
        angle3 = angle + 180.0
        deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
        proportion = [[0.0 for _ in range(nt)] for _ in range(nt)]
        distance = [[0.0 for _ in range(2)] for _ in range(nt)]
        self.U = [U0 for _ in range(nt)]
        total_deficit = [0.0 for _ in range(nt)]

        for tur in range(nt):
            distance[tur] = [distance_to_front(self.layout_x[tur], self.layout_y[tur], angle), tur]
        distance.sort()

        for turbine in range(nt):
            for num in range(turbine):
                total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0

            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
            self.U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])

            for i in range(turbine + 1, nt):

                determ = wake.determine_if_in_wake(self.layout_x[distance[turbine][1]],
                                                   self.layout_y[distance[turbine][1]],
                                                   self.layout_x[distance[i][1]], self.layout_y[distance[i][1]], k,
                                                   r0, angle3)
                proportion[distance[turbine][1]][distance[i][1]] = determ[0]

                if proportion[distance[turbine][1]][distance[i][1]] != 0.0:
                    deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[turbine][1]][
                                                                               distance[i][1]] * wake.wake_deficit(
                        Ct(self.U[distance[turbine][1]]), k, determ[1], r0)
                else:
                    deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

        # Farm efficiency
        self.profit = 0.0

        for l in range(nt):
            self.profit += power(self.U[l])
        self.efficiency = self.profit * 100.0 / (float(nt) * power(self.U[distance[0][1]]))  # same as using U0
        self.powers = [power(self.U[i]) for i in range(nt)]

        with open('speed_jensen.dat', 'w') as out:
            for i in range(len(self.U)):
                out.write('{0:f}\t{1:f}\n'.format(self.U[i], power(self.U[i]) / power(self.U[distance[0][1]])))

    def jensen_windrose(self):
        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0
        windrose_angle = self.wind_direction
        windrose_speed = self.wind_speed
        windrose_frequency = self.wind_frequency

        nt = self.nt  # Number of turbines ## Length of layout list

        direction = open('direction_power_jensen.dat', 'w', 1)

        for wind in range(len(windrose_angle)):
            U0 = windrose_speed[wind]  # Free stream wind speed
            k = 0.04  # Decay constant
            r0 = self.radius  # Turbine rotor radius
            angle = windrose_angle[wind]
            # angle2 = - 270.0 - angle  # To read windroses where N is 0 and E is 90
            angle3 = angle + 180.0
            deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
            proportion = [[0.0 for _ in range(nt)] for _ in range(nt)]
            distance = [[0.0 for _ in range(2)] for _ in range(nt)]
            U = [U0 for _ in range(nt)]
            total_deficit = [0.0 for _ in range(nt)]

            for tur in range(nt):
                distance[tur] = [distance_to_front(self.layout_x[tur], self.layout_y[tur], angle), tur]

            distance.sort()

            for turbine in range(nt):
                for num in range(turbine):
                    total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0

                total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
                U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])

                for i in range(turbine + 1, nt):

                    determ = wake.determine_if_in_wake(self.layout_x[distance[turbine][1]],
                                                       self.layout_y[distance[turbine][1]],
                                                       self.layout_x[distance[i][1]], self.layout_y[distance[i][1]], k,
                                                       r0, angle3)
                    proportion[distance[turbine][1]][distance[i][1]] = determ[0]

                    if proportion[distance[turbine][1]][distance[i][1]] != 0.0:
                        deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[turbine][1]][
                                                                                   distance[i][1]] * wake.wake_deficit(
                            Ct(U[distance[turbine][1]]), k, determ[1], r0)
                    else:
                        deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

            # Farm efficiency
            efficiency_proportion = [0.0 for _ in range(len(windrose_frequency))]
            profit_sum = 0.0
            for l in range(nt):
                profit_sum += power(U[l])
            profit = profit_sum
            efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))  # same as using U0
            efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
            self.summation += efficiency_proportion[wind]

            direction.write('{0:f} {1:f}\n'.format(angle, profit))
        direction.close()

    def ainslie_windrose(self):
        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0
        D = 2.0 * self.radius  # Diameter

        for x in range(len(self.layout_x)):
            self.layout_x[x] = (self.layout_x[x] / D)

        for x in range(len(self.layout_y)):
            self.layout_y[x] = (self.layout_y[x] / D)

        for x in range(len(self.layout_x)):
            self.layout_xD.append(self.layout_x[x] / D)

        for x in range(len(self.layout_y)):
            self.layout_yD.append(self.layout_y[x] / D)

        layout_xD = self.layout_xD
        layout_yD = self.layout_yD

        windrose_angle = self.wind_direction
        windrose_speed = self.wind_speed
        windrose_frequency = self.wind_frequency

        nt = self.nt
        # turb_data = open('turb1_power_jensen.dat', 'w', 1)
        direction = open('direction_power_ainslie.dat', 'w', 1)

        for wind in range(0, len(windrose_angle)):

            U0 = windrose_speed[wind]  # Free stream wind speed
            angle = windrose_angle[wind]
            angle3 = angle + 180.0
            wake_deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
            distance = [[0.0 for _ in range(2)] for _ in range(nt)]
            total_deficit = [0.0 for _ in range(nt)]
            U = [U0 for _ in range(nt)]

            for tur in range(nt):
                distance[tur] = [distance_to_front(layout_xD[tur], layout_yD[tur], angle), tur]

            distance.sort()

            for turbine in range(nt):
                for num in range(turbine):
                    total_deficit[distance[turbine][1]] += wake_deficit_matrix[distance[turbine][1]][
                                                               distance[num][1]] ** 2.0
                total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
                U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
                parallel_distance = [0.0 for x in range(nt)]
                perpendicular_distance = [0.0 for x in range(nt)]

                for i in range(turbine + 1, nt):
                    parallel_distance[distance[i][1]] = determine_front(angle3, layout_xD[distance[turbine][1]],
                                                                        layout_yD[distance[turbine][1]],
                                                                        layout_xD[distance[i][1]],
                                                                        layout_yD[distance[i][1]])
                    perpendicular_distance[distance[i][1]] = wake.crosswind_distance(radians(angle3),
                                                                                     layout_xD[distance[turbine][1]],
                                                                                     layout_yD[distance[turbine][1]],
                                                                                     layout_xD[distance[i][1]],
                                                                                     layout_yD[distance[i][1]])

                    if perpendicular_distance[distance[i][1]] <= 1.7 and parallel_distance[
                        distance[i][1]] > 0.0:  ## 1.7 gives same results as a bigger distance, many times faster.

                        wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = ainslie(
                            Ct(U[distance[turbine][1]]), U[distance[turbine][1]],
                            parallel_distance[distance[i][1]], perpendicular_distance[distance[i][1]])
                    else:
                        wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0
            # turb_data.write('{0:f} {1:f}\n'.format(angle, power(U[14])))
            # Farm efficiency
            profit = 0.0
            efficiency_proportion = [0.0 for _ in range(0, len(windrose_frequency))]
            for l in range(nt):
                profit += power(U[l])
            efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
            efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
            direction.write('{0:f} {1:f}\n'.format(angle, profit))
            self.summation += efficiency_proportion[wind]

    def ainslie_angle(self, wind_speed, angle):

        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0
        self.layout_xD = []
        self.layout_yD = []

        D = 2.0 * self.radius  # Diameter

        for x in range(len(self.layout_x)):
            self.layout_xD.append(self.layout_x[x] / D)

        for x in range(len(self.layout_y)):
            self.layout_yD.append(self.layout_y[x] / D)

        layout_xD = self.layout_xD
        layout_yD = self.layout_yD

        nt = len(self.layout_y)
        U0 = wind_speed  # Free stream wind speed
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
        distance = [[0.0 for _ in range(2)] for _ in range(nt)]
        total_deficit = [0.0 for _ in range(nt)]
        self.U = [U0 for _ in range(nt)]
        U = self.U

        for tur in range(nt):
            distance[tur] = [distance_to_front(layout_xD[tur], layout_yD[tur], angle), tur]

        distance.sort()

        for turbine in range(nt):
            for num in range(turbine):
                total_deficit[distance[turbine][1]] += wake_deficit_matrix[distance[turbine][1]][
                                                           distance[num][1]] ** 2.0
            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
            U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
            parallel_distance = [0.0 for x in range(nt)]
            perpendicular_distance = [0.0 for x in range(nt)]

            for i in range(turbine + 1, nt):
                parallel_distance[distance[i][1]] = determine_front(angle3, layout_xD[distance[turbine][1]],
                                                                    layout_yD[distance[turbine][1]],
                                                                    layout_xD[distance[i][1]],
                                                                    layout_yD[distance[i][1]])
                perpendicular_distance[distance[i][1]] = wake.crosswind_distance(radians(angle3),
                                                                                 layout_xD[distance[turbine][1]],
                                                                                 layout_yD[distance[turbine][1]],
                                                                                 layout_xD[distance[i][1]],
                                                                                 layout_yD[distance[i][1]])

                if perpendicular_distance[distance[i][1]] <= 1.7 and parallel_distance[
                    distance[i][1]] > 0.0:  ## 1.7 gives same results as a bigger distance, many times faster.

                    wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = ainslie(
                        Ct(U[distance[turbine][1]]), U[distance[turbine][1]],
                        parallel_distance[distance[i][1]], perpendicular_distance[distance[i][1]])
                else:
                    wake_deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0
        # turb_data.write('{0:f} {1:f}\n'.format(angle, power(U[14])))
        # Farm efficiency
        profit = 0.0
        for l in range(nt):
            profit += power(U[l])
        self.efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
        self.powers = [power(self.U[i]) for i in range(nt)]

        with open('speed_ainslie.dat', 'w') as out:
            for i in range(len(self.U)):
                out.write('{0:f}\t{1:f}\n'.format(self.U[i], power(self.U[i]) / power(self.U[distance[0][1]])))

    def larsen_angle(self, wind_speed, angle):

        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0

        U0 = wind_speed
        layout_y = self.layout_y
        layout_x = self.layout_x

        nt = len(layout_y)  # Number of turbines
        r0 = self.radius  # Turbine rotor radius
        D = 2.0 * r0
        A = pi * r0 ** 2.0
        H = 100.0  # Hub height
        ia = 0.08  # Ambient turbulence intensity according to vanluvanee. 8% on average

        def deff(U0):
            return D * sqrt((1.0 + sqrt(1.0 - Ct(U0))) / (2.0 * sqrt(1.0 - Ct(U0))))

        rnb = max(1.08 * D, 1.08 * D + 21.7 * D * (ia - 0.05))
        r95 = 0.5 * (rnb + min(H, rnb))

        def x0(U0):
            return 9.5 * D / ((2.0 * r95 / deff(U0)) ** 3.0 - 1.0)

        def c1(U0):
            return (deff(U0) / 2.0) ** (5.0 / 2.0) * (105.0 / 2.0 / pi) ** (- 1.0 / 2.0) * (Ct(U0) * A * x0(U0)) ** (
                - 5.0 / 6.0)  # Prandtl mixing length

        angle3 = angle + 180.0
        deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
        distance = [[0.0 for _ in range(2)] for _ in range(nt)]
        self.U = [U0 for _ in range(nt)]
        total_deficit = [0.0 for _ in range(nt)]

        for tur in range(nt):
            distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle), tur]
        distance.sort()

        for turbine in range(nt):
            for num in range(turbine):
                total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0
            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
            self.U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
            flag = [False for _ in range(nt)]
            proportion = [0.0 for _ in range(nt)]
            perpendicular_distance = [0.0 for _ in range(nt)]
            parallel_distance = [0.0 for _ in range(nt)]
            for i in range(turbine + 1, nt):
                proportion[distance[i][1]], flag[distance[i][1]], perpendicular_distance[distance[i][1]], \
                parallel_distance[distance[i][1]] = wake_larsen.determine_if_in_wake_larsen(
                    layout_x[distance[turbine][1]],
                    layout_y[distance[turbine][1]],
                    layout_x[distance[i][1]],
                    layout_y[distance[i][1]], A,
                    c1(self.U[distance[turbine][1]]),
                    Ct(self.U[distance[turbine][1]]), angle3,
                    r0, x0(self.U[distance[turbine][1]]))
                if parallel_distance[
                    distance[i][1]] > 0.0:  ## Add if proportion is 0, skip operation and deficit_matrix == 0.
                    if proportion[distance[i][1]] != 0.0:
                        deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[
                                                                                   distance[i][
                                                                                       1]] * wake_larsen.wake_deficit(
                            self.U[distance[turbine][1]], Ct(self.U[distance[turbine][1]]), A,
                            parallel_distance[distance[i][1]] + x0(self.U[distance[turbine][1]]),
                            perpendicular_distance[distance[i][1]], c1(self.U[distance[turbine][1]]))
                    else:
                        deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0
                else:
                    deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

        # Farm efficiency
        self.profit = 0.0

        for l in range(nt):
            self.profit += power(self.U[l])
        self.efficiency = self.profit * 100.0 / (float(nt) * power(self.U[distance[0][1]]))  # same as using U0
        self.powers = [power(self.U[i]) for i in range(nt)]

        with open('speed_larsen.dat', 'w') as out:
            for i in range(len(self.U)):
                out.write('{0:f}\t{1:f}\n'.format(self.U[i], power(self.U[i]) / power(self.U[distance[0][1]])))

    def larsen_windrose(self):

        self.efficiency = 0.0
        self.profit = []
        self.U = []
        self.powers = []
        self.summation = 0.0
        direction = open('direction_power_larsen.dat', 'w', 1)

        windrose_angle = self.wind_direction
        windrose_speed = self.wind_speed
        windrose_frequency = self.wind_frequency

        layout_y = self.layout_y
        layout_x = self.layout_x

        nt = self.nt  # Number of turbines
        r0 = self.radius  # Turbine rotor radius
        D = 2.0 * r0
        A = pi * r0 ** 2.0
        H = self.hub_height  # Hub height
        ia = 0.08  # Ambient turbulence intensity according to vanluvanee. 8% on average

        def deff(U0):
            return D * sqrt((1.0 + sqrt(1.0 - Ct(U0))) / (2.0 * sqrt(1.0 - Ct(U0))))

        rnb = max(1.08 * D, 1.08 * D + 21.7 * D * (ia - 0.05))
        r95 = 0.5 * (rnb + min(H, rnb))

        def x0(U0):
            return 9.5 * D / ((2.0 * r95 / deff(U0)) ** 3.0 - 1.0)

        def c1(U0):
            return (deff(U0) / 2.0) ** (5.0 / 2.0) * (105.0 / 2.0 / pi) ** (- 1.0 / 2.0) * (Ct(U0) * A * x0(U0)) ** (
                - 5.0 / 6.0)  # Prandtl mixing length

        for wind in range(0, len(windrose_angle)):
            U0 = windrose_speed[wind]  # Free stream wind speed
            angle = windrose_angle[wind]
            angle3 = angle + 180.0
            deficit_matrix = [[0.0 for _ in range(nt)] for _ in range(nt)]
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
                flag = [False for _ in range(nt)]
                proportion = [0.0 for _ in range(nt)]
                perpendicular_distance = [0.0 for _ in range(nt)]
                parallel_distance = [0.0 for _ in range(nt)]
                for i in range(turbine + 1, nt):
                    proportion[distance[i][1]], flag[distance[i][1]], perpendicular_distance[distance[i][1]], \
                    parallel_distance[distance[i][1]] = wake_larsen.determine_if_in_wake_larsen(
                        layout_x[distance[turbine][1]],
                        layout_y[distance[turbine][1]],
                        layout_x[distance[i][1]],
                        layout_y[distance[i][1]], A,
                        c1(U[distance[turbine][1]]),
                        Ct(U[distance[turbine][1]]), angle3,
                        r0, x0(U[distance[turbine][1]]))
                    if parallel_distance[
                        distance[i][1]] > 0.0:  ## Add if proportion is 0, skip operation and deficit_matrix == 0.
                        deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[
                                                                                   distance[i][
                                                                                       1]] * wake_larsen.wake_deficit(
                            U[distance[turbine][1]], Ct(U[distance[turbine][1]]), A,
                            parallel_distance[distance[i][1]] + x0(U[distance[turbine][1]]),
                            perpendicular_distance[distance[i][1]], c1(U[distance[turbine][1]]))
                    else:
                        deficit_matrix[distance[i][1]][distance[turbine][1]] = 0.0

            # Farm efficiency
            profit = 0.0
            efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
            for l in range(nt):
                profit += power(U[l])
            efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]]))
            efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0

            direction.write('{0:f} {1:f}\n'.format(angle, profit))
            self.summation += efficiency_proportion[wind]

        direction.close()


if __name__ == '__main__':
    a = Layout('coordinates.dat')
    # b = Layout('coordinates.dat')
    # c = Layout('coordinates.dat')
    # a.read_windrose('windrose.dat')
    # b.read_windrose('windrose.dat')
    # c.read_windrose('windrose.dat')

    a.radius = 40.0
    # b.radius = 40.0
    # c.radius = 40.0
    # a.hub_height = 70.0
    # a.larsen_windrose()
    # b.jensen_windrose()
    # c.ainslie_windrose()
    # a.larsen_angle(8.5, 180.0)
    # a.ainslie_angle(8.5, 180.0)
    a.jensen_angle(8.5, 180.0)
    print(a.profit)
