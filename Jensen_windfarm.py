__author__ = 'Sebastian Sanchez Perez Moreno' \
             's.sanchezperezmoreno@tudelft.nl'
# Jensen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
import wake
from math import sqrt, log, tan, cos, floor, ceil
from numpy import deg2rad
import time
from os import makedirs, path, chdir
# powers = ['power7', 'power5', 'power3', 'powertable', 'powerstep']
powers = ['power7']
# thrusts = ['Ct6', 'Ct4', 'Ct3', 'Cttable', 'Ctstep']
thrusts = ['Ct6']

for powertype in powers:  # Loop over power curves
    for thrusttype in thrusts:  # Loop over thrust curves
        print powertype, thrusttype
        for withdata in [True]:#, False]:  # Either measure execution time, or run once completely.
            for rose in ['360']:#, '360']:  # Either run with 30 degrees, or 360 degrees wind roses.

                # newpath = path.join('jensen_results/', powertype, thrusttype, rose)
                newpath = path.join('8mps/', powertype, thrusttype, rose)

                def power(U0):
                    #power7, power5, power3, powertable, powerstep
                    if powertype == 'power7':
                        return power7(U0)

                def Ct(U0):
                    # Ct6, Ct4, Ct3, Cttable, Ctstep
                    if thrusttype == 'Ct6':
                        return Ct6(U0)

                chdir('/home/sebasanper/PycharmProjects/area_overlap')
                layout = open('horns_rev.dat', 'r')
                if rose == '360':
                    windrose = open('horns_rev_windrose2.dat', 'r')
                elif rose == '30':
                    windrose = open('horns_rev_windrose.dat', 'r')

                if not path.exists(newpath):
                    makedirs(newpath)
                chdir(newpath)
                turb_data = open('turb1_power.dat', 'w', 1)
                direction = open('direction_power.dat', 'w', 1)
                eff = open('efficiency.dat', 'w', 1)

                def Ct6(U0):
                    if U0 < 4.0:
                        return 0.1
                    elif U0 <= 25.0:
                        return 7.3139922126945e-7 * U0 ** 6.0 - 6.68905596915255e-5 * U0 ** 5.0 + 2.3937885e-3 * U0 ** 4.0 + - 0.0420283143 * U0 ** 3.0 + 0.3716111285 * U0 ** 2.0 - 1.5686969749 * U0 + 3.2991094727
                    else:
                        return 0.1

                def power7(U0):
                    if U0 < 4.0:
                        return 0.0
                    elif U0 <= 25.0:
                        return 3.234808e-4 * U0 ** 7.0 - 0.0331940121 * U0 ** 6.0 + 1.3883148012 * U0 ** 5.0 - 30.3162345004 * U0 ** 4.0 + 367.6835557011 * U0 ** 3.0 - 2441.6860655008 * U0 ** 2.0 + 8345.6777042343 * U0 - 11352.9366182805
                    else:
                        return 0.0

                layout_x = []
                layout_y = []
                for line in layout:
                    columns = line.split()
                    layout_x.append(float(columns[0]))
                    layout_y.append(float(columns[1]))

                windrose_angle = []
                windrose_speed = []
                windrose_frequency = []
                for line in windrose:
                    columns = line.split()
                    windrose_angle.append(float(columns[0]))
                    windrose_speed.append(float(columns[1]))
                    windrose_frequency.append(float(columns[2]))

                def analysis():
                    nt = len(layout_y)  # Number of turbines ## Length of layout list
                    summation = 0.0

                    def distance_to_front(x, y, theta, r): ## TODO: Calculate distance to any front and use negative distances to order.
                        theta = deg2rad(theta)
                        return abs(x + tan(theta) * y - r / cos(theta)) / sqrt(1.0 + tan(theta) ** 2.0)

                    for wind in range(len(windrose_angle)):
                    # for wind in range(270, 271):
                        # U1 = windrose_speed[wind]  # Free stream wind speed
                        # U0 = U1 * (70.0 / 10.0) ** 0.11 # Power or log law for wind shear profile
                        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
                        U0 = 8.0
                        k = 0.04  # Decay constant
                        r0 = 40.0  # Turbine rotor radius
                        angle = windrose_angle[wind]
                        # angle = - 270.0 - angle  #To read windroses where N is 0 and E is 90.
                        angle3 = angle + 180.0
                        deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
                        proportion = [[0.0 for x in range(nt)] for x in range(nt)]
                        distance = [[0.0 for x in range(2)] for x in range(nt)]

                        U = [U0 for x in range(nt)]
                        total_deficit = [0.0 for x in range(nt)]

                        for tur in range(nt):
                            distance[tur] = [distance_to_front(layout_x[tur], layout_y[tur], angle, 100000000.0), tur]
                        distance.sort()

                        for turbine in range(nt):
                            for num in range(turbine):
                                total_deficit[distance[turbine][1]] += deficit_matrix[distance[turbine][1]][distance[num][1]] ** 2.0
                            total_deficit[distance[turbine][1]] = sqrt(total_deficit[distance[turbine][1]])
                            U[distance[turbine][1]] = U0 * (1.0 - total_deficit[distance[turbine][1]])
                            for i in range(turbine + 1, nt):
                                determ = wake.determine_if_in_wake(layout_x[distance[turbine][1]], layout_y[distance[turbine][1]], layout_x[distance[i][1]], layout_y[distance[i][1]], k, r0, angle3)
                                proportion[distance[turbine][1]][distance[i][1]] = determ[0]
                                # If statement for proportion = 0
                                deficit_matrix[distance[i][1]][distance[turbine][1]] = proportion[distance[turbine][1]][distance[i][1]] * wake.wake_deficit(Ct(U[distance[turbine][1]]), k, determ[1], r0)

                        # for turb in range(nt):
                        #     for i in range(nt):
                        #         output.write('{0:f}\t'.format(deficit_matrix[turb][i]))
                        #     output.write('\n')
                        #     output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(turb, layout_x[turb], layout_y[turb], total_deficit[turb], U[turb], int(angle), power(U[turb])))
                        # output2.write('\n')

                        # for n in range(nt):
                        #     aver[n] += power(U[n]) / 360.0
                        # ---------------------------------------TODO UNCOMMENT -----------------------------------------
                        if withdata:
                            turb_data.write('{0:f} {1:f}\n'.format(angle, power(U[2]))) # Number of turbine you want output of. Now 14.

                        # Farm efficiency
                        profit = 0.0
                        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
                        for l in range(nt):
                            profit += power(U[l])
                        efficiency = profit * 100.0 / (float(nt) * power(U[distance[0][1]])) # same as using U0
                        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
                        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
                    # ---------------------------------------TODO UNCOMMENT ------------------------------------------
                        if withdata:
                            direction.write('{0:f} {1:f}\n'.format(angle, profit))
                        summation += efficiency_proportion[wind]
                    return summation
                    # for n in range(nt):
                    #     turb_data.write('{0:f}\n'.format(aver[n]))
                # turb_data.close()
                #     # output.close()
                #     # output2.close()
                #     # draw.close()
                # direction.close()
                # layout.close()
                # windrose.close()
                    # row.close()
                    # data.close()


                eff.write('{0:f}'.format(analysis()))
                eff.close()
