__author__ = 'sebasanper'

__author__ = 'sebasanper'
# Jensen wake model with partial shadowing factor applied to horns rev. Must change Ct according to wind speed.
from math import sqrt, log
import time
import wake_ainslie
# 80 turbines
output = open('matrix.dat', 'w')
output.write('# This file has the wake deficit matrix per turbine per wind direction\n')
output2 = open('final_speed.dat', 'w')
output2.write('# This file has the deficit, wind speed and power at each turbine per wind direction.\n# Turbine number\tX-coordinate\tY-coordinate\tTotal speed deficit\tTotal wind speed\tWind direction angle\tPower produced\n')
layout = open('horns_rev.dat', 'r')
draw = open('draw_horns_rev.dat', 'w')
draw.write('# This file has the turbines affected by the wake of one turbine at one direction.\n')
# draw2 = open('drawline.dat', 'w')
turb_data = open('turb_data.dat', 'w')
direction = open('direction_efficiency.dat', 'w')
direction.write('# This file includes the efficiency of the whole farm by wind direction.\n# Wind direction angle\tFarm efficiency\n')


def Jensen(a, Nt):
    windrose = open('horns_rev_windrose2.dat', 'r')
    nt = Nt
    layout_x = [0.0 for x in range(nt)]
    layout_y = [0.0 for x in range(nt)]
    for x in range(nt):
        layout_x[x] = float(a[x][0])
        layout_y[x] = float(a[x][1])

    windrose_angle = []
    windrose_speed = []
    windrose_frequency = []
    for line in windrose:
        columns = line.split()
        windrose_angle.append(float(columns[0]))
        windrose_speed.append(float(columns[1]))
        windrose_frequency.append(float(columns[2]))

    windrose.close()
    summation = 0.0

    def Ct(U0):
        return 0.0001923077 * U0 ** 4.0 + -0.0075407925 * U0 ** 3.0 + 0.096462704 * U0 ** 2.0 - 0.5012354312 * U0 + 1.7184749184

    def power(U0):
        return -0.0071427272 * U0 ** 5.0 + 0.5981106302 * U0 ** 4.0 - 18.5218157059 * U0 ** 3.0 + 251.0929636046 * U0 ** 2.0 - 1257.8070377904 * U0 + 2043.2240149783

    # def Cp(U0):
    #     return power(U0) / (0.5 * 1.225 * pi * r0**2.0 * U0**3.0)

    # for U0 in range(4, 20):
    #     summation = 0.0
    for wind in range(0, len(windrose_speed)):
    # for wind in range(10, 11):
    #     U1 = windrose_speed[wind]  # Free stream wind speed
        # U0 = U1 * (70.0 / 10.0)**0.11 # Power or log law for wind shear profile
        # U0 = U1 * log(70.0 / 0.005) / log(10.0 / 0.005)
        U0 = 8.5
        k = 0.04  # Decay constant
        r0 = 40.0  # Turbine rotor radius
        angle = windrose_angle[wind]
        angle3 = angle + 180.0
        wake_deficit_matrix = [[0.0 for x in range(nt)] for x in range(nt)]
        for turbine in range(nt):
            flag = [False for x in range(nt)]
            proportion = [0.0 for x in range(nt)]
            for i in range(nt):
                proportion[i], flag[i] = wake_ainslie.determine_if_in_wake(layout_x[turbine], layout_y[turbine], layout_x[i], layout_y[i], k, r0, angle3)
            #     if angle == 300 and turbine == 75:
            #         if flag[i] is True and i != turbine:
            #             draw.write('{0:d}\t{1:.1f}\t{2:.1f}\t1\n'.format(i, layout_x[i], layout_y[i]))
            #         elif flag[i] is not True and i != turbine:
            #             draw.write('{0:d}\t{1:.1f}\t{2:.1f}\t0\n'.format(i, layout_x[i], layout_y[i]))
            #         elif i == turbine:
            #             draw.write('{0:d}\t{1:.1f}\t{2:.1f}\t2\n'.format(i, layout_x[i], layout_y[i]))
            # draw.write('\n')

            # Matrix with effect of each turbine <i = turbine> on every other turbine <j> of the farm
            for j in range(nt):
                if turbine != j and flag[j] is True:
                    wake_deficit_matrix[turbine][j] = proportion[j] * wake_ainslie.wake_deficit(0.81, k, wake_ainslie.distance(layout_x[turbine], layout_y[turbine], layout_x[j], layout_y[j]), r0)
                elif turbine == j:
                    wake_deficit_matrix[turbine][j] = 0.0
            #     output.write('{0:f}\t'.format(wake_deficit_matrix[j][turbine]))
            # output.write('\n')

        total_deficit = [0.0 for x in range(nt)]
        total_speed = [0.0 for x in range(nt)]
        for j in range(nt):
            for i in range(nt):
                total_deficit[j] += wake_deficit_matrix[i][j] ** 2.0
            total_deficit[j] = sqrt(total_deficit[j])
            total_speed[j] = U0 * (1.0 - total_deficit[j])
        turb_data.write('{0:f}\t{1:f}\n'.format(angle, power(total_speed[14])/923.0))
        #     output2.write('{0:d}\t{1:.1f}\t{2:.1f}\t{3:f}\t{4:f}\t{5:d}\t{6:f}\n'.format(j, layout_x[j], layout_y[j], total_deficit[j], total_speed[j], int(angle), power(total_speed[j])))
        # output2.write('\n')

        # Individual turbine efficiency
        # if angle == 180:
        #     for g in range(0, 80):
        #             turb_data.write('{0:d}\t{1:f}\n'.format(g, power(total_speed[g]) / power(total_speed[g])))

        # Farm efficiency
        profit = 0.0
        efficiency_proportion = [0.0 for x in range(0, len(windrose_frequency))]
        efficiency = 0.0
        for l in range(nt):
            profit += power(total_speed[l])
            efficiency = profit * 100.0 / (80.0 * power(U0))
        efficiency_proportion[wind] = efficiency * windrose_frequency[wind] / 100.0
        # print 'Farm efficiency with wind direction = {0:d} deg: {1:2.2f}%'.format(int(angle), efficiency)
        # print angle, efficiency
        # direction.write('{0:f}\t{1:f}\n'.format(angle, efficiency))
        summation += efficiency_proportion[wind]
    # print 'total farm efficiency is {0:f} %'.format(summation)
    #     print U0, summation
    return summation
#
turb_data.close()
output.close()
output2.close()
draw.close()
direction.close()

if __name__ == '__main__':
    start_time = time.time()
    Jensen()
print("--- %s seconds ---" % (time.time() - start_time))