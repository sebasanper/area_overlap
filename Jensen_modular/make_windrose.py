with open('windrose360.dat', 'w') as out:
    for i in range(360):
        out.write('{0:d}\t{1:f}\t{2:f}\n'.format(i, 8.0, 0.2777777))
