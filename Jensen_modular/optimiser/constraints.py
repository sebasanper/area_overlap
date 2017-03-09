from openmdao.api import Component


class Constraints(Component):
    def __init__(self):
        super(Constraints, self).__init__()
        self.add_param('x', val=[0.0, 0.0, 0.0])
        self.add_output('dx', val=[0.0, 0.0, 0.0])

    def solve_nonlinear(self, params, unknowns, resids):
        x = params['x']
        dx = []
        for i in range(len(x)):
            for j in range(i, len(x)):
                dx.append(x[i] - x[j])
        unknowns['dx'] = dx

if __name__ == '__main__':
    x = [0.0, 56.0, 884.0]
    dx = []
    for i in range(len(x)):
        print i
        for j in range(i+1, len(x)):
            print j
            dx.append(x[j] - x[i])
    print dx