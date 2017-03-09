from openmdao.api import ScipyOptimizer, Component, Group, Problem, IndepVarComp, ExecComp
from jensen import jensen_angle
from numpy import array, zeros, ones, arange


class wake(Component):
    def __init__(self):
        super(wake, self).__init__()
        self.add_param('layout', val=zeros(9,))
        self.add_output('power', val=0.0)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['power'] = jensen_angle(params['layout'], zeros(9,), 8.5, 180.0)

if __name__ == '__main__':

    top = Problem()
    root = top.root = Group()

    # Initial value of x and y set in the IndepVarComp.
    root.add('x', IndepVarComp('x', zeros(9, )))
    root.add('wake', wake())

    root.connect('x.x', 'wake.layout')

    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'SLSQP'
    top.driver.options['maxiter'] = 2000

    top.driver.add_desvar('x.x', lower=zeros(9,), upper=ones(9,) * 4480.0)
    top.driver.add_objective('wake.power')
    root.wake.deriv_options['type'] = 'fd'
    root.wake.deriv_options['form'] = 'central'

    a = []
    for i in range(9):
        for j in range(i + 1, 9):
            a.append('a' + str(i) + str(j))
            root.add(a[-1], ExecComp('r = abs(a - b)'))

    for i in range(9):
        for j in range(i + 1, 9):
            cona = 'a' + str(i) + str(j) + '.a'
            conb = 'a' + str(i) + str(j) + '.b'
            root.connect('x.x', cona, src_indices=[i])
            root.connect('x.x', conb, src_indices=[j])

    for element in a:
        const = element + '.r'
        top.driver.add_constraint(const, lower=40.0)

    top.setup()

    top['x.x'] = arange(0.0, 4100.0, 500.0)

    top.run()

    print('\n')
    print 'Minimum of %f' % (top['wake.power'])
    with open('1Dresults.dat', 'w') as out:
        for i in range(len(top['wake.layout'])):
            out.write('%f\t0.0\n' % (top['wake.layout'][i]))
