from openmdao.api import ScipyOptimizer, Component, Group, Problem, IndepVarComp, ExecComp
from jensen import jensen_angle
from numpy import array, zeros, ones, arange, concatenate
from openmdao.devtools.partition_tree_n2 import view_tree

class wake(Component):
    def __init__(self):
        super(wake, self).__init__()
        self.add_param('layout_x', val=zeros(9, ))
        self.add_param('layout_y', val=zeros(9, ))
        self.add_output('power', val=0.0)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['power'] = jensen_angle(params['layout_x'], params['layout_y'], 8.5, 180.0)


if __name__ == '__main__':
    top = Problem()
    root = top.root = Group()

    # Initial value of x and y set in the IndepVarComp.
    root.add('x', IndepVarComp('x', zeros(9, )))
    root.add('y', IndepVarComp('y', zeros(9, )))
    root.add('wake', wake())

    root.connect('x.x', 'wake.layout_x')
    root.connect('y.y', 'wake.layout_y')

    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'SLSQP'
    top.driver.options['maxiter'] = 2000

    top.driver.add_desvar('x.x', lower=zeros(9,), upper=ones(9, ) * 320.0)
    top.driver.add_desvar('y.y', lower=zeros(9,), upper=ones(9, ) * 320.0)
    top.driver.add_objective('wake.power')
    top.driver.options['disp'] = True
    root.wake.deriv_options['type'] = 'fd'
    root.wake.deriv_options['form'] = 'central'

    a = []
    for i in range(9):
        for j in range(i+1, 9):
            a.append('a' + str(i) + str(j))
            root.add(a[-1], ExecComp('r = sqrt((a - b) ** 2.0 + (c - d) ** 2.0)'))

    for i in range(9):
        for j in range(i+1, 9):
            cona = 'a' + str(i) + str(j) + '.a'
            conb = 'a' + str(i) + str(j) + '.b'
            conc = 'a' + str(i) + str(j) + '.c'
            cond = 'a' + str(i) + str(j) + '.d'
            root.connect('x.x', cona, src_indices=[i])
            root.connect('x.x', conb, src_indices=[j])
            root.connect('y.y', conc, src_indices=[i])
            root.connect('y.y', cond, src_indices=[j])

    for element in a:
        const = element + '.r'
        top.driver.add_constraint(const, lower=40.0)
    top.setup()
    # top.check_partial_derivatives()
    # top.check_total_derivatives()
    # view_tree(top)

    b = arange(0.0, 301.0, 150.0)
    top['x.x'] = concatenate((b, b, b))
    r1 = zeros(3,)
    r2 = ones(3, ) * 150.0
    r3 = ones(3,) * 300.0
    top['y.y'] = concatenate((r1, r2, r3))
    # top['y.y'] = zeros(16,)

    top.run()

    print('\n')
    print 'Minimum of %f' % (top['wake.power'])
    with open('results.dat', 'w') as out:
        for i in range(len(top['wake.layout_x'])):
            out.write('%f\t%f\n' % (top['wake.layout_x'][i], top['wake.layout_y'][i]))
