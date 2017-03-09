from openmdao.api import ScipyOptimizer, Component, Group, Problem, IndepVarComp, ExecComp
from jensen import jensen_angle


class wake(Component):
    def __init__(self):
        super(wake, self).__init__()
        self.add_param('layout_1', val=0.0)
        self.add_param('layout_2', val=0.0)
        self.add_param('layout_3', val=0.0)
        self.add_output('power', val=0.0)

    def solve_nonlinear(self, params, unknowns, resids):
        unknowns['power'] = jensen_angle([params['layout_1'], params['layout_2'], params['layout_3']], 8.5, 180.0)


if __name__ == '__main__':

    top = Problem()
    root = top.root = Group()

    # Initial value of x and y set in the IndepVarComp.
    root.add('x1', IndepVarComp('x1', 0.0))
    root.add('x2', IndepVarComp('x2', 0.0))
    root.add('x3', IndepVarComp('x3', 0.0))
    root.add('wake', wake())
    root.add('con12', ExecComp('c12 = l2-l1'))
    root.add('con13', ExecComp('c13 = l3-l1'))
    root.add('con23', ExecComp('c23 = l3-l2'))

    root.add('min1', ExecComp('min1 = l1'))
    root.add('min2', ExecComp('min2 = l2'))
    root.add('min3', ExecComp('min3 = l3'))
    root.add('max1', ExecComp('max1 = l1'))
    root.add('max2', ExecComp('max2 = l2'))
    root.add('max3', ExecComp('max3 = l3'))

    root.connect('x1.x1', 'wake.layout_1')
    root.connect('x2.x2', 'wake.layout_2')
    root.connect('x3.x3', 'wake.layout_3')

    root.connect('x1.x1', 'min1.l1')
    root.connect('x2.x2', 'min2.l2')
    root.connect('x3.x3', 'min3.l3')
    root.connect('x1.x1', 'max1.l1')
    root.connect('x2.x2', 'max2.l2')
    root.connect('x3.x3', 'max3.l3')

    root.connect('x1.x1', 'con12.l1')
    root.connect('x2.x2', 'con12.l2')
    root.connect('x1.x1', 'con13.l1')
    root.connect('x3.x3', 'con13.l3')
    root.connect('x2.x2', 'con23.l2')
    root.connect('x3.x3', 'con23.l3')

    top.driver = ScipyOptimizer()
    top.driver.options['optimizer'] = 'COBYLA'
    top.driver.options['maxiter'] = 2000

    top.driver.add_desvar('x1.x1', lower=0.0, upper=1120.0)
    top.driver.add_desvar('x2.x2', lower=40.0, upper=1120.0)
    top.driver.add_desvar('x3.x3', lower=80.0, upper=1120.0)
    top.driver.add_objective('wake.power')
    top.driver.add_constraint('con12.c12', lower=40.0)
    top.driver.add_constraint('con13.c13', lower=40.0)
    top.driver.add_constraint('con23.c23', lower=40.0)

    top.driver.add_constraint('min1.min1', lower=0.0)
    top.driver.add_constraint('min2.min2', lower=40.0)
    top.driver.add_constraint('min3.min3', lower=80.0)
    top.driver.add_constraint('max1.max1', upper=1120.0)
    top.driver.add_constraint('max2.max2', upper=1120.0)
    top.driver.add_constraint('max3.max3', upper=1120.0)

    root.wake.deriv_options['type'] = 'fd'
    # root.wake.deriv_options['form'] = 'central' #'central'

    top.setup()

    top['x1.x1'] = 0.0
    top['x2.x2'] = 40.0
    top['x3.x3'] = 1082.0

    top.run()

    print('\n')
    print('Minimum of %f found at (%f, %f, %f)' % (top['wake.power'], top['wake.layout_1'], top['wake.layout_2'], top['wake.layout_3']))