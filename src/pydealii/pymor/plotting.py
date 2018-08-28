from pprint import pprint
from pickle import loads
import numpy as np
import matplotlib.pyplot as plt
import sys


def elasticity_error_curves(errors, lambdas):
    for lvl, err in list(errors.items()):
        plt.plot(lambdas, err)
    plt.legend(['lvl {}'.format(g) for g in list(errors.keys())], bbox_to_anchor=(0.3, 0.9),
               bbox_transform=plt.gcf().transFigure)
    plt.title('reference_lvl: {}'.format(9))
    plt.xlabel('Lambda')
    plt.ylabel('|||u_H-u_h|||')
    plt.show()


if __name__ == '__main__':
    elasticity_error_curves(*loads(open(sys.argv[1]).read()))