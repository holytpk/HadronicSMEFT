''' 
Multi-dimensional polonomial parametrization.
'''

# Logger
import logging
logger = logging.getLogger(__name__)

# General imports
import operator
import numpy as np
import scipy.special
import itertools
from functools import reduce

class HyperPoly:

    @staticmethod
    def get_ndof( nvar, order ):
        ''' Compute the number of d.o.f. of the polynomial by summing up o in the formula for combinations with repetitions of order o in nvar variables'''
        return sum( [ int(scipy.special.binom(nvar + o - 1, o)) for o in xrange(order+1) ] )

    def __init__( self, order ):
        self.order       = order
        self.initialized = False

    # Initialize with data ( c1, ..., cN )
    def initialize( self, base_points, ref_point = None):

        # Let's not allow re-initialization.
        if self.initialized:
            raise RuntimeError( "Already initialized!" )

        # Make sure dimensionality of data is consistent
        if not len(set( map( len, base_points ) )) == 1:
            raise ValueError( "'base_points' are not consistent. Need a list of iterables of the same size." )

        # Length of the dataset
        self.N   = len( base_points )
        # Coordinates
        self.base_points = base_points 
        # Number of variables
        self.nvar = len( base_points[0] ) 
        # Total number of DOF
        self.ndof = HyperPoly.get_ndof( self.nvar, self.order )

        logger.debug( "Make parametrization of polynomial in %i variables to order %i" % (self.nvar, self.order ) )

        # sanity
        assert self.N == self.ndof, "Got %i base_points but was expecting %i for %i variables at order %i."%( self.N, self.ndof, self.nvar, self.order )

        # Reference point
        self.ref_point = ref_point if ref_point is not None else tuple([0 for var in xrange(self.nvar)])

        # Order of combinations (with replacements) and ascending in 'order'
        # Fill Theta, the |A|x|A| matrix in w^k = Theta^k_A w^A where w^k are the base-point weights and w^A the coefficients   
        # We expand around the reference point
        self.combination  = {}
        counter = 0
        self.Theta = np.empty( [self.ndof, self.ndof ], dtype='float64')
        for o in xrange(self.order+1):
            for comb in itertools.combinations_with_replacement( xrange(self.nvar), o ):
                self.combination[counter] = comb
                self.combination[comb]    = counter
                for k in range(self.ndof):
                    #self.Theta[k][counter] = reduce(operator.mul, [base_points[k][c]-self.ref_point[c] for c in comb], 1)
                    self.Theta[k][counter] = reduce(operator.mul, [base_points[k][c] for c in comb], 1)
                counter += 1

        self.ThetaInv = np.linalg.inv(self.Theta)

        ## Optionally expand around zero by multiplying with a matrix that translates the coefficients
        #self.expand_around_zero = np.zeros( (self.ndof,self.ndof))
        #for i_comb in range(self.ndof):
        #    comb_i = self.combination[i_comb]
        #    for o in xrange(self.order+1):
        #        for comb_j in itertools.combinations( comb_i, o ):
        #            const = reduce(operator.mul, [-self.ref_point[c] for c in comb_j], 1)
        #            remainder = list(comb_i)
        #            for c in comb_j:
        #                remainder.remove(c)
        #            remainder = tuple(remainder)                            
        #            #print ("from:", i_comb, "const comps:", comb_j, "const:", const, "contributes to:", remainder)
        #            self.expand_around_zero[self.combination[remainder],i_comb]+=const 

        # Check reference point
        if len(self.ref_point)!=self.nvar:
            logger.error('Reference point has length %i but should have length %i', len(self.ref_point), self.nvar )
            raise RuntimeError

    #def get_parametrization( self, weights, expand_around_zero=True): 
    def get_parametrization( self, weights): 
        ''' Obtain the parametrization for given weights
        '''

        #if expand_around_zero:
        #    return np.dot(np.dot(self.expand_around_zero, self.ThetaInv), np.array(weights))
        #else:
        return np.dot(self.ThetaInv, np.array(weights))
         
