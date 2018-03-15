"""Options for calculating uncorrected p-values."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang et al., All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
import sys

class PvalCalcBase(object):
    """Base class for initial p-value calculations."""

    def __init__(self, name, pval_fnc, log, test_type):
        self.log = log
        self.name = name
        self.pval_fnc = pval_fnc
        self.test_type = test_type

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """pvalues are calculated in derived classes."""
        fnc_call = "calc_pvalue({SCNT}, {STOT}, {PCNT} {PTOT})".format(
            SCNT=study_count, STOT=study_n, PCNT=pop_count, PTOT=pop_n)
        raise Exception("NOT IMPLEMENTED: {FNC_CALL} using {FNC}.".format(
            FNC_CALL=fnc_call, FNC=self.pval_fnc))


class FisherClass(PvalCalcBase):
    """From the 'fisher' package, use function, pvalue_population."""

    options = cx.OrderedDict([
        ('up', 'right_tail'),
        ('down', 'left_tail'),
        ('updown', 'two_tail')
    ])

    def __init__(self, name, log, test_type):
        import fisher
        super(FisherClass, self).__init__(name, fisher.pvalue_population, log, test_type)

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected p-values."""
        # k, n = study_true, study_tot,
        # K, N = population_true, population_tot
        # def pvalue_population(int k, int n, int K, int N): ...
        pval = self.pval_fnc(study_count, study_n, pop_count, pop_n)
        return getattr(pval, self.options[self.test_type])


class FisherScipyStats(PvalCalcBase):
    """From the scipy stats package, use function, fisher_exact."""

    options = cx.OrderedDict([
        ('up', 'greater'),
        ('down', 'less'),
        ('updown', 'two-sided')
    ])

    def __init__(self, name, log, test_type):
        from scipy import stats
        super(FisherScipyStats, self).__init__(name, stats.fisher_exact, log, test_type)

    def calc_pvalue(self, study_count, study_n, pop_count, pop_n):
        """Calculate uncorrected p-values."""
        # http://docs.scipy.org/doc/scipy-0.17.0/reference/generated/scipy.stats.fisher_exact.html
        #
        #         Atlantic  Indian                              YES       NO
        # whales     8        2    | 10 whales    study_genes    8 scnt   2    | 10 = study_n
        # sharks     1        5    |  6 sharks    not s_genes    1        5    |  6
        #         --------  ------                            --------   -----
        #            9        7      16 = pop_n     pop_genes    9 pcnt   7      16 = pop_n
        #
        # We use the preceeding table to find the p-value for whales/sharks:
        #
        # >>> import scipy.stats as stats
        # >>> oddsratio, pvalue = stats.fisher_exact([[8, 2], [1, 5]])
        #                                              a  b    c  d
        a = study_count
        b = study_n - study_count
        c = pop_count - study_count
        d = pop_n - pop_count - b
        _, p_uncorrected = self.pval_fnc( [[a, b], [c, d]], alternative=self.options[self.test_type])
        return p_uncorrected


class FisherFactory(object):
    """Factory for choosing a fisher function."""

    options = cx.OrderedDict([
        ('fisher', FisherClass),
        ('fisher_scipy_stats', FisherScipyStats),
    ])

    def __init__(self, **kws):
        self.log = kws['log'] if 'log' in kws else sys.stdout
        self.pval_fnc_name = kws["pvalcalc"] if "pvalcalc" in kws else "fisher"
        self.test_type = kws['test_type'] if 'test_type' in kws else 'updown'
        self.pval_obj = self._init_pval_obj()

    def _init_pval_obj(self):
        """Returns a Fisher object based on user-input."""
        if self.pval_fnc_name in self.options.keys():
            try:
                fisher_obj = self.options[self.pval_fnc_name](self.pval_fnc_name, self.log, self.test_type)
            except ImportError:
                print("fisher module not installed.  Falling back on scipy.stats.fisher_exact")
                fisher_obj = self.options['fisher_scipy_stats']('fisher_scipy_stats', self.log, self.test_type)

            return fisher_obj

        raise Exception("PVALUE FUNCTION({FNC}) NOT FOUND".format(FNC=self.pval_fnc_name))

    def __str__(self):
        return " ".join(self.options.keys())


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang et al., All rights reserved.
