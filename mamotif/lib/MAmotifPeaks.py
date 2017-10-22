import numpy as np
from sequence import Peak, Sequences, PeakSet


class MAnormPeak(Peak):
    """
    class for dealing with MAnorm Peak
    """
    def __init__(self, chrm, start, end, summit=None):
        Peak.__init__(self, chrm, start, end, summit)
        self.__mvalue = None  # MAnorm_Mvalue
        self.__avalue = None  # MAnorm_Avalue
        self.__pvalue = None  # MAnorm_Pvalue

    def set_mvalue(self, mvalue):
        self.__mvalue = mvalue

    def set_avalue(self, avalue):
        self.__avalue = avalue

    def set_pvalue(self, pvalue):
        self.__pvalue = pvalue

    @property
    def mvalue(self):
        if self.__mvalue is not None:
            return self.__mvalue
        else:
            print '@warning: the MAnorm mvalue is not exist!'

    @property
    def avalue(self):
        if self.__avalue is not None:
            return self.__avalue
        else:
            print '@warning: the MAnorm avalue is not exist!'

    @property
    def pvalue(self):
        if self.__pvalue is not None:
            return self.__pvalue
        else:
            print '@warning: the MAnorm pvalue is not exist!'

    def is_in_gene(self, gene):
        """
        tell a MAnormPeak is or not in a gene's promoter region
        @param gene: a instance of Gene
        """
        if self.chrm == gene.chrm:
            if self.start > gene.end or self.end < gene.start:
                return False
            else:
                return True
        else:
            return False

    def prints(self):
        """
        print sequence
        """
        print '\t'.join(
            [self.chrm, str(self.start), str(self.end), str(self.mvalue), str(self.avalue), str(self.pvalue)]
        )

    def tostring(self):
        if self.summit is not None:
            return '%s\t%d\t%d\t%d\t%f\t%f\t%s\n' % \
                   (self.chrm, self.start, self.end, self.summit - self.start,
                    self.mvalue, self.avalue, str(self.pvalue))
        else:
            return '%s\t%d\t%d\t%f\t%f\t%s\n' % (self.chrm, self.start, self.end,
                                                 self.mvalue, self.avalue, str(self.pvalue))


class MAnormPeakSet(PeakSet):
    """
    class for dealing with MAnorm peaks set
    """
    def __init__(self, seq_list=None):
        Sequences.__init__(self, seq_list)

    def set_sequences(self, peaks_list):
        # assert isinstance(peaks_list[0], MAnormPeak)  # using first element of peaks list to tell list's element type
        Sequences.set_sequences(self, peaks_list)

    def sort_peaks(self, with_val='mvalue'):
        """
        sort MAnorm peaks with mvalue or pvalue
        @param with_val: mvalue or pvalue
        """
        if with_val == 'mvalue':
            sort_indices = np.argsort(np.array([peak.mvalue for peak in self.sequences]))
            self.set_sequences(np.array(self.sequences)[sort_indices].tolist())  # sort peaks with new order
        elif with_val == 'pvalue':
            sort_indices = np.argsort(np.array([peak.pvale for peak in self.__sequences]))
            self.set_sequences(np.array(self.__sequences)[sort_indices].tolist())

    def mvalues_info(self):
        import numpy as np
        mvalues = np.array([pk.mvalue for pk in self.sequences])
        mean, std = mvalues.mean(), mvalues.std()
        return mean, std

    def modify_outlier(self):
        up_limit = self.mean + 3 * self.std
        down_limit = self.mean - 3 * self.std
        for pk in self.sequences:
            if pk.mvalue > up_limit:
                pk.set_mvalue(up_limit)
            elif pk.mvalue < down_limit:
                pk.set_mvalue(down_limit)
        pass

    def find_target_pks(self, gene, limit_pvalue=0.01):
        target_pks = []
        for pk in self.sequences:
            if pk.is_overlap(gene.promoter_zone) and pk.pvalue < limit_pvalue:
                target_pks.append(pk)
        return target_pks

    @property
    def mean(self):
        import numpy as np
        return np.array([pk.mvalue for pk in self.sequences]).mean()

    @property
    def std(self):
        import numpy as np
        return np.array([pk.mvalue for pk in self.sequences]).std()


class MotifScanPeak(Peak):
    """
    class used for dealing with MotifScan result
    """
    def __init__(self, chrm, start, end, summit=None):
        Peak.__init__(self, chrm, start, end)
        self.__motif_name = None
        self.__target_number = None

    def set_motif_info(self, motif_name, tarnum):
        self.__motif_name = motif_name
        self.__target_number = tarnum

    @property
    def target_number(self):
        if self.__target_number is not None:
            return self.__target_number
        else:
            print '@error: the target number is not exist!'

    def prints(self):
        """
        print sequence
        """
        if self.target_number is not None and isinstance(self.target_number, int):
            print '\t'.join([self.chrm, str(self.start), str(self.end), str(self.target_number)])
        else:
            Peak.prints(self)
