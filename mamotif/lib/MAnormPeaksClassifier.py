import random
from scipy import stats
from sequence import Gene, Peak, PeakSet
from MAmotifPeaks import MAnormPeakSet, MotifScanPeak


class Feature(object):
    """
    feature used for classifying MAnorm peaks
    """

    def is_exist_in_sequence(self, peak):
        pass


class FeatureMotif(Feature):
    """
    using motifscan result as a feature for classify MAnorm peaks
    """

    def __init__(self, motifscan_peak_list):
        assert isinstance(motifscan_peak_list[0], MotifScanPeak)  # motif scan peaks list is a nature of sequence
        self.__motifscan_result = motifscan_peak_list
        self.__group_motifscan_result = {}
        self.__is_grouped = False
        self.mismatch = 0

    def is_exist_in_sequence(self, peak):
        if not self.__is_grouped:
            self.group()

        for pk_motif in self.__group_motifscan_result[peak.chrm]:
            if peak == pk_motif:
                if pk_motif.target_number > 0:
                    return True
                elif pk_motif.target_number == 0:
                    return False
        peak.prints()
        print '@warning: this peak not exist in motifscan result!'
        self.mismatch += 1

    def group(self):
        for motifscan in self.__motifscan_result:
            if motifscan.chrm not in self.__group_motifscan_result.keys():
                self.__group_motifscan_result[motifscan.chrm] = []
            else:
                self.__group_motifscan_result[motifscan.chrm].append(motifscan)
        self.__is_grouped = True


class FeaturePromoter(Feature):
    """
    Using gene's promoter zone to classify peaks into promoter-overlap peaks and else
    """

    def __init__(self, gene_list):
        assert isinstance(gene_list[0], Gene)
        self.__genes = gene_list
        self.__group_genes = {}
        self.__is_grouped = False

    def is_exist_in_sequence(self, peak):
        if not self.__is_grouped:
            self.group()
        try:
            for gene in self.__group_genes[peak.chrm]:
                if peak.is_overlap(gene.promoter_zone):
                    return True
        except KeyError, e:
            print e
            return False
        return False

    def group(self):
        for gene in self.__genes:
            if gene.chrm not in self.__group_genes.keys():
                self.__group_genes[gene.chrm] = []
            else:
                self.__group_genes[gene.chrm].append(gene)
        self.__is_grouped = True


class FeatureOtherPeakOverlap(Feature):
    """
    Using other peaks bed file as a feature to classify peaks by overlap between them
    """

    def __init__(self, pk_list):
        # assert isinstance(pk_list[0], Peak)
        self.__pks = pk_list

    def is_exist_in_sequence(self, peak):
        for pk in self.__pks:
            if peak.is_overlap(pk):
                return True
        return False


class FeatureOtherPeakMatch(Feature):
    def __init__(self, pk_list):
        assert isinstance(pk_list[0], Peak)
        self.__pks = pk_list

    def is_exist_in_sequence(self, peak):
        for pk in self.__pks:
            if peak == pk:
                return True
        return False


class FeatureRandom(Feature):
    """
    Using this feature to classify peaks randomly
    """

    def is_exist_in_sequence(self, peak):
        random_num = random.randint(0, 1)
        if random_num == 0:
            return False
        else:
            return True


class FeatureMvalue(Feature):
    """
    Using mvalue to classify peaks
    """

    def __init__(self, mvalue):
        self.__mvalue = mvalue

    def is_exist_in_sequence(self, peak):
        if peak.mvalue > self.__mvalue:
            return True
        else:
            return False


class FeaturePvalue(Feature):
    """
    Using Pvalue to classify peaks
    """

    def __init__(self, pvalue):
        self.__pvalue = pvalue

    def is_exist_in_sequence(self, peak):
        if peak.pvalue < self.__pvalue:
            return True
        else:
            return False


class FeatureAbsMvalue(Feature):
    """
    Using abs of mvalue to classify peaks
    """

    def __init__(self, abs_mvalue):
        self.__mvalue = abs(abs_mvalue)

    def is_exist_in_sequence(self, peak):
        if abs(peak.mvalue) > self.__mvalue:
            return True
        else:
            return False


class Classifier(object):
    """
    class for classify pks
    """

    def __init__(self, pk_set, feature):
        assert isinstance(pk_set, PeakSet)
        self.__peaks_set = pk_set
        self.__feature = feature
        self._is_classified = False
        self.feature_yes = PeakSet()
        self.feature_no = PeakSet()

    def classify_by_feature(self):
        i = 0
        for pk in self.__peaks_set:
            if self.__feature.is_exist_in_sequence(pk):
                self.feature_yes.append(pk)
            else:
                self.feature_no.append(pk)
            i += 1
            if i % 5000 == 0:
                print '@info: Already classified %d peaks ...' % i
                pass
        self._is_classified = True

    def set_feature(self, peak_set1, peak_set2):
        self.feature_yes = peak_set1
        self.feature_no = peak_set2
        self._is_classified = True


class MAnormPeaksClassifier(Classifier):
    """
    class used for classifying MAnorm peaks into yes and no sets,
    and test result of classified result
    """

    def __init__(self, pk_set, feature):
        assert isinstance(pk_set, MAnormPeakSet)
        Classifier.__init__(self, pk_set, feature)

    def ttest_feature_classified_peaks(self):
        if not self._is_classified:
            self.classify_by_feature()
        mvalue_yes, mvalue_no = \
            [pk.mvalue for pk in self.feature_yes], [pk.mvalue for pk in self.feature_no]
        try:
            t_statistic, two_tailed_pvalue = stats.ttest_ind(mvalue_yes, mvalue_no, equal_var=False)
            if t_statistic < 0:
                pvalue_left = two_tailed_pvalue / 2
                pvalue_right = 1 - two_tailed_pvalue / 2
            else:
                pvalue_left = 1 - two_tailed_pvalue / 2
                pvalue_right = two_tailed_pvalue / 2
            return float(t_statistic), pvalue_right
        except:
            return None

    def ranksum_feature_classified_peaks(self):
        if not self._is_classified:
            self.classify_by_feature()
        mvalue_yes, mvalue_no = [pk.mvalue for pk in self.feature_yes], [pk.mvalue for pk in self.feature_no]
        try:
            z_statistic, two_side_pvalue = stats.ranksums(mvalue_yes, mvalue_no)
            if z_statistic < 0:
                pvalue_left = two_side_pvalue / 2
                pvalue_right = 1 - two_side_pvalue / 2
            else:
                pvalue_left = 1 - two_side_pvalue / 2
                pvalue_right = two_side_pvalue / 2
            return z_statistic, pvalue_right
        except:
            return None

    def kstest_feature_classified_peaks(self):
        if not self._is_classified:
            self.classify_by_feature()
        mvalue_yes, mvalue_no = [pk.mvalue for pk in self.feature_yes], [pk.mvalue for pk in self.feature_no]
        try:
            ks_statistic, two_side_pvalue = stats.ks_2samp(mvalue_yes, mvalue_no)
            if ks_statistic < 0:
                pvalue_left = two_side_pvalue / 2
                pvalue_right = 1 - two_side_pvalue / 2
            else:
                pvalue_left = 1 - two_side_pvalue / 2
                pvalue_right = two_side_pvalue / 2
            return ks_statistic, pvalue_right
        except:
            return None
