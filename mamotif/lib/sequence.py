import numpy as np


class Sequence(object):
    """
    top-level class of all second general sequence, including reads sequence, peaks sequence, genes sequence etc.
    """
    def __init__(self, chrm, start, end, strand=None):
        """
        one sequence should have information of chromosome, start, end, and strand.
        @type chrm: str
        @type end: int
        @type start: int
        @param chrm: which chromosome that the sequence belong to.
        @param start: the start point on ref genome of the sequence
        @param end: the end point on ref genome of the sequence
        @param strand: strand(+ or -) the sequence belong to.
        """

        self.__chrm = chrm.lower()
        # assert isinstance(start, int)
        self.__start = start
        # assert isinstance(end, int)
        self.__end = end
        if strand in ('+', '-'):
            self.__strand = strand
        elif strand is None:
            self.__strand = None

    @property
    def chrm(self):
        """
        @return: sequence chromosome
        """
        return self.__chrm

    @property
    def start(self):
        """
        @rtype : int
        @return: sequence start point
        """
        return self.__start

    @property
    def end(self):
        """
        @rtype : int
        @return: end point of sequence
        """
        return self.__end

    @property
    def strand(self):
        """
        @return: strand (+ or -) of sequence
        """
        return self.__strand

    @property
    def length(self):
        return self.__end - self.__start

    @property
    def midpoint(self):
        """
        @return: midpoint of sequence
        """
        return (self.__end + self.__start) / 2

    def prints(self):
        """
        print sequence
        """
        if self.__strand is None:
            print '\t'.join([self.chrm, str(self.start), str(self.end)])
        else:
            print '\t'.join([self.chrm, str(self.start), str(self.end), self.strand])

    def __eq__(self, other):
        """
        redefine ==
        """
        return self.chrm == other.chrm and \
               self.start == other.start and \
               self.end == other.end and \
               self.strand == other.strand

    def __gt__(self, other):
        """
        redefine >
        """
        return self.chrm == other.chrm and self.start > other.start


class Gene(Sequence):
    """
    gene has more information, including tss etc.
    """
    def __init__(self, name, name2, chrm, start, end, strand):
        Sequence.__init__(self, chrm, start, end, strand)
        self.__name = name
        self.name2 = name2  # gene official symbol

    @property
    def tss(self):
        """
        get transcription start sites
        @return: transcription start sites
        """
        if self.strand == '+':
            return self.start
        elif self.strand == '-':
            return self.end
        else:
            return None

    @property
    def promoter_zone(self, upstream=2000, downstream=2000):
        if self.strand == '+':
            promoter_start = self.tss - upstream
            promoter_end = self.tss + downstream
            return Sequence(self.chrm, promoter_start, promoter_end, strand=self.strand)
        elif self.strand == '-':
            promoter_start = self.tss - downstream
            promoter_end = self.tss + upstream
            return Sequence(self.chrm, promoter_start, promoter_end, strand=self.strand)

    @property
    def name(self):
        """
        name of gene
        """
        return self.__name

    def prints(self):
        if self.strand is None:
            print '\t'.join([self.name, self.__chrm, str(self.__start), str(self.__end)])
        else:
            print '\t'.join([self.name, self.__chrm, str(self.__start), str(self.__end), self.strand])


class Peak(Sequence):
    """
    peak is some kind of binding site of genome, sometimes we need to know overlap between peaks
    """
    def __init__(self, chrm, start, end, summit=None):  # the summit is relative position of start
        Sequence.__init__(self, chrm, start, end)
        if summit is not None:
            self.__summit = start + summit
        else:
            self.__summit = summit

    @property
    def summit(self):
        if self.__summit is None:
            return self.midpoint  # when summit is not known, take the midpoint as summit.
        return self.__summit

    def is_overlap(self, other_seq):
        """
        whether overlap with another peak
        @param other_seq: another peak
        """
        if self.chrm == other_seq.chrm:
            if self.end < other_seq.start or other_seq.end < self.start:
                return False
            else:
                return True
        else:
            return False

    def dist_with_gene(self, gene):
        """
        calculate the distance between peak summit and gene tss
        @param gene: Gene instance
        """
        if self.chrm == gene.chrm:
            return abs(self.summit - gene.tss)
        else:
            print '@warning: %s is not the same chromosome with this peak!' % gene.name

    def prints(self):
        if self.summit is not None:
            return '%s\t%d\t%d\t%d' % (self.chrm, self.start, self.end, self.summit - self.start)
        else:
            return '%s\t%d\t%d' % (self.chrm, self.start, self.end)

    def tostring(self):
        if self.summit is not None:
            return '%s\t%d\t%d\t%d\n' % (self.chrm, self.start, self.end, self.summit - self.start)
        else:
            return '%s\t%d\t%d\n' % (self.chrm, self.start, self.end)

    def __eq__(self, other):
        return self.chrm.lower() == other.chrm.lower() and \
               self.start == other.start and \
               self.end == other.end


class Sequences(object):
    """
    this class define a lot basic operations of sequences list
    """
    def __init__(self, seq_list=None):
        self.__sequences = []
        self.__size = 0
        if seq_list is not None:
            self.set_sequences(seq_list)

    def set_sequences(self, sequence_list):
        """
        set a list of sequences for operation
        @type sequence_list: list
        @param sequence_list: list of sequences
        """
        self.__sequences = sequence_list
        self.__size = len(self.__sequences)

    def append(self, a_seq):
        """
        add one sequence
        """
        self.__sequences.append(a_seq)
        self.__size += 1

    def __add__(self, other):
        seq_list = self.__sequences + other.sequences
        seqs = Sequences()
        seqs.set_sequences(seq_list)
        return seqs

    @property
    def chromosomes_set(self):
        """
        chromosomes group
        @return: list of chromosomes group
        """
        chroms = []
        for seq in self.__sequences:
            if seq.chrm not in chroms:
                chroms.append(seq.chrm)
        return chroms

    @property
    def sequences(self):
        return self.__sequences

    def sort_sequences(self):
        """
        sort peaks by start
        """
        grps = self.group_by_chromosomes()
        sorted_sequences = []
        chrms = grps.keys()
        num_chrms = []
        alpha_chrms = []
        for chrm in chrms:
            m = chrm.replace('chr', '')
            try:
                m = int(m)
                num_chrms.append(m)
            except:
                alpha_chrms.append(m)
        num_chrms = np.array(num_chrms)
        num_chrms.sort()
        alpha_chrms = np.array(alpha_chrms)
        alpha_chrms.sort()
        sorted_chrms = \
            ['chr' + str(num) for num in num_chrms.tolist()] + ['chr' + alpha for alpha in alpha_chrms.tolist()]
        for chrm in sorted_chrms:
            chrm_seqs = grps[chrm]
            chrm_idx = np.array([seq.start for seq in chrm_seqs]).argsort()
            sorted_sequences += [chrm_seqs[i] for i in chrm_idx]
        self.set_sequences(sorted_sequences)

    def group_by_chromosomes(self):
        """
        group sequence by chromosomes
        @return: dict (keys are chromosomes name)
        """
        groups = {}
        for chrm in self.chromosomes_set:
            groups[chrm] = []
        for sq in self.__sequences:
            groups[sq.chrm].append(sq)
        return groups

    def group_by_strand(self):
        """
        @return: dict(keys are - and +)
        """
        if self.__sequences[0].strand is not None:
            grp = {'-': [], '+': []}
            for seq in self.__sequences:
                grp[seq.strand].append(seq)
            return grp
        else:
            print 'the sequence do not have information of strand!'

    @property
    def size(self):
        """
        @return: size of sequences list
        """
        return self.__size

    def intersect(self, other):
        intersect_seqlist = []
        for seq1 in self.sequences:
            for seq2 in other.sequences:
                if seq1 == seq2:
                    intersect_seqlist.append(seq1)
                    break
        inter_seqs = Sequences()
        inter_seqs.set_sequences(intersect_seqlist)
        return inter_seqs

    def __iter__(self):
        return iter(self.sequences)

    def __getitem__(self, item):
        if isinstance(item, str):
            return self.group_by_chromosomes()[item]
        elif isinstance(item, int):
            return self.sequences[item]


class GeneSet(Sequences):
    """
    Using for dealing with gene set
    """

    def __init__(self, seq_list=None):
        Sequences.__init__(self, seq_list)
        self.__group = None
        self.__is_grouped = False
        self.__group_tss = None

    def set_sequences(self, sequence_list):
        assert isinstance(sequence_list[0], Gene)  # using first element to tell gene list
        Sequences.set_sequences(self, sequence_list)

    def get_gene_by_name(self, name):
        for gene in self.sequences:
            if name.lower() == gene.name.lower():
                gene.prints()
                return gene
        print 'Can not find %s' % name

    def find_target_gene(self, peak, is_nearest=False):
        if not self.__is_grouped:
            self.__group = self.group_by_chromosomes()
            self.__is_grouped = True
            self.__group_tss = \
                {key: np.array([gene.tss for gene in self.__group[key]]) for key in self.__group.keys()}
        try:
            genes_chr = self.__group[peak.chrm]
        except:
            print 'This peak do not have target gene in this gene set!'
            return []

        if is_nearest:
            distances = abs(self.__group_tss[peak.chrm] - peak.summit)
            idx_min = distances.argmin()
            # return [genes_chr[idx_min]]
            if distances[idx_min] < 100000:
                return [genes_chr[idx_min]]
            else:
                return []

        return [gene for gene in genes_chr if peak.is_overlap(gene.promoter_zone)]

    def get_all_tss(self):
        return np.array([gene.tss for gene in self.sequences])


class PeakSet(Sequences):
    """
    peaks set
    """

    def __init__(self, seq_list=None):
        Sequences.__init__(self, seq_list)

    def set_sequences(self, pks_list):
        # assert isinstance(pks_list[0], Peak)
        Sequences.set_sequences(self, pks_list)
