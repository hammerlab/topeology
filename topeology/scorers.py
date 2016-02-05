from .pmbec import PMBEC

class Scorer(object):
    def __init__(self):
        self.aa = PMBEC()
        self.gap_penalty = abs(self.aa.calculate_min_int())

class CSSWLScorer(Scorer):
    """Use the Complete-Striped-Smith-Waterman-Library library to align seq_a and seq_b."""

    def score_multiple(self, df, col_a, col_b):
        return df.apply(
            lambda row: self.score(row[col_a], row[col_b]),
            axis=1)

    def score(self, seq_a, seq_b):
        # StripedSmithWaterman expects str vs. unicode
        seq_a = str(trim_seq(seq_a))
        seq_b = str(trim_seq(seq_b))

        # From https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/
        # blob/master/README.md:
        # "Note: When SSW open a gap, the gap open penalty alone is applied."
        from skbio.alignment import StripedSmithWaterman
        query = StripedSmithWaterman(
            seq_a, protein=True,
            gap_open_penalty=self.gap_penalty, gap_extend_penalty=self.gap_penalty,
            substitution_matrix=self.aa.as_int_dict())
        return query(seq_b)["optimal_alignment_score"]

class SeqAlignScorer(Scorer):
    """Use the seq-align library, via the pmbecalign C extension, to align seq_a and seq_b."""

    def __init__(self):
        Scorer.__init__(self)

        # Gap penalty = gap_open + gap_extend * length of gap, so we"ll just rely on gap_extend
        from pmbecalign import pmbec_init
        pmbec_init(self.gap_penalty, self.aa.as_int_list())

    def score_multiple(self, df, col_a, col_b):
        new_df = df.copy()
        trimmed_col_a = "trimmed_" + col_a
        trimmed_col_b = "trimmed_" + col_b
        new_df[trimmed_col_a] = new_df[col_a].apply(trim_seq)
        new_df[trimmed_col_b] = new_df[col_b].apply(trim_seq)

        from pmbecalign import pmbec_score_multiple
        return pmbec_score_multiple(list(new_df[trimmed_col_a]),
                                    list(new_df[trimmed_col_b]))

def trim_seq(seq):
    return seq[2:-1]
