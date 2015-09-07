#include "postgres.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "fmgr.h"
#include "utils/builtins.h"
#include "string_buffer/string_buffer.h"
#include "smith_waterman.h"
#include <zlib.h>
#include "alignment_macros.h"
#include "alignment_scoring.h"
#include "alignment_scoring_load.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

// Copied from https://github.com/noporpoise/seq-align/blob/05893831281ac26cfe299777416ef28f29f908f1
// /src/alignment_scoring.c
void add_mutations(scoring_t* scoring, const char *str, const int *scores,
                           char use_match_mismatch) {
  size_t i, j, len = strlen(str);
  char a, b;
  int score;

  for (i = 0; i < len; i++) {
    a = scoring->case_sensitive ? str[i] : tolower(str[i]);
    for (j = 0; j < len; j++) {
      b = scoring->case_sensitive ? str[j] : tolower(str[j]);
      score = ARR_LOOKUP(scores, len, i, j);

      scoring_add_mutation(scoring, a, b, score);
    }
  }

  scoring->use_match_mismatch = use_match_mismatch;
}

static char AMINO_ACIDS[] = "ARNDCQEGHILKMFPSTWYVBZX*";

// Generated from pmbec.py
static const int PMBEC[576] =
  { 100,   -4,  -39,   -7,    6,  -20,   -1,   19,  -32,   -1,  -18,
    -3,  -29,  -18,   22,   17,   24,  -50,  -27,   34, -100, -100,
    -100, -100,   -4,  100,  -14,  -24,  -24,    3,  -31,  -11,   38,
    -34,  -30,   66,  -10,  -20,  -10,    8,  -17,  -20,    3,  -36,
    -100, -100, -100, -100,  -39,  -14,  100,   13,    2,   25,    4,
    -3,   20,  -29,   -7,  -17,    5,  -11,   -5,   11,   -2,   17,
    -1,  -27, -100, -100, -100, -100,   -7,  -24,   13,  100,    8,
    0,   50,   -6,  -10,  -14,  -18,  -31,  -15,   -9,   18,   -7,
    4,    2,  -10,  -12, -100, -100, -100, -100,    6,  -24,    2,
    8,  100,  -23,   13,    8,   -5,    3,  -14,  -28,  -12,   21,
    -5,  -17,  -19,   29,    4,  -12, -100, -100, -100, -100,  -20,
    3,   25,    0,  -23,  100,   18,    3,   -5,  -28,   -1,  -15,
    1,  -43,   19,   21,   17,   -6,  -24,  -19, -100, -100, -100,
    -100,   -1,  -31,    4,   50,   13,   18,  100,   -4,  -20,   -7,
    2,  -37,  -12,   -3,   10,  -12,   -4,   -5,  -13,   -5, -100,
    -100, -100, -100,   19,  -11,   -3,   -6,    8,    3,   -4,  100,
    -16,  -11,  -23,   -8,  -19,   -2,    5,   22,    7,  -16,  -12,
    -9, -100, -100, -100, -100,  -32,   38,   20,  -10,   -5,   -5,
    -20,  -16,  100,  -37,  -19,   32,    0,   -3,  -25,    5,  -14,
    6,   16,  -46, -100, -100, -100, -100,   -1,  -34,  -29,  -14,
    3,  -28,   -7,  -11,  -37,  100,   42,  -14,   26,   23,  -12,
    -29,  -17,    1,   -9,   42, -100, -100, -100, -100,  -18,  -30,
    -7,  -18,  -14,   -1,    2,  -23,  -19,   42,  100,   -8,   47,
    24,  -21,  -36,  -21,   -4,   -7,   27, -100, -100, -100, -100,
    -3,   66,  -17,  -31,  -28,  -15,  -37,   -8,   32,  -14,   -8,
    100,   -4,  -18,  -16,    2,  -17,  -26,   -2,  -18, -100, -100,
    -100, -100,  -29,  -10,    5,  -15,  -12,    1,  -12,  -19,    0,
    26,   47,   -4,  100,   16,  -39,  -19,  -30,   12,    7,   -7,
    -100, -100, -100, -100,  -18,  -20,  -11,   -9,   21,  -43,   -3,
    -2,   -3,   23,   24,  -18,   16,  100,  -30,  -42,  -40,   28,
    35,   -7, -100, -100, -100, -100,   22,  -10,   -5,   18,   -5,
    19,   10,    5,  -25,  -12,  -21,  -16,  -39,  -30,  100,   -3,
    9,  -20,  -31,    9, -100, -100, -100, -100,   17,    8,   11,
    -7,  -17,   21,  -12,   22,    5,  -29,  -36,    2,  -19,  -42,
    -3,  100,   55,  -25,  -29,   -5, -100, -100, -100, -100,   24,
    -17,   -2,    4,  -19,   17,   -4,    7,  -14,  -17,  -21,  -17,
    -30,  -40,    9,   55,  100,  -16,  -25,   31, -100, -100, -100,
    -100,  -50,  -20,   17,    2,   29,   -6,   -5,  -16,    6,    1,
    -4,  -26,   12,   28,  -20,  -25,  -16,  100,   37,  -16, -100,
    -100, -100, -100,  -27,    3,   -1,  -10,    4,  -24,  -13,  -12,
    16,   -9,   -7,   -2,    7,   35,  -31,  -29,  -25,   37,  100,
    -22, -100, -100, -100, -100,   34,  -36,  -27,  -12,  -12,  -19,
    -5,   -9,  -46,   42,   27,  -18,   -7,   -7,    9,   -5,   31,
    -16,  -22,  100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100, -100, -100, -100, -100, -100, -100, -100,
    -100, -100, -100, -100 };

PG_FUNCTION_INFO_V1(alignment_score);

Datum
alignment_score(PG_FUNCTION_ARGS)
{
  text *seq_a_text = PG_GETARG_TEXT_P(0);
  text *seq_b_text = PG_GETARG_TEXT_P(1);
  char *seq_a = text_to_cstring(seq_a_text);
  char *seq_b = text_to_cstring(seq_b_text);

  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  // Default match/mismatch scores; ignorable because matrix covers all scores
  int match = -100;
  int mismatch = -100;

  // Don't ever use the defaults specified above, since the matrix covers all scores
  char use_match_mismatch = 0;

  // Gap penalty is equal to min(PMBEC); see pmbec_min() in pmbec.py
  int gap_open = -50;
  int gap_extend = -50;

  // No special treatment for gaps at the start
  char no_start_gap_penalty = 0;
  char no_end_gap_penalty = 0;

  // No forcing of gaplessness in either peptide
  char no_gaps_in_a = 0;
  char no_gaps_in_b = 0;

  // Substitutions are allowed
  char no_mismatches = 0;

  char case_sensitive = 0;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);

  add_mutations(&scoring, AMINO_ACIDS, PMBEC, use_match_mismatch);

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  int max_score = INT_MIN;
  while(smith_waterman_fetch(sw, result)) {
    if(result->score > max_score) {
      max_score = result->score;
    }
  }

  smith_waterman_free(sw);
  alignment_free(result);

  PG_RETURN_INT32(max_score);
}
