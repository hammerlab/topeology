#include "postgres.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fmgr.h"
#include "utils/builtins.h"
#include "string_buffer/string_buffer.h"
#include "smith_waterman.h"
#include <zlib.h>
#include "alignment_scoring_load.h"

#ifdef PG_MODULE_MAGIC
PG_MODULE_MAGIC;
#endif

PG_FUNCTION_INFO_V1(alignment_score);

Datum
alignment_score(PG_FUNCTION_ARGS)
{
  text* seq_a_text = PG_GETARG_TEXT_P(0);
  text* seq_b_text = PG_GETARG_TEXT_P(1);
  char* seq_a = text_to_cstring(seq_a_text);
  char* seq_b = text_to_cstring(seq_b_text);

  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  int match = 1;
  int mismatch = -2;
  int gap_open = -4;
  int gap_extend = -1;
  char no_start_gap_penalty = 1;
  char no_end_gap_penalty = 1;
  char no_gaps_in_a = 0, no_gaps_in_b = 0;
  char no_mismatches = 0;
  char case_sensitive = 0;

  scoring_t scoring;
  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
  gzFile sub_matrix_file = gzopen(PMBEC_FILE_NAME, "r");
  align_scoring_load_matrix(sub_matrix_file, PMBEC_FILE_NAME, &scoring, case_sensitive);
  gzclose(sub_matrix_file);
  scoring.use_match_mismatch = 0;

  smith_waterman_align(seq_a, seq_b, &scoring, sw);

  int max_score = -10000;
  while(smith_waterman_fetch(sw, result)) {
    if(result->score > max_score) {
      max_score = result->score;
    }
  }

  smith_waterman_free(sw);
  alignment_free(result);

  PG_RETURN_INT32(max_score);
}
