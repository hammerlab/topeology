// Copyright (c) 2016. Mount Sinai School of Medicine
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <Python.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include "string_buffer/string_buffer.h"
#include "smith_waterman.h"
#include <zlib.h>
#include "alignment_macros.h"
#include "alignment_scoring.h"
#include "alignment_scoring_load.h"

// Copied from seq-align source
void add_mutations(scoring_t* scoring, const char *str, int *scores,
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
static scoring_t scoring;
static bool scoring_initialized = false;

static char *
string_from_pyobject(PyObject *obj) {
#if PY_MAJOR_VERSION >= 3
  return PyUnicode_AsUTF8(obj);
#else
  return PyBytes_AsString(obj);
#endif
}

// Copied from http://stackoverflow.com/a/13942236
static PyObject *
make_list(double array[], size_t size) {
  PyObject *l = PyList_New(size);
  size_t i;
  for (i = 0; i != size; ++i) {
    PyList_SET_ITEM(l, i, PyFloat_FromDouble(array[i]));
  }
  return l;
}

static bool is_amino_acids(char *seq) {
  size_t len_seq = strlen(seq);
  int i;
  for (i = 0; i < len_seq; i++) {
    // Is this letter is an amino acid?
    if (memchr(AMINO_ACIDS, toupper(seq[i]), sizeof(AMINO_ACIDS)) == NULL) {
      return false;
    }
  }

  return true;
}

int pmbec_score_int(scoring_t *scoring, char *seq_a, char *seq_b) {
  sw_aligner_t *sw = smith_waterman_new();
  alignment_t *result = alignment_create(256);

  smith_waterman_align(seq_a, seq_b, scoring, sw);

  int max_score = INT_MIN;
  while(smith_waterman_fetch(sw, result)) {
    if(result->score > max_score) {
      max_score = result->score;
    }
  }

  smith_waterman_free(sw);
  alignment_free(result);

  return max_score;
}

static PyObject *
pmbec_norm_score(PyObject *self, PyObject *args) {
  if (!scoring_initialized) {
    PyErr_SetString(PyExc_ValueError, "pmbec_init must be called before aligning sequences");
    return NULL;
  }

  char *seq_a, *seq_b;
  if (!PyArg_ParseTuple(args, "ss", &seq_a, &seq_b))
    return NULL;

  if (strlen(seq_a) != strlen(seq_b)) {
    PyErr_SetString(PyExc_ValueError, "Sequence lengths must be equal");
    return NULL;
  }

  if (!is_amino_acids(seq_a) || !is_amino_acids(seq_b)) {
    PyErr_SetString(PyExc_ValueError, "Sequences must only comprise amino acids");
    return NULL;
  }

  int val = pmbec_score_int(&scoring, seq_a, seq_b);

  if (val < 0) {
    return Py_None;
  }

  float val_norm = ((float) val / 100.0) / (float) strlen(seq_a);
  return PyFloat_FromDouble(val_norm);
}

static PyObject *
pmbec_score(PyObject *self, PyObject *args) {
  if (!scoring_initialized) {
    PyErr_SetString(PyExc_ValueError, "pmbec_init must be called before aligning sequences");
    return NULL;
  }

  char *seq_a, *seq_b;
  if (!PyArg_ParseTuple(args, "ss", &seq_a, &seq_b))
    return NULL;

  if (!is_amino_acids(seq_a) || !is_amino_acids(seq_b)) {
    PyErr_SetString(PyExc_ValueError, "Sequences must only comprise amino acids");
    return NULL;
  }

  int val = pmbec_score_int(&scoring, seq_a, seq_b);

  if (val < 0) {
    return Py_None;
  }

  return PyFloat_FromDouble((float) val / 100.0);
}

static PyObject *
pmbec_score_multiple(PyObject *self, PyObject *args) {
  PyObject *array_seq_a, *array_seq_b;
  char **seq_a_values, **seq_b_values;
  double *output;
  int num_seqs_a, num_seqs_b;
  int i;

  if (!PyArg_ParseTuple(args, "O!O!", &PyList_Type, &array_seq_a,
                        &PyList_Type, &array_seq_b)) {
    return NULL;
  }

  num_seqs_a = PyList_Size(array_seq_a);
  num_seqs_b = PyList_Size(array_seq_b);
  if (num_seqs_a != num_seqs_b) {
    PyErr_SetString(PyExc_ValueError, "...");
    return NULL;
  }

  seq_a_values = calloc(num_seqs_a, sizeof(char*));
  seq_b_values = calloc(num_seqs_b, sizeof(char*));
  for (i = 0; i < num_seqs_a; i++) {
    seq_a_values[i] = (char*) string_from_pyobject(PyList_GetItem(array_seq_a, i));
    seq_b_values[i] = (char*) string_from_pyobject(PyList_GetItem(array_seq_b, i));
  }

  output = calloc(num_seqs_b, sizeof(double));
  for (i = 0; i < num_seqs_a; i++) {
    output[i] = pmbec_score_int(&scoring, seq_a_values[i], seq_b_values[i]) / 100.0;
  }

  free(seq_a_values);
  free(seq_b_values);

  PyObject* output_list = make_list(output, num_seqs_a);

  free(output);

  return output_list;
}

static PyObject *
pmbec_init(PyObject *self, PyObject *args) {
  PyObject *array;
  int *pmbec_values;
  int pmbec_count, gap_extend, i;

  if (!PyArg_ParseTuple(args, "iO!", &gap_extend, &PyList_Type, &array)) {
    return NULL;
  }

  // Accept only a positive value for the gap extension
  if (gap_extend < 0) {
    PyErr_SetString(PyExc_ValueError, "Gap extension penalty must be >= 0");
    return NULL;
  }

  pmbec_count = PyList_Size(array);
  if (pmbec_count != 576) {
    PyErr_SetString(PyExc_ValueError, "Substitution matrices must be of length 576");
    return NULL;
  }

  pmbec_values = calloc(pmbec_count, sizeof(int));
  for (i = 0; i < pmbec_count; i++) {
    pmbec_values[i] = (int) PyLong_AsLong(PyList_GetItem(array, i));
  }

  // Default match/mismatch scores; ignorable because matrix covers all scores
  int match = -100;
  int mismatch = -100;

  // Don't ever use the defaults specified above, since the matrix covers all scores
  char use_match_mismatch = 0;

  // Gap penalty = gap_open + gap_extend * length of gap, so we'll just rely on gap_extend
  // without any opening-specific penalties
  gap_extend = gap_extend * -1; // seq-align expects negative gap penalties
  int gap_open = 0;

  // No special treatment for gaps at the start
  char no_start_gap_penalty = 0;
  char no_end_gap_penalty = 0;

  // No forcing of gaplessness in either peptide
  char no_gaps_in_a = 0;
  char no_gaps_in_b = 0;

  // Substitutions are allowed
  char no_mismatches = 0;

  char case_sensitive = 0;

  scoring_init(&scoring, match, mismatch, gap_open, gap_extend,
               no_start_gap_penalty, no_end_gap_penalty,
               no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);

  add_mutations(&scoring, AMINO_ACIDS, pmbec_values, use_match_mismatch);

  scoring_initialized = true;

  free(pmbec_values);

  return Py_None;
}

static PyMethodDef PmbecAlignMethods[] = {
  {"pmbec_init", pmbec_init, METH_VARARGS,
   "Initialize Smith-Waterman alignment scoring using the PMBEC matrix."},
  {"pmbec_score", pmbec_score, METH_VARARGS,
   "Calculate the Smith-Waterman alignment score of two sequences using the PMBEC matrix."},
  {"pmbec_score_multiple", pmbec_score_multiple, METH_VARARGS,
   "Calculate the Smith-Waterman alignment score of two sequence lists using the PMBEC matrix."},
  {"pmbec_norm_score", pmbec_norm_score, METH_VARARGS,
   "Calculate the length-normalized Smith-Waterman alignment score of two sequences using the PMBEC matrix."},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef PmbecAlign = {
  PyModuleDef_HEAD_INIT,
  "pmbecalign",
  "",
  -1,
  PmbecAlignMethods
};
#endif

#if PY_MAJOR_VERSION >= 3
PyObject *
PyInit_pmbecalign(void) {
  return PyModule_Create(&PmbecAlign);
}
#else
void
initpmbecalign(void) {
  Py_InitModule("pmbecalign", PmbecAlignMethods);
}
#endif
