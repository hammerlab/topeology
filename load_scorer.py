import subprocess as sp
from os import path, environ
from shutil import copy2
import psycopg2

def load():
    seq_align_path = environ['SEQ_ALIGN_DIR']
    pg_include_dir_server = environ['PG_INCLUDE_DIR_SERVER']
    pg_lib_dir = environ['PG_LIB_DIR']

    file_template = "alignment_score.%s"
    seq_align_libs = path.join(seq_align_path, "libs")
    seq_align_src = path.join(seq_align_path, "src")
    bit_array = path.join(seq_align_libs, "bit_array/bit_array.o")

    compile_command = ["cc", "-g", file_template % "c",
                       "-o", file_template % "o",
                       "-Wall", "-Wextra",
                       "-fpic", "-c",
                       "-I", seq_align_libs,
                       "-I", pg_include_dir_server,
                       "-I", seq_align_src]
    print(sp.list2cmdline(compile_command))
    sp.check_output(compile_command)

    link_command = ["cc", "-g", "-bundle", "-flat_namespace", "-undefined", "suppress",
                    "-o", file_template % "so",
                    file_template % "o", bit_array,
                    "-I", pg_include_dir_server,
                    "-I", seq_align_src,
                    "-L", seq_align_libs,
                    "-L", seq_align_src,
                    "-lalign", "-lpthread", "-lz"]
    sp.check_output(link_command)

    copy2(file_template % "so", pg_lib_dir)
    pg_lib_dir_file_path = path.join(pg_lib_dir, file_template % "so")

    conn = psycopg2.connect("dbname='tavi' user='tavi' host='localhost'")
    cur = conn.cursor()
    cur.execute("LOAD '%s'" % pg_lib_dir_file_path)
    fun = ("CREATE FUNCTION alignment_score(text, text) RETURNS integer "
           "AS '%s', 'alignment_score' LANGUAGE C STRICT;") % pg_lib_dir_file_path
    cur.execute(fun)
    conn.commit()

load()
