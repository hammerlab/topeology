import subprocess as sp
from os import path, environ
from shutil import copy2
import psycopg2
import platform

def handle_false(value):
    # Ensure that false in config isn't interpreted as True
    if not value or value.lower() == 'false':
        return False
    return True

DEBUG = handle_false(environ.get('DEBUG', False))

def load():
    seq_align_path = environ['SEQ_ALIGN_DIR']
    pg_include_dir_server = environ['PG_INCLUDE_DIR_SERVER']
    pg_lib_dir = environ['PG_LIB_DIR']

    file_template = "alignment_score.%s"
    seq_align_libs = path.join(seq_align_path, "libs")
    seq_align_src = path.join(seq_align_path, "src")
    bit_array = path.join(seq_align_libs, "bit_array/bit_array.o")

    platform_uname = platform.uname()
    compiler = "gcc"
    os_specific_compile_args = ["-fpic"]
    os_specific_link_args = ["-shared"]
    if "Darwin" in platform_uname:
        compiler = "cc"
        os_specific_link_args = ["-bundle", "-flat_namespace", "-undefined", "suppress"]
        os_specific_compile_args = None

    compile_command = [compiler,
                       file_template % "c",
                       "-o", file_template % "o",
                       "-Wall", "-Wextra",
                       "-fpic",
                       "-c",
                       "-I", seq_align_libs,
                       "-I", pg_include_dir_server,
                       "-I", seq_align_src]
    if os_specific_compile_args:
        compile_command.extend(os_specific_compile_args)
    if DEBUG:
        compile_command.append("-g")
        print(sp.list2cmdline(compile_command))
    sp.check_call(compile_command)

    link_command = [compiler,
                    file_template % "o",
                    "-o", file_template % "so",
                    bit_array,
                    "-I", pg_include_dir_server,
                    "-I", seq_align_src,
                    "-L", seq_align_libs,
                    "-L", seq_align_src,
                    "-lalign", "-lpthread", "-lz"]
    if os_specific_link_args:
        link_command.extend(os_specific_link_args)
    if DEBUG:
        link_command.append("-g")
        print(sp.list2cmdline(link_command))
    sp.check_call(link_command)

    copy2(file_template % "so", pg_lib_dir)
    pg_lib_dir_file_path = path.join(pg_lib_dir, file_template % "so")

    conn = psycopg2.connect(environ['DATABASE_URI'])
    cur = conn.cursor()
    cur.execute("LOAD '%s'" % pg_lib_dir_file_path)
    try:
        cur.execute("DROP FUNCTION alignment_score(text, text)")
    except:
        pass
    create_fn_sql = ("CREATE FUNCTION alignment_score(text, text) RETURNS integer "
                     "AS '%s', 'alignment_score' LANGUAGE C STRICT;") % pg_lib_dir_file_path
    cur.execute(create_fn_sql)
    conn.commit()

load()
