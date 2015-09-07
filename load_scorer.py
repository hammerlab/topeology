import subprocess as sp
from os import path, environ
from shutil import copy2
import psycopg2
import platform
from urlparse import urlparse

def handle_false(value):
    # Ensure that false in config isn't interpreted as True
    if not value or value.lower() == 'false':
        return False
    return True

DEBUG = handle_false(environ.get('DEBUG', False))

def pg_conn(db_uri):
    db_uri_parsed = urlparse(db_uri)
    username = db_uri_parsed.username
    password = db_uri_parsed.password
    database = db_uri_parsed.path[1:]
    hostname = db_uri_parsed.hostname

    # Weirdly, need to do this to avoid requiring a password. Handle this better.
    if hostname == 'localhost':
        hostname = None

    return psycopg2.connect(
        database=database,
        user=username,
        host=hostname,
        password=password)

def load():
    seq_align_path = environ['SEQ_ALIGN_DIR']
    pg_include_dir_server = environ['PG_INCLUDE_DIR_SERVER']
    pg_lib_dir = environ['PG_LIB_DIR']

    file_template = "alignment_score.%s"
    seq_align_libs = path.join(seq_align_path, "libs")
    seq_align_src = path.join(seq_align_path, "src")
    bit_array = path.join(seq_align_libs, "bit_array/bit_array.o")
    string_buffer = path.join(seq_align_libs, "string_buffer/string_buffer.o")

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
                       "-D PMBEC_FILE_NAME=\"%s\"" % environ['PMBEC_FILE_NAME'],
                       "-o", file_template % "o",
                       "-Wall", "-Wextra",
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
                    string_buffer,
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

    # Using copy2 requires sudo-ing the whole file, which interferes with DB roles,
    # so for now I'm just using cp.
    sp.check_call(["sudo",
                   "cp",
                   file_template % "so",
                   pg_lib_dir])
    pg_lib_dir_file_path = path.join(pg_lib_dir, file_template % "so")

    conn = pg_conn(environ['DATABASE_URI'])
    cur = conn.cursor()
    cur.execute("LOAD '%s'" % pg_lib_dir_file_path)
    cur.execute("DROP FUNCTION IF EXISTS alignment_score(text, text)")
    create_fn_sql = ("CREATE FUNCTION alignment_score(text, text) RETURNS integer "
                     "AS '%s', 'alignment_score' LANGUAGE C STRICT;") % pg_lib_dir_file_path
    cur.execute(create_fn_sql)
    conn.commit()

load()
