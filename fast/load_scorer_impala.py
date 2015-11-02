import subprocess as sp
from os import path, environ
from shutil import copy2
import platform
from urlparse import urlparse
from impala.dbapi import connect

def handle_false(value):
    # Ensure that false in config isn't interpreted as True
    if not value or value.lower() == 'false':
        return False
    return True

DEBUG = handle_false(environ.get('DEBUG', False))

def run_command(command_str, error_ok=False):
    try:
        sp.check_call(command_str.split(" "))
    except subprocess.CalledProcessError as e:
        if error_ok:
            return
        raise e

def load():
    seq_align_path = environ['SEQ_ALIGN_DIR']
    impala_udf_devel_dir = environ['IMPALA_UDF_DEVEL_DIR']

    file_template = "impala_alignment_score.%s"
    seq_align_libs = path.join(seq_align_path, "libs")
    seq_align_src = path.join(seq_align_path, "src")
    bit_array = path.join(seq_align_libs, "bit_array/bit_array.o")
    string_buffer = path.join(seq_align_libs, "string_buffer/string_buffer.o")

    platform_uname = platform.uname()
    compiler = "g++"
    os_specific_compile_args = ["-emit-llvm", "-O3", "-c", "-fPIC"]
    os_specific_link_args = ["-shared"]

    compile_command = [compiler,
                       file_template % "cc",
                       "-o", file_template % "o",
                       "-I", seq_align_libs,
                       "-I", seq_align_src,
                       "-I", impala_udf_devel_dir]
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

    hdfs_dir = environ['IMPALA_UDF_HDFS_DIR']
    file_on_hdfs = path.join(hdfs_dir, file_template % "so")
    run_command("hdfs dfs -rm -skipTrash %s" % file_on_hdfs, error_ok=True)
    run_command("hdfs dfs -copyFromLocal %s %s" % (file_template % "so", hdfs_dir))

    impala_host = environ['IMPALA_HOST']
    impala_port = environ['IMPALA_PORT']
    impala_db = environ['IMPALA_DB']
    conn = connect(host=impala_host, port=impala_port, database=impala_db)
    cursor = conn.cursor()
    cursor.execute("drop function if exists alignment_score(string, string)")
    cursor.execute(("create function alignment_score(string, string) returns int location "
                    "'%s' "
                    "symbol='alignment_score_impala'") % file_on_hdfs)

load()
