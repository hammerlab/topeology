export DATABASE_URI='postgres://user@localhost'
export PG_CONFIG='pg_config'
export PG_INCLUDE_DIR_SERVER=`$PG_CONFIG --includedir-server`
export PG_LIB_DIR=`$PG_CONFIG --libdir`
export SEQ_ALIGN_DIR='/path/to/seq-align'
export IMPALA_UDF_DEVEL_DIR='/path/to/udf-devel'
export IMPALA_UDF_HDFS_DIR='/hdfs/path/to/udf-dir'
export IMPALA_HOST='localhost'
export IMPALA_PORT=21050
export IMPALA_DB='db_name'
export DEBUG=False
