export DATABASE_URI='postgres://user@localhost'
export PG_CONFIG='pg_config'
export PG_INCLUDE_DIR_SERVER=`$PG_CONFIG --includedir-server`
export PG_LIB_DIR=`$PG_CONFIG --libdir`
export SEQ_ALIGN_DIR='/path/to/seq-align'
export DEBUG=False

