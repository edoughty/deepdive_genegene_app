#! /bin/bash

psql -U $PGUSER -h $PGHOST -c "TRUNCATE TABLE docids;" $PGDATABASE
