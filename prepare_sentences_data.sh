#! /bin/bash

psql -U $PGUSER -h $PGHOST -c "TRUNCATE TABLE sentences;" $PGDATABASE
