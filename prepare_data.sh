#! /usr/bin/env bash

export PGUSER=${PGUSER:-`whoami`}
export PGPASSWORD=${PGPASSWORD:-}
export PGDATABASE=${PGDATABASE}

dropdb -U $PGUSER -h $PGHOST -p $PGPORT $PGDATABASE
createdb -U $PGUSER -h $PGHOST -p $PGPORT $PGDATABASE

psql -U $PGUSER -h $PGHOST -p $PGPORT -c "drop schema if exists public cascade; create schema public;" $PGDATABASE


psql -U $PGUSER -h $PGHOST -p $PGPORT -c "CREATE TABLE genegene(
	docid text,
	mid1 text,
	mid2 text,
	mention1 text,
	mention2 text,
	is_correct boolean,
    features text[],
	sentence text,
	id bigint);" $PGDATABASE

psql -U $PGUSER -h $PGHOST -p $PGPORT -c "CREATE TABLE docids(
	id bigint,
	docid text,
	folder text);" $PGDATABASE

psql -U $PGUSER -h $PGHOST -p $PGPORT -c "CREATE TABLE documents(
	id bigint,
	docid text,
	document text);" $PGDATABASE

psql -U $PGUSER -h $PGHOST -p $PGPORT -c "CREATE TABLE sentences(
	id bigint,
	docid text,
	sentid text,
	sentence text,
	text text);" $PGDATABASE

