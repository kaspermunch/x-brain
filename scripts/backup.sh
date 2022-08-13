#!/bin/bash
 
# The source path to backup. Can be local or remote.
#SOURCE=$HOME/ChrXh2/people/kmt/x-brain
SOURCE=x-brain

# Where to store the incremental backups
DESTBASE=`dirname $SOURCE`/BACKUP
 
# Where to store today's backup
DEST="$DESTBASE/$(date +%Y-%m-%d)"
# Where to find yesterday's backup
YESTERDAY="$DESTBASE/$(date -d yesterday +%Y-%m-%d)/"
 
# Use yesterday's backup as the incremental base if it exists
if [ -d "$YESTERDAY" ]
then
	OPTS="--link-dest $YESTERDAY"
fi
 
# Run the rsync
rsync -av $OPTS "$SOURCE" "$DEST"
