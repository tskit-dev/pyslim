#!/bin/bash
# CC-BY https://answers.launchpad.net/inkscape/+question/48164
# $1 is the file to extract layers from
# $2: are the layers to always include
#
# needs xmlstarlet
#
# Notes:
#  - doesn't deal with spaces in object names

if ! command -v "xmlstarlet" &>/dev/null || ! command -v 'inkscape' &>/dev/null;
then
    echo >&2 "can't find xmlstarlet and/or inkscape, not regenerating SVG-derived figures"
    exit 0
fi

TMPFILE=$(mktemp /tmp/output.XXXXXXXXX).svg
cp $1 $TMPFILE

YESLAYERS=${*:2}
ALLLAYERS=$(xmlstarlet sel -t -m "//*[@inkscape:groupmode=\"layer\"]" -v "concat(@inkscape:label,' ')" $TMPFILE|sort -Vr)
NOLAYERS=$(comm -1 -3 <(echo $YESLAYERS|tr ' ' '\n'|sort) <(echo $ALLLAYERS|tr ' ' '\n'|sort))

for layer in $NOLAYERS
do
    id=$(xmlstarlet sel -t -m "//*[@inkscape:label=\"$layer\"]" -v "@id" $TMPFILE)
    xmlstarlet ed -S -L -d "//*[@id=\"$id\"]" $TMPFILE
done

# --export-area-drawing will crop to visible
inkscape --export-area-drawing --export-type="svg" --export-filename=- $TMPFILE

[ -f $TMPFILE ] && rm $TMPFILE
