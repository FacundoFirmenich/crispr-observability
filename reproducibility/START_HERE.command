#!/bin/sh
set -eu
cd -- "$(dirname -- "$0")"
exec ./crispr-reproduce --menu
