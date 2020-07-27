#!/usr/local/bin/gawk -f
{if ( ($1 =="AMRMultiGrid::") )
{if ( ($2 =="iteration") )
print gensub(/,/,"",1,$4) " " gensub(/,/,"",1,$8)}}