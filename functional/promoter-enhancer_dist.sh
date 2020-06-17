#!/bin/bash
# $1 = gene_list.txt
# $2 = enhancer.txt
# $3 = output.txt
# make sure format is refseq and isoform identifiers .1 etc have been removed

dos2unix $1 &>/dev/null
REF_EXON=~/genome/mm10_bt2/refgene-limited.bed
list=""
firstline=1
while read gene; do
    exons="$(grep $gene'\s' $REF_EXON)"
    [[ $exons == '' ]] && continue
    first="$(echo "$exons" | head -n1)"
    last="$(echo "$exons" | tail -n1)"

    first=($first)
    last=($last)
    strand=${first[5]}
    
    if [ "$strand" == "+"  ]
    then	
	tss=${first[1]}
	line=${first[@]}
    else
	tss=${last[2]}
	line=${last[@]}
    fi
    line=($line)
    addline=${line[0]}$'\t'$tss$'\t'$((tss+1))$'\t'${line[@]:3:5}
    if [ $firstline -eq 0 ]; then addline=$'\n'$addline; fi
    firstline=0
    list="$list$addline"
done < $1

sorted="$(echo "$list" | sort -k1,1 -k2,2n -)"
closest="$(echo "$sorted" | bedtools closest -a - -b $2)"

hitlist=""
firstline=1
while read hit; do
    hit=($hit)
    distance=$((${hit[1]}-${hit[7]}/2-${hit[8]}/2))
    addhit=${hit[3]}$'\t'${distance#-}
    if [ $firstline -eq 0 ]; then addhit=$'\n'$addhit; fi
    firstline=0
    hitlist="$hitlist$addhit"
done <<< "$closest"

echo "$hitlist" > $3
