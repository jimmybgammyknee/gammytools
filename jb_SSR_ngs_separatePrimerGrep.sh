#!/bin/bash -l

# Jimmy (jimmymbreen@gmail.com)	-	20141002
# Take Merged FastQ file and split by F/R primer and create histogram of seq lengths

if [ "$#" != "2" ]; then
        echo "Usage: jb_SSR_ngs_separatePrimerGrep.sh [Merged fastQ.gz] [primer list]"
        exit 0
fi

## Primer list must be in the format:
## SampleName	Forward_seq	Reverse_seq	Est_Amplicon_Size
## (separared by tabs)

# File IO and naming
file=$1
primers=$2
fname=$(basename $file .fastq.gz)

gammytools="/Users/jbreen/programs/gammytools"
fqCon=$gammytools/jb_fastq2fasta.py
plotH=$gammytools/jb_plot_histogram.py

# Create fasta from fastq and Remove newline in fastafile
zcat $file | python $fqCon | awk '!/^>/ { printf "%s", $0; n = "\n" }  /^>/ { print n $0; n = "" } END { printf "%s", n }' > "$fname"_new.fa

while read line
 do
        # Pipe in variables
        Fseq=$(echo $line | awk '{print $2}')
        Rseq_rev=$(echo $line | awk '{print $3}' | tr ACGT TGCA | rev)
        Rseq=$(echo $line | awk '{print $3}')
        Fseq_rev=$(echo $line | awk '{print $2}' | tr ACGT TGCA | rev)
        name=$(echo $line | awk '{print $1}')
#       regex=$(echo $line | awk '{print $4}')

        # Look for 100% matches in Forward direction and remove primer sequence
        grep "^"$Fseq".*"$Rseq_rev"$" "$fname"_new.fa | \
                sed -e "s/^$Fseq//g" -e "s/$Rseq_rev$//g"  >> "$fname"_"$name".seq.list

        # Give a summary of the number of files found
        echo -e "\n"$name" Forward Reads" >> "$fname"_summary_list.txt
        wc -l "$fname"_"$name".seq.list >> "$fname"_summary_list.txt

        # Look for 100% matches in Reverse direction and remove primer sequence
        grep "^"$Rseq".*"$Fseq_rev"$" "$fname"_new.fa | \
                sed -e "s/^$Rseq//g" -e "s/$Fseq_rev$//g" | tr ACGT TGCA | rev  >> "$fname"_"$name".seq.list

        # Give a summary of the number of files found
        echo -e "\nTotal Sequences for "$name" (including Reverse Reads)" >> "$fname"_summary_list.txt
        wc -l "$fname"_"$name".seq.list >> "$fname"_summary_list.txt

       # Remove sequences that dont contain the SSR
#       awk -v ref="$regex" 'match($0, ref) {print $0}' "$fname"_"$name".seq.list > "$fname"_"$name".seqflt.list
#       echo -e "\nTotal Sequences for "$name" that actually contain the SSR" >> "$fname"_summary_list.txt
#       wc -l "$fname"_"$name".seqflt.list >> "$fname"_summary_list.txt

        # Print fasta and align using clustal-omega
        #awk '/^/{print ">seq"(++i)}!/>/' "$fname"_"$name".seqflt.list > "$fname"_"$name".fasta
        awk '/^/{print ">seq"(++i)}!/>/' "$fname"_"$name".seq.list > "$fname"_"$name".fasta
	clustal-omega-1.2.0-macosx -i "$fname"_"$name".fasta --output-order=tree-order | \
                awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > "$fname"_"$name".aligned.fa

        # Make length histogram
        #awk '{print length}' "$fname"_"$name".seqflt.list | sort | \
	awk '{print length}' "$fname"_"$name".seq.list | sort | \
                uniq -c | sort -k2,2n -r | awk '{print $2"\t"$1}' > "$fname"_"$name".seq.list.histogram

        # Plot histogram with separate code
        python $plotH "$fname"_"$name".seq.list.histogram > "$fname"_"$name".seq.list.histogram.pdf

        # Clean up after yourself
        rm "$fname"_"$name".seq.list
        rm "$fname"_"$name".seq.list.histogram
        #rm "$fname"_"$name".seqflt.list
done < $primers

rm "$fname"_new.fa
