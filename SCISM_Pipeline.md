\#Log in and then go to where our files are

    cd /home/shared/inbre-group1/Protist_metagenomes

Now we need to count how many raw reads there are in every gzip file
====================================================================

    for i in *.gz
    do
    zcat $1 | wc -l
    done

The above output for each file needs to be divided over 4 to get the actual number of reads in each .gz raw file
================================================================================================================

Results are in a google sheet: <a href="https://docs.google.com/spreadsheets/d/1prtZ__sc0rpe5pEABET-1FkiZi-zJ5vN156N1KvVYmQ/edit#gid=0" class="uri">https://docs.google.com/spreadsheets/d/1prtZ__sc0rpe5pEABET-1FkiZi-zJ5vN156N1KvVYmQ/edit#gid=0</a>
======================================================================================================================================================================================================================================================

Let\`s see our table
====================

    cat Read_data.csv

    ## file_name,raw_file_count_lines,raw_reads
    ## Korr_all_R1.fastq.gz,169031356,42257839
    ## Korr_all_R2.fastq.gz,169031356,42257839
    ## Para_neg_R1.fastq.gz,317970048,79492512
    ## Para_neg_R2.fastq.gz,317970048,79492512
    ## Para_R1.fastq.gz,188380316,47095079
    ## Para_R2.fastq.gz,188380316,47095079
    ## Sip_R1.fastq.gz,110517900,27629475
    ## Sip_R2.fastq.gz,110517900,27629475
    ## SK_all_R1.fastq.gz,222086924,55521731
    ## SK_all_R2.fastq.gz,222086924,55521731

Read\_Data table
================

<table>
<thead>
<tr class="header">
<th>file_name</th>
<th>raw_file_count_lines</th>
<th>raw_reads</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Korr_all_R1.fastq.gz</td>
<td>169031356</td>
<td>42257839</td>
</tr>
<tr class="even">
<td>Korr_all_R2.fastq.gz</td>
<td>169031356</td>
<td>42257839</td>
</tr>
<tr class="odd">
<td>Para_neg_R1.fastq.gz</td>
<td>317970048</td>
<td>79492512</td>
</tr>
<tr class="even">
<td>Para_neg_R2.fastq.gz</td>
<td>317970048</td>
<td>79492512</td>
</tr>
<tr class="odd">
<td>Para_R1.fastq.gz</td>
<td>188380316</td>
<td>47095079</td>
</tr>
<tr class="even">
<td>Para_R2.fastq.gz</td>
<td>188380316</td>
<td>47095079</td>
</tr>
<tr class="odd">
<td>Sip_R1.fastq.gz</td>
<td>110517900</td>
<td>27629475</td>
</tr>
<tr class="even">
<td>Sip_R2.fastq.gz</td>
<td>110517900</td>
<td>27629475</td>
</tr>
<tr class="odd">
<td>SK_all_R1.fastq.gz</td>
<td>222086924</td>
<td>55521731</td>
</tr>
<tr class="even">
<td>SK_all_R2.fastq.gz</td>
<td>222086924</td>
<td>55521731</td>
</tr>
</tbody>
</table>
