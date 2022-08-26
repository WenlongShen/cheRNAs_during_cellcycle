
#!/usr/bin/awk -f


{
if ($3=="exon")
    {
        split($9, a, " ");
        ID=a[4]; 
        L[ID]+=$5-$4+1
    } 
}
END{
    for(i in L)
        {print i"\t"L[i]}
    }