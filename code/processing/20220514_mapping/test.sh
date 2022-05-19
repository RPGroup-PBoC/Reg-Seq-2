
echo "test"
awk 'NR%4==1 {
        split($0, a, " ");
    }
    NR%4==2 {
        seq=$0;
    }
    NR%4==3 {
        Dict[a[1]] = seq;
        print Dict[a[1]]
    }
    
    ' test_R1.fastq