# GERP++
* [Home Page](http://mendel.stanford.edu/sidowlab/downloads/gerp/index.html)


The code from GERP++ Home Page is not available, but I found a copy in Github [tvkent/GERPplusplus](https://github.com/tvkent/GERPplusplus).
You can also find it [here](gerp++KRT).

After downloading it, unzip it and run `make`. **(Please note that only available in Linux, can't in Windows)**

How to run?

`-t` tree file
`-f` fasta file
```shell
gerpcol -t NP_000006.2.blast.dnd -f NP_000006.2.blast.aln.fa -e NP_000006.2
```

### Update 2022.10.28
I found the code from github can't run.

### Update 2022.11.14
Correct the code from github, and run successfully. But it only can calculate scores for `ATCG`.

1. Error Source
    
    When reading MSA in fasta sequences, it will read one more special char (maybe `\n`),
    which makes it can't map the MSA to tree file.
2. How to solve?
   1. Update 71th line in `gerp__KRT` from `return _title;` to `return _title.substr(0, _title.length() - 1);`
        ![](https://xdcat-tuchuang.oss-cn-beijing.aliyuncs.com/img/20221114215142.png)
   2. Run `make`
   3. Use it!!

3. Example
    1. The origin code from github `gerp++KRT`
   2. The modified code `gerp++KRT_patch`
      * gerpcol executor `gerpcol`
      * MSA Fasta input `NP_000006.2.blast.aln.fa`
      * Tree input `NP_000006.2.blast.dnd`
      * Result `NP_000006.2.blast.aln.fa.rates`
4. **Note**
* My input examples come from ClustalOmega, and fasta file is generate by biopython.
* **Make tree file in one line without space string.**
* It will only calculate the results for `ATCG`, can't apply to protein seq.


